#!/usr/bin/env python3
"""Build validated native references, then resume the full-cell SpaGRN notebook."""
import copy,fcntl,json,os,subprocess,sys,tempfile,time,traceback
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'.T577_spaGRN_work'
RES=ROOT/'GRN_resource/macaque'
PYTHON='/data/lvyongji/envs/conda_envs/spagrn/bin/python'
JUPYTER_ROOT=Path('/tmp/T577_native_workflow_jupyter')

def now():return datetime.now(timezone.utc).isoformat(timespec='seconds')
def atomic_json(path,value):
    path=Path(path)
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=str(path.parent))
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(value,f,ensure_ascii=False,indent=1);f.write('\n');f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    except BaseException:
        if os.path.exists(tmp):os.unlink(tmp)
        raise

def status(stage,**details):
    row=dict(stage=stage,updated_utc=now(),pid=os.getpid(),**details)
    atomic_json(WORK/'native_workflow_status.json',row)
    print(json.dumps(row),flush=True)

def execute_notebook(path,kernel,preserve_completed=False):
    import nbformat
    from nbclient import NotebookClient
    spec=JUPYTER_ROOT/'kernels'/kernel
    spec.mkdir(parents=True,exist_ok=True)
    atomic_json(spec/'kernel.json',dict(argv=[PYTHON,'-m','ipykernel_launcher','-f','{connection_file}'],display_name='Python ('+kernel+')',language='python'))
    os.environ['JUPYTER_PATH']=str(JUPYTER_ROOT)
    book=nbformat.read(str(path),as_version=4)
    nbformat.validate(book)
    class CheckpointClient(NotebookClient):
        def process_message(self,msg,cell,cell_index):
            result=super().process_message(msg,cell,cell_index)
            if time.monotonic()-self.last_checkpoint>60:
                atomic_json(path,self.nb)
                self.last_checkpoint=time.monotonic()
            return result
    client=CheckpointClient(book,kernel_name=kernel,timeout=None,startup_timeout=300,resources={'metadata':{'path':str(ROOT)}},allow_errors=False)
    client.last_checkpoint=time.monotonic()
    try:
        with client.setup_kernel(cwd=str(ROOT)):
            for index,cell in enumerate(list(book.cells)):
                if cell.cell_type!='code':continue
                restore=preserve_completed and index<=9 and cell.get('execution_count') is not None
                saved=copy.deepcopy(cell) if restore else None
                runtime=copy.deepcopy(cell) if restore else cell
                status('running_notebook_cell',notebook=path.name,cell_index=index,preserve_prior_output=restore)
                try:
                    client.execute_cell(runtime,index)
                finally:
                    if restore:book.cells[index]=saved
                    atomic_json(path,book)
                print('CELL_COMPLETE notebook={} index={} time={}'.format(path.name,index,now()),flush=True)
    except BaseException:
        atomic_json(path,book)
        raise
    if preserve_completed:
        book.metadata['spagrn_analysis'].update(state='complete',native_rankings_ready=True,completed_outputs_preserved_during_kernel_rehydration=True,remaining_stages=[])
    else:
        book.metadata['native_reference'].update(rankings_ready=True,state='complete')
    nbformat.validate(book)
    atomic_json(path,book)
    status('notebook_complete',notebook=path.name)

def select_gpu():
    # Recheck live GPU availability after the hours-long reference build.
    for attempt in range(120):
        result=subprocess.run(['nvidia-smi','--query-gpu=index,name,memory.free,memory.used','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True)
        rows=[]
        for line in result.stdout.splitlines():
            index,name,free,used=[s.strip() for s in line.split(',')]
            rows.append(dict(index=int(index),name=name,memory_free_mib=int(free),memory_used_mib=int(used)))
        eligible=[r for r in rows if r['memory_free_mib']>=12288]
        if eligible:
            chosen=max(eligible,key=lambda r:(r['memory_free_mib'],-r['index']))
            record=dict(selected_gpu=chosen,selection_utc=now(),minimum_free_mib=12288,all_gpus=rows)
            atomic_json(WORK/'gpu_selection.json',record)
            os.environ['CUDA_VISIBLE_DEVICES']=str(chosen['index'])
            return record
        status('waiting_for_available_gpu',attempt=attempt+1,maximum_attempts=120)
        time.sleep(30)
    raise RuntimeError('No GPU with at least 12 GiB free within the 1-hour selection budget; no existing GPU job was interrupted')

def main():
    WORK.mkdir(exist_ok=True)
    lock=open('/tmp/T577_native_workflow.lock','a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',NUMBA_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONUNBUFFERED='1',MPLCONFIGDIR='/tmp/T577_spagrn_mpl',IPYTHONDIR='/tmp/T577_spagrn_ipython',JUPYTER_RUNTIME_DIR='/tmp/T577_native_workflow_runtime',NUMBA_CACHE_DIR='/tmp/T577_spagrn_numba')
    for key in ['MPLCONFIGDIR','IPYTHONDIR','JUPYTER_RUNTIME_DIR','NUMBA_CACHE_DIR']:Path(os.environ[key]).mkdir(parents=True,exist_ok=True)
    atomic_json(WORK/'native_workflow_process.json',dict(pid=os.getpid(),started_utc=now(),runner=str(Path(__file__).resolve()),python=sys.executable))
    try:
        manifest=json.loads((RES/'reference_manifest.json').read_text())
        if manifest.get('native_rankings_ready') is True:
            status('verifying_completed_native_reference')
            import hashlib
            def checksum(path):
                digest=hashlib.sha256()
                with open(path,'rb') as handle:
                    for chunk in iter(lambda:handle.read(8<<20),b''):digest.update(chunk)
                return digest.hexdigest()
            records=manifest['native_reference_validation']
            assert {r['path'] for r in records}==set(manifest['ranking_files']) and len(records)==2
            for record in records:
                assert record['gene_count']==28592 and record['motif_count']==10249
                assert checksum(record['path'])==record['sha256'], 'Completed native ranking checksum changed'
            annotation=json.loads((RES/'build/tf_annotation_preparation.json').read_text())
            for pathkey,hashkey in [('motif_annotation_file','motif_annotation_sha256'),('tf_list_file','tf_list_sha256')]:
                assert manifest[pathkey]==annotation[pathkey]
                assert checksum(manifest[pathkey])==annotation[hashkey]
            status('completed_native_reference_reused',ranking_files=manifest['ranking_files'])
        else:
            status('building_native_macaque_reference')
            execute_notebook(ROOT/'T577_SpaGRN_macaque_reference.ipynb','spagrn-t577-reference')
            manifest=json.loads((RES/'reference_manifest.json').read_text())
            assert manifest['native_rankings_ready'] is True
        status('selecting_live_gpu')
        gpu=select_gpu()
        status('resuming_T577_SpaGRN',gpu=gpu['selected_gpu'])
        execute_notebook(ROOT/'T577_SpaGRN.ipynb','spagrn-t577-native',preserve_completed=True)
        summary=json.loads((ROOT/'T577_SpaGRN.summary.json').read_text())
        assert summary['status']=='complete' and summary['all_cells_preserved'] is True
        assert summary['n_cells']==288801 and summary['n_genes']==15541
        status('complete',summary=summary)
        return 0
    except BaseException as exc:
        status('failed',error_type=type(exc).__name__,error=str(exc))
        traceback.print_exc()
        return 1

if __name__=='__main__':sys.exit(main())
