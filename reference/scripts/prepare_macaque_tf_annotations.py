"""Transfer available public HGNC motif annotations using Ensembl 91 orthology.
This T577-specific TF list uses unambiguous native loci and one-to-one orthologs.
Full native gene background and all public motifs are retained in the rankings.
"""
from pathlib import Path
import csv,io,zipfile,json,hashlib,os
from collections import Counter

ROOT=Path(__file__).resolve().parent
MAPPING=ROOT/'T577_macaque_human_orthologs_Ensembl91.tsv'
ARCHIVE=ROOT/'v10nr_clust_public.zip'
SNAPSHOT='v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl'
OUT=ROOT/'build'

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(8<<20),b''):h.update(chunk)
    return h.hexdigest()

def main():
    aliases={r['native_gene_id']:r for r in csv.DictReader(open(OUT/'gene_aliases.tsv'),delimiter='\t')}
    eligible={}
    reasons=Counter()
    for r in csv.DictReader(open(MAPPING),delimiter='\t'):
        if not r['motif_ids_exact_pwm_overlap']:
            continue
        if r['homology_type']!='ortholog_one2one':
            reasons['non_one_to_one_homology']+=1;continue
        if r['macaque_symbol_locus_ambiguous']!='false':
            reasons['ambiguous_native_locus']+=1;continue
        native=r['macaque_gene_id']
        assert native in aliases
        if aliases[native]['ranking_gene']!=r['T577_gene_symbol']:
            reasons['native_alias_not_input_symbol']+=1;continue
        key=r['human_gene_name']
        pair=(r['human_gene_id'],r['macaque_gene_id'],r['T577_gene_symbol'])
        if key in eligible:
            assert pair==(eligible[key]['human_gene_id'],eligible[key]['macaque_gene_id'],eligible[key]['T577_gene_symbol']),('ambiguous transfer',key)
        eligible[key]=r
    with zipfile.ZipFile(ARCHIVE) as z:
        motifs={Path(n).stem:n for n in z.namelist() if n.endswith('.cb')}
        assert len(motifs)==10249
        motif_dir=Path(os.environ.get('T577_MOTIF_SCRATCH','/tmp/T577_macaque_motifs'));motif_dir.mkdir(parents=True,exist_ok=True)
        for motif,member in motifs.items():
            path=motif_dir/(motif+'.cb')
            data=z.read(member)
            if path.exists():assert path.read_bytes()==data
            else:path.write_bytes(data)
        (OUT/'motif_ids.txt').write_text('\n'.join(sorted(motifs))+'\n')
        with z.open(SNAPSHOT) as raw:
            rows=csv.DictReader(io.TextIOWrapper(raw),delimiter='\t')
            fields=list(rows.fieldnames)
            extra=['human_snapshot_gene_name','human_snapshot_orthologous_identity','human_snapshot_orthologous_gene_name','human_snapshot_orthologous_species','macaque_gene_id','human_gene_id','ensembl_compara_release','ensembl_homology_type','macaque_percent_identity','human_percent_identity']
            output=OUT/'motifs-v10nr-public.macaca_fascicularis.Ensembl91.one2one.tbl'
            tf_names=set();annotated_motifs=set();n=0
            with open(output,'w') as f:
                writer=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t',lineterminator='\n');writer.writeheader()
                for r in rows:
                    human=r['gene_name']
                    if r['#motif_id'] not in motifs or human not in eligible:continue
                    orth=eligible[human]
                    r.update(human_snapshot_gene_name=human,human_snapshot_orthologous_identity=r['orthologous_identity'],human_snapshot_orthologous_gene_name=r['orthologous_gene_name'],human_snapshot_orthologous_species=r['orthologous_species'])
                    r['gene_name']=orth['T577_gene_symbol']
                    r['orthologous_gene_name']=human
                    r['orthologous_species']='homo_sapiens'
                    # Compara fields are percentages, pySCENIC expects an identity fraction.
                    r['orthologous_identity']=str(min(float(orth['macaque_percent_identity']),float(orth['human_percent_identity']))/100)
                    r['description']+='; transferred to Macaca fascicularis using Ensembl Compara 91 one-to-one orthology '+orth['human_gene_id']+' -> '+orth['macaque_gene_id']
                    r.update(macaque_gene_id=orth['macaque_gene_id'],human_gene_id=orth['human_gene_id'],ensembl_compara_release='91',ensembl_homology_type=orth['homology_type'],macaque_percent_identity=orth['macaque_percent_identity'],human_percent_identity=orth['human_percent_identity'])
                    writer.writerow(r);tf_names.add(r['gene_name']);annotated_motifs.add(r['#motif_id']);n+=1
    assert len(tf_names)==len(eligible)==1151, (len(tf_names),len(eligible))
    tf_file=OUT/'T577_macaca_fascicularis.Ensembl91.one2one.tfs.txt'
    tf_file.write_text('\n'.join(sorted(tf_names))+'\n')
    report=dict(status='ready_for_native_ranking_build',n_public_motifs=len(motifs),n_tfs=len(tf_names),n_annotation_rows=n,n_annotated_motifs=len(annotated_motifs),excluded_motif_candidate_rows=dict(reasons),policy='T577 input symbols with unique native GTF loci, Ensembl91 ortholog_one2one, available public PWM exact-ID overlap; human TF experimental motifs transferred by conserved orthology; native binding not directly experimentally validated',mapping_source=str(MAPPING),mapping_sha256=sha256(MAPPING),archive_sha256=sha256(ARCHIVE),snapshot_member=SNAPSHOT,motif_annotation_file=str(output),motif_annotation_sha256=sha256(output),tf_list_file=str(tf_file),tf_list_sha256=sha256(tf_file),motif_ids_file=str(OUT/'motif_ids.txt'),motif_dir=str(motif_dir))
    (OUT/'tf_annotation_preparation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
