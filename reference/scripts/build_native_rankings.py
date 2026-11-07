#!/usr/bin/env python3
"""Build native gene rankings with official cisTarget scoring and per-motif checkpoints."""

import argparse
import concurrent.futures as futures
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading
import time


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    os.replace(temporary, path)


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--motifs-dir", required=True)
    parser.add_argument("--motif-list", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--source-dir", default=str(Path(__file__).parent / "source/create_cisTarget_databases"))
    parser.add_argument("--deps-dir", required=True, help="Isolated directory containing flatbuffers; does not modify the SpaGRN environment.")
    parser.add_argument("--gene-regex", default=r"#[0-9]+$")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--motif-timeout", type=float, default=1200, help="Deadline per cbust invocation; a chunked motif is bounded by total-timeout.")
    parser.add_argument("--total-timeout", type=float, default=86400)
    parser.add_argument("--seed", type=int, default=20261014)
    parser.add_argument("--complex-pwm-threshold", type=int, default=8)
    parser.add_argument("--chunk-bases", type=int, default=10000000)
    args = parser.parse_args()
    if args.workers < 1 or args.motif_timeout <= 0 or args.total_timeout <= 0 or args.chunk_bases < 1 or args.complex_pwm_threshold < 0:
        parser.error("workers, chunk size and timeouts must be positive; PWM threshold must be nonnegative")

    def total_deadline(signum, frame):
        raise TimeoutError("Total build deadline exceeded")

    signal.signal(signal.SIGALRM, total_deadline)
    signal.setitimer(signal.ITIMER_REAL, max(0.001, args.total_timeout - (time.monotonic() - started)))

    # Keep optional dependencies and JIT/cache writes outside the inference environment.
    sys.dont_write_bytecode = True
    sys.path[:0] = [str(Path(args.deps_dir).resolve()), str(Path(args.source_dir).resolve())]
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/spagrn_native_cistarget_numba_cache")
    import numpy as np
    import pandas as pd
    from clusterbuster import get_motif_id_to_filename_and_nbr_motifs_dict, run_cluster_buster_for_motif
    from cistarget_db import CisTargetDatabase, DatabaseTypes, MotifOrTrackIDs, MotifsOrTracksType, RegionOrGeneIDs

    source = Path(args.source_dir).resolve()
    fasta = Path(args.fasta).resolve()
    motif_list = Path(args.motif_list).resolve()
    prefix = Path(args.output_prefix).resolve()
    checkpoint = prefix.with_name(prefix.name + ".checkpoints")
    checkpoint.mkdir(parents=True, exist_ok=True)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cbust = Path(__file__).parent / "cbust"
    if not cbust.is_file() or not os.access(cbust, os.X_OK):
        raise RuntimeError("Downloaded cbust binary is missing or not executable: " + str(cbust))
    gene_ids = RegionOrGeneIDs.get_region_or_gene_ids_from_fasta(str(fasta), args.gene_regex)
    genes = list(gene_ids.ids)
    if not genes:
        raise RuntimeError("FASTA does not define a gene background")
    motif_files, pwm_counts = get_motif_id_to_filename_and_nbr_motifs_dict(str(Path(args.motifs_dir).resolve()), str(motif_list))
    if not motif_files:
        raise RuntimeError("No motif files were selected")
    fasta_sha = sha256(fasta)
    chunks = []
    if any(count > args.complex_pwm_threshold for count in pwm_counts.values()):
        chunk_dir = checkpoint / "fasta_chunks"
        chunk_dir.mkdir(exist_ok=True)
        chunk_meta_path = chunk_dir / "chunks.json"
        chunk_config = {"fasta_sha256": fasta_sha, "max_bases": args.chunk_bases, "gene_regex": args.gene_regex}
        if chunk_meta_path.exists():
            chunk_meta = json.loads(chunk_meta_path.read_text())
            if chunk_meta["inputs"] != chunk_config:
                raise RuntimeError("Existing FASTA chunks belong to different inputs; use a new output prefix")
            chunks = chunk_meta["chunks"]
            for chunk in chunks:
                if sha256(chunk_dir / chunk["filename"]) != chunk["sha256"]:
                    raise RuntimeError("Existing FASTA chunk hash mismatch: " + chunk["filename"])
            if chunk_meta["genes"] != genes:
                raise RuntimeError("FASTA chunks and original gene background differ")
        else:
            chunk_genes = set()
            handle = None
            chunk_bases = 0
            chunk_records = 0
            chunk_path = None

            def finish_chunk():
                if handle is not None:
                    handle.close()
                    os.replace(str(chunk_path) + ".pending", chunk_path)
                    chunks.append({"filename": chunk_path.name, "sha256": sha256(chunk_path), "bases": chunk_bases, "records": chunk_records})

            def write_record(header, lines):
                nonlocal handle, chunk_bases, chunk_records, chunk_path
                bases = sum(len(line.strip()) for line in lines)
                if bases == 0:
                    raise RuntimeError("Empty FASTA record: " + header.strip())
                if handle is not None and chunk_bases + bases > args.chunk_bases:
                    finish_chunk()
                    handle = None
                if handle is None:
                    chunk_path = chunk_dir / ("chunk_%04d.fa" % len(chunks))
                    handle = open(str(chunk_path) + ".pending", "w")
                    chunk_bases = 0
                    chunk_records = 0
                handle.write(header)
                handle.writelines(lines)
                chunk_bases += bases
                chunk_records += 1
                chunk_genes.add(re.sub(args.gene_regex, "", header[1:].split()[0]))

            header, lines = None, []
            with open(fasta) as original:
                for line in original:
                    if line.startswith(">"):
                        if header is not None:
                            write_record(header, lines)
                        header, lines = line, []
                    elif header is None and line.strip():
                        raise RuntimeError("Sequence data precedes the first FASTA header")
                    else:
                        lines.append(line)
                if header is not None:
                    write_record(header, lines)
            finish_chunk()
            if chunk_genes != set(genes):
                raise RuntimeError("FASTA chunks and original gene background differ")
            atomic_json(chunk_meta_path, {"inputs": chunk_config, "chunks": chunks, "genes": genes})
    genes_sha = hashlib.sha256("\n".join(genes).encode()).hexdigest()
    identity = {
        "fasta_sha256": fasta_sha, "motif_list_sha256": sha256(motif_list),
        "cbust_sha256": sha256(cbust), "gene_ids_sha256": genes_sha,
        "gene_regex": args.gene_regex, "seed": args.seed,
        "scoring_parameters": {"f": 4, "c": 0.0, "r": 10000, "b": 0, "t": 1, "mask": False},
        "source_sha256": {name: sha256(source / name) for name in ["clusterbuster.py", "cistarget_db.py", "orderstatistics.py", "feather_v1_or_v2.py"]},
        "motifs": [{"id": mid, "sha256": sha256(path), "pwm_count": pwm_counts[mid]} for mid, path in sorted(motif_files.items())],
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__},
        "chunking": {"complex_pwm_threshold": args.complex_pwm_threshold, "chunk_bases": args.chunk_bases, "chunks": chunks},
    }
    identity_sha = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    identity_path = checkpoint / "input_identity.json"
    if identity_path.exists():
        previous = json.loads(identity_path.read_text())
        if previous.get("identity_sha256") != identity_sha:
            raise RuntimeError("Existing checkpoints belong to different inputs or parameters; use a new output prefix")
    else:
        atomic_json(identity_path, {"identity_sha256": identity_sha, "inputs": identity})

    def checkpoint_path(mid):
        return checkpoint / (hashlib.sha256(mid.encode()).hexdigest() + ".npz")

    def read_checkpoint(mid):
        with np.load(checkpoint_path(mid), allow_pickle=False) as saved:
            if saved["motif_id"].item() != mid or saved["identity_sha256"].item() != identity_sha:
                raise RuntimeError("Checkpoint identity mismatch: " + mid)
            scores = saved["scores"]
            if scores.dtype != np.float32 or scores.shape != (len(genes),) or not np.isfinite(scores).all() or (scores < 0).any():
                raise RuntimeError("Invalid checkpoint score vector: " + mid)
            return scores

    completed = []
    pending = []
    for mid in motif_files:
        if checkpoint_path(mid).exists():
            read_checkpoint(mid)
            completed.append(mid)
        else:
            pending.append(mid)
    print(json.dumps({"stage": "inputs_ready", "genes": len(genes), "motifs": len(motif_files), "resumed": len(completed), "pending": len(pending), "workers": args.workers, "identity_sha256": identity_sha}), flush=True)

    # The official scorer uses Popen.communicate() without a deadline. This wrapper
    # adds a subprocess deadline and kills only cbust process groups created here.
    original_popen = subprocess.Popen
    active = {}
    lock = threading.Lock()
    stop = threading.Event()

    class DeadlinePopen(original_popen):
        def __init__(self, *positional, **keywords):
            with lock:
                if stop.is_set():
                    raise RuntimeError("Build stopped before starting another motif")
                keywords["start_new_session"] = True
                super().__init__(*positional, **keywords)
                active[self.pid] = self

        def communicate(self, *positional, **keywords):
            keywords["timeout"] = args.motif_timeout
            try:
                return super().communicate(*positional, **keywords)
            except subprocess.TimeoutExpired as exc:
                try:
                    os.killpg(self.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                super().communicate(timeout=10)
                raise RuntimeError("cbust invocation exceeded its deadline") from exc
            finally:
                with lock:
                    active.pop(self.pid, None)

    def cancel_active():
        stop.set()
        with lock:
            for pid in list(active):
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

    def score_one(mid):
        if stop.is_set():
            raise RuntimeError("Build stopped")
        motif_started = time.monotonic()
        chunked = pwm_counts[mid] > args.complex_pwm_threshold
        sequence_files = [checkpoint / "fasta_chunks" / item["filename"] for item in chunks] if chunked else [fasta]
        scores = np.zeros(len(genes), dtype=np.float32)
        partial_path = checkpoint_path(mid).with_suffix(".partial.npz")
        last_completed = -1
        if chunked and partial_path.exists():
            with np.load(partial_path, allow_pickle=False) as saved:
                if saved["motif_id"].item() != mid or saved["identity_sha256"].item() != identity_sha:
                    raise RuntimeError("Partial motif checkpoint identity mismatch: " + mid)
                scores = saved["aggregate_scores"]
                last_completed = int(saved["last_completed_chunk"].item())
                if scores.dtype != np.float32 or scores.shape != (len(genes),) or not np.isfinite(scores).all() or (scores < 0).any() or not (-1 <= last_completed < len(sequence_files)):
                    raise RuntimeError("Invalid partial motif checkpoint: " + mid)
        for chunk_index, sequence_file in enumerate(sequence_files):
            if chunk_index <= last_completed:
                continue
            if stop.is_set():
                raise RuntimeError("Build stopped")
            returned_mid, frame = run_cluster_buster_for_motif(str(cbust), str(sequence_file), motif_files[mid], mid, args.gene_regex)
            if returned_mid != mid or not frame.index.is_unique or not set(frame.index).issubset(set(genes)):
                raise RuntimeError("Unexpected motif ID or unknown/duplicate gene IDs in cbust output: " + mid)
            returned_scores = frame["crm_score"]
            if not np.isfinite(returned_scores.to_numpy()).all() or (returned_scores < 0).any():
                raise RuntimeError("Nonfinite or negative CRM score: " + mid)
            # The official scorer omits zero-hit genes. Across chunks, taking
            # the maximum preserves its maximum-over-TSS gene aggregation.
            np.maximum(scores, returned_scores.reindex(genes, fill_value=0).to_numpy(dtype=np.float32), out=scores)
            if chunked:
                temporary_partial = partial_path.with_name(partial_path.name + ".tmp-" + str(threading.get_ident()))
                with open(temporary_partial, "wb") as handle:
                    np.savez_compressed(handle, motif_id=mid, identity_sha256=identity_sha, aggregate_scores=scores, last_completed_chunk=chunk_index)
                os.replace(temporary_partial, partial_path)
                print(json.dumps({"stage": "motif_chunk_completed", "motif": mid, "chunk": chunk_index + 1, "chunks": len(sequence_files)}), flush=True)
        destination = checkpoint_path(mid)
        temporary = destination.with_name(destination.name + ".tmp-" + str(threading.get_ident()))
        with open(temporary, "wb") as handle:
            np.savez_compressed(handle, motif_id=mid, identity_sha256=identity_sha, scores=scores)
        os.replace(temporary, destination)
        read_checkpoint(mid)
        if chunked:
            partial_path.unlink()
        return {"motif": mid, "seconds": time.monotonic() - motif_started, "genes_with_positive_scores": int((scores > 0).sum())}

    executor = futures.ThreadPoolExecutor(max_workers=args.workers)
    jobs = {}
    subprocess.Popen = DeadlinePopen
    try:
        for mid in pending:
            jobs[executor.submit(score_one, mid)] = mid
        for future in futures.as_completed(jobs, timeout=max(0.001, args.total_timeout - (time.monotonic() - started))):
            report = future.result()
            completed.append(report["motif"])
            report.update(stage="motif_completed", completed=len(completed), total=len(motif_files))
            print(json.dumps(report), flush=True)
    except BaseException:
        cancel_active()
        for future in jobs:
            future.cancel()
        raise
    finally:
        executor.shutdown(wait=True)
        subprocess.Popen = original_popen
    if len(completed) != len(motif_files):
        raise RuntimeError("Not every requested motif completed; rankings will not be written")
    if time.monotonic() - started > args.total_timeout:
        raise TimeoutError("Total build deadline exceeded before rankings")

    motif_ids = MotifOrTrackIDs(list(motif_files), MotifsOrTracksType.MOTIFS)
    scores_db = CisTargetDatabase.create_db(DatabaseTypes.SCORES_DB_MOTIFS_VS_GENES, gene_ids, motif_ids, order="F")
    for mid in motif_ids.ids:
        scores_db.update_scores_for_motif_or_track(mid, pd.Series(read_checkpoint(mid), index=genes))
    rankings = scores_db.convert_scores_db_to_rankings_db(seed=args.seed)
    expected_ranks = np.arange(len(genes), dtype=rankings.dtype)
    for mid in motif_ids.ids:
        ranks = rankings.df[mid].to_numpy()
        if not np.array_equal(np.sort(ranks), expected_ranks):
            raise RuntimeError("Rankings are not a complete zero-based permutation: " + mid)
        if (np.diff(scores_db.df[mid].to_numpy()[np.argsort(ranks)]) > 0).any():
            raise RuntimeError("Ranking order disagrees with descending motif scores: " + mid)
    del scores_db
    if time.monotonic() - started > args.total_timeout:
        raise TimeoutError("Total build deadline exceeded during ranking validation")
    staged_prefix = str(prefix) + ".pending"
    staged = Path(rankings.transpose().write_db(db_prefix=staged_prefix, version=2))
    if time.monotonic() - started > args.total_timeout:
        raise TimeoutError("Total build deadline exceeded during Feather serialization")
    final = Path(str(prefix) + ".genes_vs_motifs.rankings.feather")
    os.replace(staged, final)
    manifest = {
        "success": True, "identity_sha256": identity_sha, "gene_count": len(genes),
        "motif_count": len(motif_ids), "motif_ids": list(motif_ids.ids),
        "seed": args.seed, "ranking_convention": "Motif rows and gene columns; rank zero is highest CRM score; ties broken with recorded seed",
        "zero_hit_policy": "Missing positive cbust hits have gene score zero",
        "rankings_file": str(final), "rankings_sha256": sha256(final),
        "checkpoint_directory": str(checkpoint), "elapsed_seconds": time.monotonic() - started,
    }
    atomic_json(prefix.with_name(prefix.name + ".build.json"), manifest)
    print(json.dumps({"stage": "complete", **{key: value for key, value in manifest.items() if key != "motif_ids"}}), flush=True)


if __name__ == "__main__":
    main()
