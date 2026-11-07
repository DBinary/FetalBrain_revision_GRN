"""Prepare native transcript-TSS FASTA from Ensembl 91, with explicit gene aliases.
Coordinates: GTF 1-based inclusive -> 0-based TSS; windows include the TSS base.
No gene is removed from the background for a missing or ambiguous gene symbol.
"""
from pathlib import Path
from collections import Counter, defaultdict
import gzip, re, json, csv, hashlib, os

ROOT = Path(__file__).resolve().parent
GTF = ROOT / 'Macaca_fascicularis.Macaca_fascicularis_5.0.91.gtf.gz'
GENOME = ROOT / 'Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa.gz'
OUT = ROOT / 'build'
WINDOWS = [('10kbp_up_10kbp_down', 10000, 10000), ('500bp_up_100bp_down', 500, 100)]

def fasta_records(path):
    name, pieces = None, []
    with gzip.open(path, 'rb') as handle:
        for line in handle:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(pieces)
                name, pieces = line[1:].split()[0].decode(), []
            else:
                pieces.append(line.strip())
        if name is not None:
            yield name, b''.join(pieces)

def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    completion = OUT / 'promoter_preparation.json'
    if completion.exists():
        meta = json.loads(completion.read_text())
        assert meta['source_genome_sha256'] == sha256(GENOME)
        assert meta['source_gtf_sha256'] == sha256(GTF)
        for record in meta['windows']:
            assert sha256(ROOT / record['fasta_relative_path']) == record['sha256']
        print(json.dumps({'reused_verified_promoters': True, 'genes': meta['n_genes']}, ensure_ascii=False), flush=True)
        return
    genes, transcripts = {}, defaultdict(set)
    with gzip.open(GTF, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] not in ('gene', 'transcript'):
                continue
            attr = dict(re.findall(r'(\S+) "([^"]*)";', fields[8]))
            gene = attr['gene_id']
            if fields[2] == 'gene':
                assert gene not in genes, ('duplicate native gene ID', gene)
                genes[gene] = {'gene_name': attr.get('gene_name', ''), 'biotype': attr.get('gene_biotype', ''), 'chromosome': fields[0], 'strand': fields[6]}
            else:
                assert fields[6] in ('+', '-')
                tss = int(fields[3]) - 1 if fields[6] == '+' else int(fields[4]) - 1
                transcripts[gene].add((fields[0], tss, fields[6]))
    assert set(genes) == set(transcripts), ('gene/transcript mismatch', len(set(genes) ^ set(transcripts)))
    names = Counter(record['gene_name'] for record in genes.values() if record['gene_name'])
    aliases = {gene: (record['gene_name'] if record['gene_name'] and names[record['gene_name']] == 1 and not any(c.isspace() for c in record['gene_name']) and '#' not in record['gene_name'] else gene) for gene, record in genes.items()}
    assert len(set(aliases.values())) == len(genes)
    assert all(not any(c.isspace() for c in alias) and '#' not in alias for alias in aliases.values())
    by_chromosome = defaultdict(list)
    for gene, locations in transcripts.items():
        for i, (chromosome, tss, strand) in enumerate(sorted(locations), 1):
            by_chromosome[chromosome].append((gene, aliases[gene] + '#' + str(i), tss, strand))
    for locations in by_chromosome.values():
        locations.sort(key=lambda record: (record[2], record[1]))
    handles = {label: open(OUT / (label + '.transcript_tss.fa.partial'), 'wb') for label, _, _ in WINDOWS}
    stats = {label: {'sequences': 0, 'bases': 0, 'clipped_to_sequence_boundary': 0, 'genes': set()} for label, _, _ in WINDOWS}
    seen = set()
    complement = bytes.maketrans(b'ACGTRYKMSWBDHVNacgtrykmswbdhvn', b'TGCAYRMKSWVHDBNtgcayrmkswvhdbn')
    try:
        for chromosome, sequence in fasta_records(GENOME):
            if chromosome not in by_chromosome:
                continue
            seen.add(chromosome)
            for gene, header, tss, strand in by_chromosome[chromosome]:
                assert 0 <= tss < len(sequence), (gene, chromosome, tss, len(sequence))
                for label, upstream, downstream in WINDOWS:
                    lo, hi = (tss-upstream, tss+downstream+1) if strand == '+' else (tss-downstream, tss+upstream+1)
                    start, end = max(0, lo), min(len(sequence), hi)
                    seq = sequence[start:end]
                    if strand == '-':
                        seq = seq.translate(complement)[::-1]
                    assert seq and len(seq) == end-start
                    handles[label].write(b'>' + header.encode() + b'\n' + seq.upper() + b'\n')
                    stats[label]['sequences'] += 1
                    stats[label]['bases'] += len(seq)
                    stats[label]['clipped_to_sequence_boundary'] += int(lo != start or hi != end)
                    stats[label]['genes'].add(gene)
            print('Extracted', chromosome, len(by_chromosome[chromosome]), 'unique gene/TSS locations', flush=True)
        assert seen == set(by_chromosome), ('GTF chromosomes absent from FASTA', sorted(set(by_chromosome)-seen))
    finally:
        for handle in handles.values():
            handle.close()
    with open(OUT / 'gene_aliases.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['native_gene_id','ranking_gene','gene_name','gene_biotype','chromosome','strand','n_unique_tss','alias_policy'])
        for gene in sorted(genes):
            rec = genes[gene]
            policy = 'unambiguous_native_gene_name' if aliases[gene] != gene else ('missing_gene_name_use_native_id' if not rec['gene_name'] else ('ambiguous_gene_name_use_native_id' if names[rec['gene_name']] > 1 else 'non_fasta_gene_name_use_native_id'))
            writer.writerow([gene,aliases[gene],rec['gene_name'],rec['biotype'],rec['chromosome'],rec['strand'],len(transcripts[gene]),policy])
    manifest = {'assembly':'Macaca_fascicularis_5.0','ensembl_release':91, 'n_genes':len(genes),'n_distinct_gene_tss_locations':sum(map(len,transcripts.values())), 'no_gene_names':sum(not rec['gene_name'] for rec in genes.values()),'ambiguous_gene_names':sum(n>1 for n in names.values()),'source_genome_sha256':sha256(GENOME),'source_gtf_sha256':sha256(GTF),'coordinate_convention':'0-based TSS; [tss-upstream,tss+downstream+1) on plus strand; reversed on minus strand; clip to sequence boundaries; identical gene/TSS locations deduplicated','gene_alias_policy':'use unique native gene_name; otherwise retain native Ensembl stable gene ID; gene background includes all GTF gene IDs','windows':[]}
    for label, upstream, downstream in WINDOWS:
        stat = stats[label]
        assert stat.pop('genes') == set(genes)
        path = OUT / (label + '.transcript_tss.fa')
        os.replace(str(path)+'.partial', path)
        manifest['windows'].append(dict(name=label,upstream_bp=upstream,downstream_bp=downstream,fasta_relative_path=str(path.relative_to(ROOT)),sha256=sha256(path),bytes=path.stat().st_size,**stat))
    completion.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2), flush=True)

if __name__ == '__main__':
    main()
