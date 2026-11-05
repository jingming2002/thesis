#!/usr/bin/env python3
"""Build current sequence-manifest and model-selection tables from final trees."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path('/lustre/BIF/nobackup/wang500/thesis1')
OUT = ROOT / 'results/avg_phylogeny_supplement_20261104'
OLD = ROOT / 'results/avg_phylogeny_supplement_20261103/Sequence_manifest.tsv'

RUNS = [
    ('Fhs_high5', 'Fhs', 'results/fhs_ncbi_top5_20261103/high5'),
    ('Fhs_expanded52', 'Fhs', 'results/fhs_ncbi_top5_20261103/expanded52'),
    ('FolD', 'FolD', 'results/amg_family_representative_trees_20261025'),
    ('NrdG', 'NrdG', 'results/amg_family_representative_trees_20261025'),
    *[(g, g, f'results/network5_rebuild_20261026/{g}')
      for g in ('FolE', 'GlyA', 'IscS', 'IscU_like', 'MetK')],
]

def read_tsv(path: Path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path: Path, rows, fields):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def fasta(path: Path):
    out, key, seq = {}, None, []
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if key is not None: out[key] = ''.join(seq)
                key = line[1:].strip().split()[0]; seq = []
            else: seq.append(line.strip())
    if key is not None: out[key] = ''.join(seq)
    return out

def tree_tips(path: Path):
    s = path.read_text().strip()
    return [x.strip().strip("'\"") for x in re.findall(r'(?:\(|,)([^():,;]+):', s)]

def read_clusters(path: Path | None):
    members = {}
    if path is None or not path.is_file(): return members
    cluster = ''
    with path.open() as fh:
        for line in fh:
            if line.startswith('>Cluster '): cluster = line.strip()[1:]
            else:
                if '>' in line and '...' in line:
                    member = line.split('>', 1)[1].split('...', 1)[0]
                    members[member] = cluster
    return members

def source_class(s):
    if s in ('XVIR', 'Viral candidate'): return 'Viral candidate'
    if s in ('NCBI', 'NCBI RefSeq'): return 'NCBI RefSeq'
    if s == 'PIBAC' or s == 'Pig gut reference': return 'Pig-gut reference (PiBAC-associated)'
    return 'Other reference'

old_rows = read_tsv(OLD)
old_by_tip = {(r['Gene_or_protein_family'], r['Tree_tip_ID']): r for r in old_rows}
manifest, summary = [], []
for run_id, family, rel in RUNS:
    d = ROOT / rel
    if family in ('FolD', 'NrdG'):
        mapping_path = d / f'{family}.tree_id_mapping.tsv'
        full_path = d / f'{family}.tree_input.faa'
        raw_aln = d / f'{family}.mafft.faa'
        trimmed = d / f'{family}.trimmed.faa'
        report = d / f'{family}.IQTREE.iqtree'
        tree = d / f'{family}.IQTREE.treefile'
        map_rows = read_tsv(mapping_path)
        def getv(r, *keys): return next((r[k] for k in keys if k in r and r[k] != ''), '')
        cluster_root = ROOT / ('results/folD_all125_local_blast_20261024' if family == 'FolD' else 'results/nrdg_all612_local_blast_20261024')
        pig_cluster_path = ROOT / 'results/folD_nrdG_pibac_refs_20261024' / ('FolD.PIBAC_nr70_cov80.faa.clstr' if family == 'FolD' else 'NrdG_putative.PIBAC_nr70_cov80.faa.clstr')
        viral_clusters = read_clusters(cluster_root / ('FolD_nr70_cov80.faa.clstr' if family == 'FolD' else 'NrdG_nr70_cov80.faa.clstr'))
        pig_clusters = read_clusters(pig_cluster_path)
    elif family == 'Fhs':
        mapping_path = d / 'tree_id_mapping.tsv'; full_path = d / 'tree_input.faa'
        raw_aln = d / 'Fhs.mafft.faa'; trimmed = d / 'Fhs.trimmed.faa'
        report = d / 'Fhs.IQTREE.iqtree'; tree = d / 'Fhs.IQTREE.treefile'
        map_rows = read_tsv(mapping_path)
        def getv(r, *keys): return next((r[k] for k in keys if k in r and r[k] != ''), '')
        viral_clusters = {}
        pig_clusters = read_clusters(ROOT / 'results/fhs_high5_uniform_local_besthit_20261103/references/PIBAC_nr70_cov80.faa.clstr')
    else:
        mapping_path = d / 'tree_id_mapping.tsv'; full_path = d / 'tree_input.faa'
        raw_aln = d / f'{family}.mafft.faa'; trimmed = d / f'{family}.trimmed.faa'
        report = d / f'{family}.IQTREE.iqtree'; tree = d / f'{family}.IQTREE.treefile'
        map_rows = read_tsv(mapping_path)
        def getv(r, *keys): return next((r[k] for k in keys if k in r and r[k] != ''), '')
        viral_clusters = read_clusters(d / 'viral_nr70_cov80.faa.clstr')
        pig_clusters = read_clusters(d / 'PIBAC_nr70_cov80.faa.clstr')
    seqs = fasta(full_path)
    mapping = {getv(r, 'tree_id'): r for r in map_rows}
    ref_descriptions = {}
    ncbi_fasta = d / 'NCBI_full_length.faa'
    if family == 'Fhs': ncbi_fasta = ROOT / 'results/fhs_ncbi_top5_20261103/references/NCBI_full_length.faa'
    if ncbi_fasta.is_file():
        with ncbi_fasta.open() as fh:
            for line in fh:
                if line.startswith('>'):
                    h = line[1:].strip()
                    ref_descriptions[h.split()[0].split('.',1)[0]] = h
    pig_taxonomy = {}
    pig_selection = d / 'PIBAC_selection.tsv'
    if family == 'Fhs': pig_selection = ROOT / 'results/fhs_ncbi_top5_20261103/references/PIBAC_selection.tsv'
    if pig_selection.is_file():
        for pr in read_tsv(pig_selection):
            pid = getv(pr, 'protein_id')
            lineage = getv(pr, 'gtdb_taxonomy')
            if pid and lineage: pig_taxonomy[pid] = lineage
    fhs_inventory = {}
    if family == 'Fhs':
        inv = ROOT / 'results/fhs_ncbi_top5_20261103/candidate_inventory.tsv'
        fhs_inventory = {r['protein_id']: r for r in read_tsv(inv)}
    viral_qc = {}
    viral_qc_path = d / f'{family}.viral_selection.tsv'
    if family == 'Fhs': viral_qc_path = None
    if viral_qc_path and viral_qc_path.is_file():
        for vr in read_tsv(viral_qc_path):
            pid = getv(vr, 'protein_id')
            if pid: viral_qc[pid] = vr
    if len(mapping) != len(map_rows): raise ValueError(f'{run_id}: duplicate tree_id in mapping')
    tips, seqids = tree_tips(tree), list(seqs)
    if set(tips) != set(seqs) or len(tips) != len(seqs):
        raise ValueError(f'{run_id}: tree/FASTA tips differ: fasta={len(seqs)} tree={len(tips)}')
    if set(mapping) != set(seqs): raise ValueError(f'{run_id}: mapping/tree_input IDs differ')
    for tip in seqids:
        m = mapping[tip]
        src = source_class(getv(m, 'source'))
        oid = getv(m, 'original_protein_id', 'original_id')
        legacy = old_by_tip.get((family, tip), {})
        votu, genome, organism = '', '', getv(m, 'organism_name') or 'Not recorded'
        tax = getv(m, 'taxon_label', 'viral_family', 'phylum') or 'Not recorded'
        if src == 'Viral candidate':
            votu = oid.split('|', 1)[0] if '|' in oid else oid
            genome = votu
            organism = 'Not applicable (viral sequence; vOTU label used)'
        elif src.startswith('Pig-gut'):
            genome = oid.split('__', 1)[0]
            tax = pig_taxonomy.get(oid, getv(m, 'taxon_label', 'phylum') or legacy.get('Available_taxonomic_assignment', 'Not recorded'))
            organism = legacy.get('Organism_name', 'Not recorded')
        elif src == 'NCBI RefSeq':
            genome = oid
            descr = ref_descriptions.get(oid.split('.',1)[0], '')
            organism_match = re.search(r'\[([^\]]+)\]\s*$', descr)
            organism = organism_match.group(1) if organism_match else legacy.get('Organism_name', 'Not recorded')
            tax = getv(m, 'taxon_label', 'phylum') or legacy.get('Available_taxonomic_assignment', 'Not recorded')
        else:
            genome = oid
            organism = legacy.get('Organism_name', 'Not recorded')
            tax = getv(m, 'taxon_label', 'phylum') or legacy.get('Available_taxonomic_assignment', 'Not recorded')
        confidence = getv(m, 'viral_origin_confidence') or legacy.get('Viral_origin_confidence', '')
        invrow = fhs_inventory.get(oid, {})
        if family == 'Fhs' and src == 'Viral candidate':
            confidence = invrow.get('viral_origin_confidence', confidence)
        if src == 'Viral candidate' and not confidence:
            confidence = 'Not used as a selection filter; category not recorded here'
        cluster_id = 'Not applied'
        if src == 'Viral candidate' and viral_clusters:
            cluster_id = viral_clusters.get(oid, 'No exact cluster membership found')
        elif src.startswith('Pig-gut') and pig_clusters:
            cluster_id = pig_clusters.get(oid, 'No exact cluster membership found')
        integrity = legacy.get('ORF_integrity', '')
        if src == 'Viral candidate' and family == 'Fhs':
            partial = invrow.get('partial', 'Not recorded')
            internal = invrow.get('internal_stop', 'Not recorded')
            integrity = f'partial={partial}; internal_stop={internal}'
        elif src == 'Viral candidate' and oid in viral_qc:
            integrity = viral_qc[oid].get('orf_integrity', 'Not recorded')
        elif src.startswith('Pig-gut'):
            integrity = 'partial=00; no internal stop (reference inclusion criteria)'
        elif src == 'NCBI RefSeq':
            integrity = 'Full-length RefSeq protein retrieved; viral ORF status not applicable'
        elif src == 'Other reference':
            integrity = legacy.get('ORF_integrity', 'Curated/reference protein; viral ORF status not applicable')
        manifest.append({
            'Analysis_ID': run_id, 'Gene_or_protein_family': family,
            'Tree_tip_ID': tip, 'Original_FASTA_header': invrow.get('original_header', legacy.get('Original_FASTA_header', oid or tip)),
            'Sequence_source': src, 'Accession_with_version_or_protein_ID': oid,
            'vOTU_ID_or_reference_genome_ID': genome or 'Not recorded',
            'Organism_name': organism, 'Available_taxonomic_assignment': tax,
            'Protein_length_aa': len(seqs[tip].rstrip('*')),
            'Source_FASTA_file': str(full_path),
            'CDHIT_cluster_ID': cluster_id,
            'Is_outgroup': 'No', 'Viral_origin_confidence': confidence,
            'ORF_integrity': integrity or 'Not assessed in this manifest',
            'Notes': getv(m, 'provenance', 'evidence', 'note') or legacy.get('Notes', '')
        })
    counts = {}
    for m in map_rows:
        k = source_class(getv(m, 'source'))
        counts[k] = counts.get(k, 0) + 1
    report_text = report.read_text()
    model_m = re.search(r'^Best-fit model according to BIC:\s*(\S+)', report_text, re.M)
    input_m = re.search(r'^Input data:\s*(\d+) sequences with (\d+) amino-acid sites', report_text, re.M)
    if not model_m or not input_m: raise ValueError(f'{run_id}: cannot read IQ-TREE model/alignment metadata')
    if len(fasta(raw_aln)) != len(seqs) or len(fasta(trimmed)) != len(seqs): raise ValueError(f'{run_id}: alignment count mismatch')
    prelen = len(next(iter(fasta(raw_aln).values())))
    trimlen = len(next(iter(fasta(trimmed).values())))
    if int(input_m.group(1)) != len(tips) or int(input_m.group(2)) != trimlen: raise ValueError(f'{run_id}: IQ-TREE report differs from alignment')
    if family == 'Fhs': query_count, query_scope = 5, 'Five high viral-origin-confidence Fhs queries; reused local BLAST search. Expanded52 tree adds 47 non-query viral tips.'
    elif family == 'FolD': query_count, query_scope = 86, 'All 86 CD-HIT viral representatives in the saved local BLAST search; only selected viral-tip hits were assembled into this tree.'
    elif family == 'NrdG': query_count, query_scope = 386, 'All 386 CD-HIT viral representatives in the saved local BLAST search; only selected viral-tip hits were assembled into this tree.'
    else:
        query_count = len(fasta(d / 'viral_nr70_cov80.faa'))
        query_scope = 'All viral proteins after 70% identity / 80% bidirectional-coverage CD-HIT, used as NCBI BLAST queries.'
    if family == 'Fhs':
        mafft_version, mafft_settings, trimal_version = '7.525', '--auto --thread 64', '1.5.rev1'
        iqtree_settings = '-m MFP -B 1000 -alrt 1000 --bnni -T AUTO --threads-max 64 -seed 20261025'
    elif family in ('FolD', 'NrdG'):
        mafft_version, mafft_settings, trimal_version = 'Not recorded', '--auto --thread 16', 'Not recorded'
        iqtree_settings = '-m MFP -B 1000 -alrt 1000 --bnni -T 16 -seed 20261025'
    else:
        mafft_version, mafft_settings, trimal_version = '7.525', '--auto --thread 64', 'Not recorded'
        iqtree_settings = '-m MFP -B 1000 -alrt 1000 --bnni -T AUTO --threads-max 64 -seed 20261025'
    summary.append({
        'Analysis_ID': run_id, 'Gene_family': family, 'Final_tree_file': str(tree),
        'BLAST_query_count': query_count, 'BLAST_query_scope': query_scope,
        'Viral_candidate_tips': counts.get('Viral candidate',0),
        'Pig_gut_reference_tips': counts.get('Pig-gut reference (PiBAC-associated)',0),
        'NCBI_reference_tips': counts.get('NCBI RefSeq',0),
        'Other_reference_tips': counts.get('Other reference',0),
        'Total_input_sequences': len(seqs), 'Tree_tips': len(tips),
        'Alignment_length_before_trimming_aa': prelen,
        'Alignment_length_after_trimming_aa': trimlen,
        'ModelFinder_criterion': 'BIC (MFP)', 'Selected_model': model_m.group(1),
        'MAFFT_version_and_settings': f'{mafft_version}; {mafft_settings}',
        'trimAl_version_and_settings': f'{trimal_version}; -automated1',
        'IQ_TREE_version_and_settings': f'3.1.2; {iqtree_settings}',
        'Rooting_or_outgroup': 'Unrooted; no -o specified',
        'Alignment_tree_IDs_exact_match': 'Yes'
    })

OUT.mkdir(parents=True, exist_ok=True)
manifest_fields = ['Analysis_ID','Gene_or_protein_family','Tree_tip_ID','Original_FASTA_header','Sequence_source','Accession_with_version_or_protein_ID','vOTU_ID_or_reference_genome_ID','Organism_name','Available_taxonomic_assignment','Protein_length_aa','Source_FASTA_file','CDHIT_cluster_ID','Is_outgroup','Viral_origin_confidence','ORF_integrity','Notes']
summary_fields = ['Analysis_ID','Gene_family','Final_tree_file','BLAST_query_count','BLAST_query_scope','Viral_candidate_tips','Pig_gut_reference_tips','NCBI_reference_tips','Other_reference_tips','Total_input_sequences','Tree_tips','Alignment_length_before_trimming_aa','Alignment_length_after_trimming_aa','ModelFinder_criterion','Selected_model','MAFFT_version_and_settings','trimAl_version_and_settings','IQ_TREE_version_and_settings','Rooting_or_outgroup','Alignment_tree_IDs_exact_match']
write_tsv(OUT/'Table_S3_Sequence_manifest.tsv', manifest, manifest_fields)
write_tsv(OUT/'Table_S4_Model_selection.tsv', summary, summary_fields)
print(f'Wrote {len(manifest)} sequence rows and {len(summary)} tree rows to {OUT}')
