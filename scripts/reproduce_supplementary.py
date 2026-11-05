#!/usr/bin/env python3
"""Rebuild S1/S2/S3/S4/S5 from archived source tables, leaving originals unchanged."""
from pathlib import Path
import argparse, shutil
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--output',default='reproduced/supplementary');args=ap.parse_args()
dest=Path(args.output).resolve();dest.mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(ROOT/p,sep='\t',dtype=str,keep_default_na=False)
def save(d,name):d.to_csv(dest/name,sep='\t',index=False,na_rep='NA')
f=read('results/curation/AVG_final_checkv90_complete.tsv.gz')
m=read('figures/figure2/components/panel_A/inputs/functional_mapping.tsv')
v=read('data/contig_master_index.tsv.gz')
bg=read('results/phatyp/eligible_background.tsv')
cols=['Protein','Contig','AVG_category','Viral_Origin_Confidence','Function','reviewed_gene_name','review_status','review_reason','KEGG_hmm_id','Pfam_hmm_id','PHROG_hmm_id','orf_start','orf_end','strand','partial','protein_length_aa','orf_qc_decision','orf_qc_reason','within_100bp_of_edge','contig_length','checkv_completeness','checkv_completeness_method','virus_score','n_hallmarks']
s=f[cols].merge(m[['Protein','Display_module','Mapping_evidence_or_exclusion','Lifestyle']],on='Protein',validate='one_to_one')
s=s.merge(v[['Contig','PhaTYP_pred','PhaTYP_score','PhaGCN3_class','PhaGCN3_family','host_prediction','host_confidence','host_methods']],on='Contig',validate='many_to_one')
s['In_lifestyle_background']=s.Contig.isin(set(bg.Contig))
save(s,'Table_S1_AVG_catalogue.tsv')
save(v[v.Contig.isin(set(f.Contig))],'Table_S1_vOTU_metadata.tsv')
save(bg,'Table_S1_lifestyle_background.tsv')
counts=read('figures/figure2/components/panel_A/functional_module_protein_counts.tsv')
definitions=[];annotations=[]
for _,r in counts.iterrows():
    sub=m[(m.AVG_category==r['Class'])&(m.Display_module==r['Module'])]
    definitions.append([r['Class'],r['Module'],len(sub),sub.Contig.nunique(),sub.Function.nunique(),' | '.join(sorted(set(sub.Mapping_evidence_or_exclusion)))])
    for fn,g in sub.groupby('Function',sort=True):annotations.append([r['Class'],r['Module'],fn,len(g),g.Contig.nunique(),' | '.join(sorted(set(g.KEGG_hmm_id)-{'NA',''}))])
save(pd.DataFrame(definitions,columns=['AVG_class','Module','Proteins','vOTUs','Distinct_Function_annotations','Recorded_assignment_basis']),'Table_S2_module_definitions.tsv')
save(pd.DataFrame(annotations,columns=['AVG_class','Module','Function','Proteins','vOTUs','Observed_KEGG_IDs']),'Table_S2_annotation_mapping.tsv')
shutil.copyfile(ROOT/'figures/figure2/components/panel_A/inputs/functional_mapping.tsv',dest/'Table_S2_protein_assignments.tsv')
shutil.copyfile(ROOT/'figures/figure2/components/panel_B/all_module_tests.tsv',dest/'Table_S2_lifestyle_tests.tsv')
for source,target in [('Table_S3_Sequence_manifest.tsv','Table_S3_phylogenetic_sequences.tsv'),('Table_S4_Model_selection.tsv','Table_S4_model_selection.tsv')]:
    d=read('results/phylogeny/server_audit/'+source);save(d[d.Analysis_ID!='Fhs_expanded52'],target)
for votu,gene in [(1650,'Fhs'),(17468,'FolD'),(8476,'NrdG')]:
    shutil.copyfile(ROOT/f'figures/figure4/composite/sources/vOTU{votu}_gene_context.tsv',dest/f'Table_S5_{gene}_vOTU{votu}_context.tsv')
print('Recreated 12 supplementary TSV files in '+str(dest))
