#!/usr/bin/env python3
"""Validate catalogue, modules, Fisher tests, tree tips and archived checksums."""
from pathlib import Path
import json, hashlib, math, csv
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from Bio import Phylo, SeqIO

ROOT=Path(__file__).resolve().parents[1]
def read(p):return pd.read_csv(ROOT/p,sep='\t',keep_default_na=False,dtype=str)
checks=[]
def check(name,ok,details=''):
    checks.append({'check':name,'passed':bool(ok),'details':str(details)})

s1=read('supplementary/Table_S1_AVG_catalogue.tsv')
check('Final catalogue unique proteins',len(s1)==s1.Protein.nunique()==11725)
check('Final AVG-positive vOTUs',s1.Contig.nunique()==5384)
check('AVG class counts',s1.AVG_category.value_counts().to_dict()=={'AMG':3028,'APG':1808,'AReG':6889})
check('Distinct Function labels per class',s1.groupby('AVG_category').Function.nunique().to_dict()=={'AMG':229,'APG':101,'AReG':158})
bg=read('supplementary/Table_S1_lifestyle_background.tsv')
check('Eligible background',len(bg)==bg.Contig.nunique()==9137)
check('Background lifestyle denominators',bg.PhaTYP_Pred.str.lower().value_counts().to_dict()=={'virulent':6061,'temperate':3076})
sel=s1[s1.Contig.isin(bg.Contig)]
check('Lifestyle-eligible final proteins and vOTUs',len(sel)==10604 and sel.Contig.nunique()==4860)
lc=sel[['Contig','Lifestyle']].drop_duplicates().Lifestyle.value_counts().to_dict()
check('AVG-positive lifestyle counts',lc=={'Virulent':2811,'Temperate':2049},lc)
mapping=read('supplementary/Table_S2_protein_assignments.tsv')
mods=read('supplementary/Table_S2_module_definitions.tsv')
check('Module count',len(mods)==32)
check('Module-assignment protein keys',set(mapping.Protein)==set(s1.Protein))
check('Mapped protein count',sum(~mapping.Display_module.isin(['NA','']))==7979)
tests=read('supplementary/Table_S2_lifestyle_tests.tsv')
ps=[]; count_issues=[];stat_issues=[]
for _,row in tests.iterrows():
    subset=mapping[(mapping.AVG_category==row.Class)&(mapping.Display_module==row.Module)]
    ids=set(subset.Contig)
    v=int(bg.loc[bg.PhaTYP_Pred.str.lower()=='virulent','Contig'].isin(ids).sum())
    t=int(bg.loc[bg.PhaTYP_Pred.str.lower()=='temperate','Contig'].isin(ids).sum())
    if v!=int(row.Virulent_k) or t!=int(row.Temperate_k):count_issues.append(row.Module)
    table=[[v,6061-v],[t,3076-t]]
    p=fisher_exact(table,alternative='two-sided').pvalue;ps.append(p)
    est=odds_ratio(table,kind='conditional');ci=est.confidence_interval(confidence_level=.95)
    observed=[float(row.OR),float(row.CI_lower),float(row.CI_upper),float(row.p)]
    expected=[est.statistic,ci.low,ci.high,p]
    # R fisher.test and SciPy use different numerical root-finding tolerances.
    # R-to-R equality is checked separately in validation/R_statistics_comparison.tsv.
    if not np.allclose(observed,expected,rtol=1e-3,atol=1e-12):stat_issues.append(row.Module)
    m=mods[(mods.AVG_class==row.Class)&(mods.Module==row.Module)].iloc[0]
    check('Module totals: '+row.Module,len(subset)==int(m.Proteins) and subset.Contig.nunique()==int(m.vOTUs))
check('Module lifestyle counts reconstructed from catalogue',not count_issues,count_issues)
check('Fisher P and conditional OR/exact CI independently reproduced (cross-engine rtol=0.001)',not stat_issues,stat_issues)
ps=np.array(ps);order=np.argsort(ps);q=np.empty(len(ps));q[order]=np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1].clip(0,1)
check('BH adjustment across all 32 modules',np.allclose(q,tests.q.astype(float),rtol=1e-6,atol=1e-14))
s3=read('supplementary/Table_S3_phylogenetic_sequences.tsv');s4=read('supplementary/Table_S4_model_selection.tsv')
check('Eight main phylogenetic analyses',len(s4)==8 and len(s3)==547)
for _,r in s4.iterrows():
    aid=r.Analysis_ID;folder=ROOT/'results/phylogeny'/aid
    trees=list(folder.glob('*.treefile'))
    tr=Phylo.read(trees[0],'newick');tips={x.name for x in tr.get_terminals()};sm=s3[s3.Analysis_ID==aid]
    check(aid+' tree IDs and manifest',tips==set(sm.Tree_tip_ID) and len(tips)==int(r.Tree_tips))
    align=list(folder.glob('*.trimmed.faa'))
    if align:
        rec=list(SeqIO.parse(align[0],'fasta'))
        check(aid+' trimmed alignment IDs and length',set(x.id for x in rec)==tips and {len(x.seq) for x in rec}=={int(r.Alignment_length_after_trimming_aa)})
check('Fhs main tree has five viral candidates',int(s4.loc[s4.Analysis_ID=='Fhs_high5','Viral_candidate_tips'].iloc[0])==5)
manifest=ROOT/'FILE_MANIFEST.tsv'
if manifest.exists():
    failed=[]
    for r in csv.DictReader(manifest.open(),delimiter='\t'):
        p=ROOT/r['Path'];h=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else ''
        if h!=r['SHA256']:failed.append(r['Path'])
    check('File checksums',not failed,failed)
report={'passed':all(x['passed'] for x in checks),'checks':checks,'limitations':['Current Fhs_high5 FASTA, alignments and IQ-TREE logs are missing; its tree and sequence manifest are present.','Initial 40,314-candidate protein-level source is absent; only archived summary counts are available.','Figure 3 exact manuscript PNG is preserved; the archived R drawing is a reconstruction.']}
p=ROOT/'validation/report.json';p.parent.mkdir(exist_ok=True);p.write_text(json.dumps(report,indent=2))
print(json.dumps({'passed':report['passed'],'checks':len(checks),'failures':[x for x in checks if not x['passed']]},indent=2))
raise SystemExit(0 if report['passed'] else 1)
