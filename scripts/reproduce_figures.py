#!/usr/bin/env python3
"""Re-render figures in an isolated copy; preserve all archived manuscript assets."""
from pathlib import Path
import argparse, subprocess, shutil, sys, json, os
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--figures',nargs='+',choices=['1','2','3','4'],default=['1','2','3','4']);ap.add_argument('--output',default='reproduced');ap.add_argument('--rscript',default='Rscript');args=ap.parse_args()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
logs=[]
def run(cmd,cwd):
    r=subprocess.run([str(x) for x in cmd],cwd=cwd,text=True,capture_output=True)
    logs.append({'command':[str(x) for x in cmd],'cwd':str(cwd),'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
    (out/'run_log.json').write_text(json.dumps(logs,indent=2))
    if r.returncode:raise RuntimeError(f'Command failed: {cmd}\n{r.stderr}')
for n in args.figures:
    dest=out/f'figure{n}';shutil.copytree(ROOT/f'figures/figure{n}',dest,dirs_exist_ok=True)
    if n=='1':run([args.rscript,'plot_overview_with_lifestyle.R'],dest)
    elif n=='2':
        run([args.rscript,'plot_repertoire.R','.'],dest/'components/panel_B')
        run([args.rscript,'plot_protein_counts.R','.'],dest/'components/panel_A')
        run([sys.executable,'combine_panels.py'],dest)
    elif n=='3':run([args.rscript,'draw_pathway.R','.'],dest/'editable_reconstruction')
    else:
        run([args.rscript,'plot_AVG_trees.R'],dest/'trees')
        for v in ['1650','17468','8476']:run([sys.executable,f'draw_vOTU{v}_context.py'],dest/'composite/sources')
        shutil.copyfile(dest/'trees/figures/AVG3_midpoint_combined.pdf',dest/'composite/sources/AVG3_midpoint_combined.pdf')
        run([sys.executable,'combine_phylogeny_context.py'],dest/'composite')
    import fitz
    expected={'1':'AVG_overview_with_lifestyle_R.pdf','2':'AVG_functional_composition_and_lifestyle_AB.pdf','3':'editable_reconstruction/one_carbon_pathway.pdf','4':'composite/AVG_phylogeny_and_genomic_context.pdf'}[n]
    pdf=dest/expected
    if not pdf.exists() or pdf.stat().st_size==0:raise RuntimeError(f'Missing/empty output: {pdf}')
    with fitz.open(pdf) as document:
        if len(document)!=1:raise RuntimeError(f'Expected one-page figure: {pdf}')
        document[0].get_pixmap(matrix=fitz.Matrix(.2,.2),alpha=False)
    print('Rendered and checked figure '+n,flush=True)
