"""Compact overview and nrdG neighbourhood, matching the supplied Fhs layout.
Run from any directory: python draw_vOTU8476_context.py
"""
from pathlib import Path
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Polygon,ConnectionPatch

ROOT=Path(__file__).resolve().parent
with (ROOT/'vOTU8476_gene_context.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
    for k in ['start_position','end_position','strandedness','gene_position']:r[k]=int(r[k])
LENGTH=68827
target=next(r for r in rows if r['gene_position']==27)
assert target['start_position']==35928 and target['ko_id']=='K04068'
origin=target['start_position'];lo=-6.;hi=6.5
hallmarks=[r for r in rows if r['virsorter_category'] in ('0','3')]
assert [r['gene_position'] for r in hallmarks]==[4,15,18]
COL={'target':'#873E99','hallmark':'#4479B0','viral':'#D79348','annotated':'#859BAA','unknown':'#DEDEDE'}
def category(r):
    if r['gene_position']==27:return 'target'
    if r['virsorter_category'] in ('0','3'):return 'hallmark'
    if r['gene_position']==16:return 'viral'
    if r['gene_position'] in (19,35):return 'annotated'
    return 'unknown'

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(13,2.65),facecolor='white')
over=fig.add_axes([.13,.66,.82,.28])
over.set(xlim=(0,110),ylim=(-.50,1.0));over.axis('off')
over.add_patch(Rectangle((0,-.035),LENGTH/1000,.07,color='#D5D9DC',lw=0))
over.text(-1.5,0,'vOTU_8476',ha='right',va='center',fontsize=12,fontweight='bold')
over.text(LENGTH/1000+.7,0,'68.8 kb',ha='left',va='center',fontsize=8,color='#666666')
win=(origin/1000+lo,origin/1000+hi)
over.add_patch(Rectangle((win[0],-.10),win[1]-win[0],.20,facecolor='#F0E8F2',edgecolor=COL['target'],lw=.6,zorder=1))
for r in hallmarks+[target]:
    a=(r['start_position']-1)/1000;b=r['end_position']/1000
    over.add_patch(Rectangle((a,-.075),b-a,.15,facecolor=COL[category(r)],lw=0,zorder=2))
for r,x in zip(hallmarks,[4,24,33]):
    mid=(r['start_position']+r['end_position'])/2000
    over.annotate('hallmark',(mid,.08),(x,.58),ha='center',fontsize=8,color=COL['hallmark'],
                  arrowprops={'arrowstyle':'-','lw':.55,'color':COL['hallmark']})
mid=(target['start_position']+target['end_position'])/2000
over.text(mid,.16,'nrdG',ha='center',va='bottom',fontstyle='italic',fontweight='bold',color=COL['target'],fontsize=9)
over.plot([90,100],[-.19,-.19],lw=1,color='#666666')
over.text(95,-.29,'10 kb',ha='center',va='top',fontsize=8,color='#555555')

ax=fig.add_axes([.13,.25,.82,.38]);ax.set(xlim=(lo,hi),ylim=(-.28,.67));ax.set_yticks([])
for s in ['top','right','left']:ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#B9BEC2');ax.spines['bottom'].set_linewidth(.65)
ax.tick_params(axis='x',labelsize=8,length=2,color='#AAAAAA')
ax.set_xticks([-6,-4,-2,0,2,4,6]);ax.set_xlabel('Position relative to nrdG (kb)',fontsize=9,labelpad=5)
ax.axvline(0,ls=':',color='#D6D1D9',lw=.6,zorder=0)
for r in rows:
    a=(r['start_position']-origin)/1000;b=(r['end_position']+1-origin)/1000
    if b<=lo or a>=hi:continue
    s=max(a,lo);e=min(b,hi);h=.072;head=min(.13,.45*(e-s))
    clipped=(r['strandedness']==1 and b>hi) or (r['strandedness']==-1 and a<lo)
    if clipped:pts=[(s,-h),(e,-h),(e,h),(s,h)]
    elif r['strandedness']==1:pts=[(s,-h),(e-head,-h),(e,0),(e-head,h),(s,h)]
    else:pts=[(e,-h),(s+head,-h),(s,0),(s+head,h),(e,h)]
    ax.add_patch(Polygon(pts,facecolor=COL[category(r)],edgecolor='#767676',lw=.5,zorder=3))

labels={16:('Tape measure',-.0),18:('Portal\n(hallmark)',0),19:('dUTPase',0),27:('nrdG',0),35:('RNA polymerase-like',0)}
for r in rows:
    if r['gene_position'] not in labels:continue
    a=(r['start_position']-origin)/1000;b=(r['end_position']+1-origin)/1000
    x=(max(lo,a)+min(hi,b))/2;lab=labels[r['gene_position']][0]
    y=.14
    if r['gene_position']==16:x=lo+.55
    ax.text(x,y,lab,ha='center',va='bottom',fontsize=9,
            color=COL['target'] if r['gene_position']==27 else (COL['hallmark'] if category(r)=='hallmark' else '#333333'),
            fontweight='bold' if r['gene_position']==27 else 'normal',fontstyle='italic' if r['gene_position']==27 else 'normal')
for x,xx in zip(win,[lo,hi]):
    fig.add_artist(ConnectionPatch((x,-.105),(xx,.63),coordsA=over.transData,coordsB=ax.transData,color='#D3CDD7',lw=.65,zorder=0))
for ext in ['pdf','svg','png']:fig.savefig(ROOT/('vOTU8476_NrdG_context.'+ext),dpi=300,facecolor='white')
plt.close(fig)
evidence=[{k:r[k] for k in ['gene_position','start_position','end_position','strandedness','virsorter_category','viral_hit','pfam_hits','vogdb_hits']} for r in hallmarks]
(ROOT/'hallmark_evidence.json').write_text(json.dumps(evidence,indent=2))
print('Saved compact context; VirSorter hallmark ORFs: 4, 15, 18.')
