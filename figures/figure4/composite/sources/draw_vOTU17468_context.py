"""FolD cluster C12: compact, coordinate-grounded genome context.
Run: python draw_vOTU17468_context.py
Both panels are reversed to orient folD left-to-right.
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, ConnectionPatch

ROOT=Path(__file__).resolve().parent
with (ROOT/'vOTU17468_gene_context.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
    for k in ['start_position','end_position','strandedness','gene_position']:r[k]=int(r[k])
LENGTH=48558
target=next(r for r in rows if r['gene_position']==9)
assert (target['start_position'],target['end_position'],target['strandedness'])==(8606,9421,-1)
assert 'PF00763' in target['pfam_hits'] and 'PF02882' in target['pfam_hits']
# Convert one-based inclusive CDS positions to reversed zero-based half-open intervals.
for r in rows:
    r['a']=LENGTH-r['end_position'];r['b']=LENGTH-r['start_position']+1
    r['direction']=-r['strandedness']
origin=(target['a']+target['b'])/2
lo,hi=-5.15,5.15
hallmarks=[r for r in rows if r['virsorter_category'] in ('0','3')]
assert len(rows)==61 and len(hallmarks)==7
COL={'target':'#873E99','hallmark':'#4479B0','onecarbon':'#D79348','annotated':'#859BAA','other':'#DEDEDE'}
def category(r):
    n=r['gene_position']
    if n==9:return 'target'
    if r['virsorter_category'] in ('0','3'):return 'hallmark'
    if n==6:return 'onecarbon'
    if n in (3,5,11,12,15):return 'annotated'
    return 'other'

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(13,2.65),facecolor='white')
over=fig.add_axes([.13,.66,.82,.28]);over.set(xlim=(0,75),ylim=(-.5,1));over.axis('off')
over.add_patch(Rectangle((0,-.035),LENGTH/1000,.07,color='#D5D9DC',lw=0))
over.text(-1,0,'vOTU_17468',ha='right',va='center',fontsize=12,fontweight='bold')
over.text(LENGTH/1000+.6,0,'48.6 kb',ha='left',va='center',fontsize=8,color='#666666')
win=(origin/1000+lo,origin/1000+hi)
over.add_patch(Rectangle((win[0],-.10),win[1]-win[0],.20,facecolor='#F0E8F2',edgecolor=COL['target'],lw=.6,zorder=1))
for r in hallmarks+[target]:
    over.add_patch(Rectangle((r['a']/1000,-.075),(r['b']-r['a'])/1000,.15,facecolor=COL[category(r)],lw=0,zorder=2))
# Group leader lines avoid repeated overlapping labels for clustered hallmarks.
groups=[([51,52,53,55],5.0),([48],14.5),([33,35],25.2)]
for ids,x in groups:
    over.text(x,.62,'hallmark',ha='center',va='bottom',fontsize=8,color=COL['hallmark'])
    for r in hallmarks:
        if r['gene_position'] in ids:
            mid=(r['a']+r['b'])/2000
            over.plot([mid,mid,x],[.08,.25,.60],color=COL['hallmark'],lw=.55)
over.text(origin/1000,.17,'folD',ha='center',va='bottom',fontstyle='italic',fontweight='bold',color=COL['target'],fontsize=9)
over.plot([62,72],[-.19,-.19],lw=1,color='#666666')
over.text(67,-.29,'10 kb',ha='center',va='top',fontsize=8,color='#555555')

ax=fig.add_axes([.13,.25,.82,.38]);ax.set(xlim=(lo,hi),ylim=(-.28,.67));ax.set_yticks([])
for s in ['top','right','left']:ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#B9BEC2');ax.spines['bottom'].set_linewidth(.65)
ax.tick_params(axis='x',labelsize=8,length=2,color='#AAAAAA')
ax.set_xticks([-4,-2,0,2,4]);ax.set_xlabel('Position relative to folD midpoint (kb)',fontsize=9,labelpad=5)
ax.axvline(0,ls=':',color='#D6D1D9',lw=.6,zorder=0)
for r in rows:
    a=(r['a']-origin)/1000;b=(r['b']-origin)/1000
    if b<=lo or a>=hi:continue
    s=max(a,lo);e=min(b,hi);h=.072;head=min(.13,.45*(e-s))
    clipped=(r['direction']==1 and b>hi) or (r['direction']==-1 and a<lo)
    if clipped:pts=[(s,-h),(e,-h),(e,h),(s,h)]
    elif r['direction']==1:pts=[(s,-h),(e-head,-h),(e,0),(e-head,h),(s,h)]
    else:pts=[(e,-h),(s+head,-h),(s,0),(s+head,h),(e,h)]
    ax.add_patch(Polygon(pts,facecolor=COL[category(r)],edgecolor='#767676',lw=.5,zorder=3))
labels={3:'Exonuclease',5:'DNA polymerase',6:'ThyX',9:'folD',11:'DNA helicase',15:'DNA ligase'}
for r in rows:
    if r['gene_position'] not in labels:continue
    a=(r['a']-origin)/1000;b=(r['b']-origin)/1000
    x=(max(lo,a)+min(hi,b))/2
    if r['gene_position']==3:x=hi-.38
    ax.text(x,.14,labels[r['gene_position']],ha='center',va='bottom',fontsize=9,
        color=COL['target'] if r['gene_position']==9 else '#333333',
        fontweight='bold' if r['gene_position']==9 else 'normal',fontstyle='italic' if r['gene_position']==9 else 'normal')
for x,xx in zip(win,[lo,hi]):
    fig.add_artist(ConnectionPatch((x,-.105),(xx,.63),coordsA=over.transData,coordsB=ax.transData,color='#D3CDD7',lw=.65,zorder=0))
for ext in ['pdf','svg','png']:fig.savefig(ROOT/('vOTU17468_FolD_context.'+ext),dpi=300,facecolor='white')
plt.close(fig)
print('Saved vOTU_17468: 61 CDS, FolD ORF 9, 7 positioned VirSorter hallmarks.')
