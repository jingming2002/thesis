"""Arrange supplied vector trees with matching gene-family genomic examples.
Requires PyMuPDF. Run from any directory with the supplied sources folder.
"""
from pathlib import Path
import fitz

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'sources'
tree=fitz.open(SRC/'AVG3_midpoint_combined.pdf')
assert len(tree)==1 and abs(tree[0].rect.width-1944)<1
W=1440; ROW=550; TOP=18; BOTTOM=52
out=fitz.open();page=out.new_page(width=W,height=TOP+3*ROW+BOTTOM)
specs=[('Fhs',0,648,'vOTU1650_Fhs_context'),
       ('FolD',648,1288,'vOTU17468_FolD_context'),
       ('NrdG',1288,1944,'vOTU8476_NrdG_context')]
for i,(gene,x0,x1,stem) in enumerate(specs):
    y=TOP+i*ROW
    page.insert_text((18,y+20),f'{chr(65+2*i)}  {gene}',fontname='hebo',fontsize=18,color=(.12,.16,.19))
    page.show_pdf_page(fitz.Rect(0,y+34,570,y+ROW-6),tree,0,
        clip=fitz.Rect(x0,44,x1,648),keep_proportion=True)
    # Keep the overview and its expanded neighbourhood together.
    with fitz.open(SRC/f'{stem}.pdf') as context:
        page.show_pdf_page(fitz.Rect(568,y+138,1425,y+381),context,0,
            clip=fitz.Rect(15,8,904,177),keep_proportion=True)
    page.insert_text((580,y+20),chr(66+2*i),fontname='hebo',fontsize=18,color=(.12,.16,.19))

# Reuse the original phylogenetic legend unchanged, with a readable scale.
page.show_pdf_page(fitz.Rect(180,TOP+3*ROW+8,1320,TOP+3*ROW+34),tree,0,
    clip=fitz.Rect(335,662,1665,684),keep_proportion=True)
out.set_metadata({'title':'Phylogenetic relationships and genomic contexts of Fhs, FolD and NrdG',
 'subject':'A–B Fhs; C–D FolD; E–F NrdG. Genomic examples represent gene families, not necessarily displayed tree tips.'})
out.save(ROOT/'AVG_phylogeny_and_genomic_context.pdf',garbage=4,deflate=True)
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(ROOT/'AVG_phylogeny_and_genomic_context.png')
# Convert the assembled vector PDF to a portable SVG; glyph paths preserve fonts.
(ROOT/'AVG_phylogeny_and_genomic_context.svg').write_text(page.get_svg_image(text_as_path=True))
page.get_pixmap(matrix=fitz.Matrix(1500/W,1500/W),alpha=False).save(ROOT/'preview.png')
assert len(out)==1
print('Created one-page vector PDF, SVG and 300-dpi PNG.')
