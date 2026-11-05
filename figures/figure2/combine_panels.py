from pathlib import Path
import fitz

root = Path(__file__).resolve().parent
a = fitz.open(root / 'components/panel_A/AVG_functional_protein_counts.pdf')
b = fitz.open(root / 'components/panel_B/AVG_repertoire_redesign.pdf')
width = 1224  # 17 inches; both component PDFs share the same physical width.
clip_a = fitz.Rect(0, 45, width, a[0].rect.height - 9)
clip_b = fitz.Rect(0, 39, width, b[0].rect.height - 20)
header = 39
gap = 14
height = 2 * header + clip_a.height + clip_b.height + gap + 12
doc = fitz.open()
page = doc.new_page(width=width, height=height)
ink = (0.13, 0.19, 0.24)
page.insert_text((18, 26), 'A', fontsize=22, fontname='hebo', color=ink)
page.insert_text((49, 25), 'Functional composition of auxiliary viral genes', fontsize=17, fontname='hebo', color=ink)
page.show_pdf_page(fitz.Rect(0, header, width, header + clip_a.height), a, 0, clip=clip_a)
by = header + clip_a.height + gap
page.insert_text((18, by + 26), 'B', fontsize=22, fontname='hebo', color=ink)
page.insert_text((49, by + 25), 'Functional repertoires across viral lifestyles', fontsize=17, fontname='hebo', color=ink)
page.show_pdf_page(fitz.Rect(0, by + header, width, by + header + clip_b.height), b, 0, clip=clip_b)
doc.set_metadata({'title': 'Pig gut AVG functional composition and lifestyle associations', 'subject': 'A: protein counts; B: carriage frequencies and lifestyle odds ratios', 'creator': 'R-generated plots composed with PyMuPDF'})
doc.save(root / 'AVG_functional_composition_and_lifestyle_AB.pdf', garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(1600/width,1600/width),alpha=False).save(root / 'preview.png')
page.get_pixmap(dpi=300,alpha=False).save(root / 'AVG_functional_composition_and_lifestyle_AB.png')
(root / 'AVG_functional_composition_and_lifestyle_AB.svg').write_text(page.get_svg_image(), encoding='utf-8')
assert len(page.get_images()) == 0, 'Expected vector-only PDF panels'
print(f'Combined vector figure: {width/72:.2f} x {height/72:.2f} inches')
