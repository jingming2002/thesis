# Integrated phylogeny and genomic-context figure

Panels A–B: Fhs phylogeny and vOTU_1650 context.
Panels C–D: FolD phylogeny and vOTU_17468 context.
Panels E–F: NrdG phylogeny and vOTU_8476 context.

The supplied AVG3_midpoint_combined.pdf is the authoritative tree source. No topology, rooting, branch support, labels or clade collapsing was changed. Each row pairs a gene-family phylogeny with a representative genomic context. vOTU_17468 and vOTU_8476 are not labelled tips in the supplied trees; pairing therefore indicates the gene family, not a direct tip-to-genome mapping. No linking lines are drawn.

The original phylogenetic colour legend is retained at the bottom. In the genomic contexts, purple denotes the target gene and blue denotes hallmark annotations. Genomic scales differ and are shown individually. Panel D is oriented with folD left-to-right. Context panels preserve their previously verified annotations and coordinates.

Fhs hallmark labels derive from the original PHROG-based source classification; FolD and NrdG hallmark labels derive from DRAM-v/VirSorter categories 0/3. The normalized source_hallmark_flag in the Fhs table is a source classification flag, not a VirSorter score.

PDF retains vector artwork and selectable text. SVG retains vector shapes with font glyphs converted to paths for portability. PNG is 300 dpi. The bundle contains the supplied tree PDF, context source figures, context data and plotting scripts, and the assembly script. The supplied PDF does not contain the original tree-building input files or R plotting script, so those are not included here.

Reproduce the assembly with PyMuPDF: `python combine_phylogeny_context.py`.
