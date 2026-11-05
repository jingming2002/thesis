# FolD genomic context: vOTU_17468

Compact layout matching the Fhs and NrdG context figures. The gene order, lengths and directions match the supplied FolD cluster C12 screenshot. Matching also uses the genome-level virus score (0.9811, rounded to 0.981) and geNomad hallmark count (8).

- Contig: vOTU_17468|k79_636348@FCR.HF3_300|72
- Genome length: 48,558 bp; CheckV completeness: 100%; predicted lifestyle: virulent.
- FolD: DRAM-v ORF 9, positions 8,606–9,421 (one-based inclusive), reverse strand; PF02882.26 and PF00763.29. Viral database hit: tetrahydrofolate dehydrogenase/cyclohydrase. VOG hit: bifunctional protein FolD.
- Both panels are reversed to orient folD left-to-right. Local coordinates are relative to the folD midpoint; the viewport spans −5.15 to +5.15 kb. CDS intersecting the viewport boundary are clipped and drawn without an artificial arrowhead at the clipped end.
- Purple: folD; blue: located hallmark CDS; ochre: ThyX; blue-grey: selected annotated DNA-associated proteins; light grey: other CDS. Light grey does not necessarily mean unannotated.
- The full contig has 61 DRAM-v CDS. Seven positioned hallmark CDS (ORFs 33, 35, 48, 51, 52, 53, 55) are identified by VirSorter categories 0/3 in the DRAM-v annotations. The geNomad summary reports eight hallmarks, but corresponding per-gene coordinates were not available here. These distinct annotation sets are not treated as interchangeable. All seven available positioned hallmarks are plotted, with grouped labels to avoid crowding.

Sources: context_source/masters/annotations.tsv and life_source/extracted/AVG_lifestyle_comparison/inputs/eligible_background.tsv. The complete contig annotations and hallmark evidence are included in this bundle.

Suggested caption: Genomic context of folD in vOTU_17468. The upper panel shows the complete genome, with folD in purple and VirSorter hallmark genes in blue. The lower panel expands the folD neighbourhood; arrows indicate CDS positions and orientations. Both panels are oriented with folD transcribed from left to right.

Reproduce with Python and matplotlib: `python draw_vOTU17468_context.py`.
