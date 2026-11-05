# AVG repertoire redesign: prevalence and lifestyle association

Design preview prepared 26 October 2026 from the existing AVG_lifestyle_comparison bundle (Library id libfile_2886fccb5c048191a581f8e3ab907864).

## Files and reproduction

- AVG_repertoire_redesign.pdf / .svg: vector figure.
- AVG_repertoire_redesign.png: 300 dpi figure; preview.png: smaller display copy.
- plot_repertoire.R: complete analysis and drawing, using R stats and grid; svglite is optional.
- inputs/: source protein-level module mapping and eligible viral background.
- revised_mapping_audit.tsv, mapping_changes.tsv: annotation-level revisions.
- unique_module_carriers.tsv: distinct module/vOTU associations.
- all_module_tests.tsv: all 32 tested and displayed modules.
- plotted_data.tsv: exact values and labels for the 32 displayed rows.
- caption_and_methods.txt, sessionInfo.txt, verification.json.

Run: Rscript plot_repertoire.R /absolute/path/to/AVG_repertoire_redesign

## What is displayed

Each class shows all mapped modules (12 AMG, 9 APG, 11 AReG). Within each class, rows are ordered by that pooled count. A viral genome can carry multiple modules, so percentages are not compositional and do not sum to 100%.

The prevalence denominator is all eligible vOTUs in each predicted lifestyle, including genomes without a mapped AVG: 6,061 virulent and 3,076 temperate. Unlike the earlier host-order bubble figure, this analysis does not require a host prediction and is not restricted to its 3,896 displayed host-assigned vOTUs.

Left: percentage of vOTUs carrying at least one gene in the module. Orange = virulent; teal = temperate. The small vertical offset separates closely overlapping points. The light connecting line is a visual aid, not a paired design.

Right: conditional maximum-likelihood odds ratio, virulent relative to temperate, and two-sided exact 95% confidence interval from R fisher.test. OR > 1 indicates higher carriage odds among predicted virulent genomes. Color indicates the direction of the estimate; filled markers indicate BH q < 0.05, open markers q >= 0.05. Intervals are unadjusted 95% intervals, not simultaneous/BH-adjusted intervals. The q column reports the multiple-testing-adjusted result.

Each 2x2 table contains lifestyle by module carriage/non-carriage. BH correction includes all 32 mapped modules, including low-count modules. This differs from the older lifestyle table, which excluded modules with fewer than 20 pooled carriers, and from the host-order figure's 64 omnibus tests. Tests describe unadjusted catalogue-level associations; host taxonomy, genome length and viral relatedness are not controlled for here.

## Classification revisions and scope

This is a new visualization preview, not a replacement for the previous figures or Supplementary Table S1. Other module definitions are inherited from the existing mapping; this is not an independent reannotation or a complete re-evaluation of the boundary between auxiliary and core viral genes.

76 proteins were omitted from the old toxin-antitoxin module because their labels were nonspecific toxin calls or did not provide explicit TA-family support. These include generic toxin, YdaS/YdaT, hemolysin/leukocidin, toxin transport/secretion and papain-fold toxin labels. Their omission means insufficient support for this particular module, not proven absence of TA activity. The retained TA module contains 543 candidate proteins across the full annotation set; 437 eligible lifestyle-assigned vOTUs carry at least one retained candidate (176 virulent; 261 temperate). These calls are annotation-based candidates, not validated paired operons or an orphan-antitoxin analysis.

The user's newer counts of 172 toxins and 628 antitoxins are NOT used: their underlying protein/vOTU-level assignment table was not available in these inputs. Reconcile that table before treating the present TA row as final. Antitoxin-related entries originally classified as AReG are not moved into APG in this preview.

Seven proteins annotated as translation repressors were omitted from the transcriptional repression module. All excluded proteins remain present in the audit file and underlying catalogue. Neither change alters original AVG category labels or total catalogue size.

Cobalamin, sporulation, CRISPR-Cas and other module labels denote associated annotations, not evidence for complete pathways or experimentally validated phenotypes. Gene-level top-ten identities and counts remain in Supplementary Table S1; the present module counts are different quantities.
