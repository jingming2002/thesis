# Functional modules and protein counts by AVG class

The data table contains all 32 named modules from the revised repertoire mapping: 12 AMG, 9 APG and 11 AReG modules. This version displays all 32 modules. Coverage counts below refer to all 32 modules. Rows are sorted within each class by decreasing protein count. Bar labels are exact protein counts. Horizontal axes have different ranges to show within-class patterns clearly; compare printed values for comparisons across classes. All bars use the same color, without a redundant color legend.

Counts use the full 11,725-protein catalogue, including candidates without an eligible lifestyle prediction. Each distinct Protein identifier is counted once in its assigned module. These are protein records, not sequence-deduplicated clusters, sequencing abundance, or unique vOTU counts. The count table also supplies distinct vOTU counts for reference.

The displayed named modules cover 2,177 of 3,028 AMG proteins, 1,631 of 1,808 APG proteins and 4,171 of 6,889 AReG proteins. Other/unresolved entries are recorded in class_count_summary.tsv and the full input audit; they are not plotted as an additional biological pathway. Thus the displayed bars do not sum to the full class totals.

Module-level counts combine all constituent annotations and should not be substituted for individual named-gene counts in Supplementary Table S1. For example, cobalamin biosynthesis-related includes more than CobS. Functional labels summarize annotation associations, not demonstrated complete pathways or experimentally validated biological effects.

This figure inherits the revised module mapping from AVG_repertoire_redesign: 76 nonspecific or unsupported toxin annotations are omitted from the TA module, and seven translation repressor annotations are omitted from transcriptional repression. The full audit retains all records and revision reasons. The TA count of 543 is an annotation-based candidate subset, not the newer 172-toxin/628-antitoxin classification, which still requires the corresponding per-protein assignment table before integration. No claim about orphan status or pairing is made here. No original AVG category labels or catalogue totals were changed.

## Reproduce

Rscript plot_protein_counts.R /absolute/path/to/AVG_functional_protein_counts

Tested with R 4.3.3, ggplot2 3.4.4 and grid 4.3.3. Optional SVG export uses svglite 2.1.3. The exact runtime record is in sessionInfo.txt. No significance testing is used for this descriptive count figure.

Files include vector PDF/SVG, 300 dpi PNG, preview PNG, R code, full protein-level mapping, module counts and class coverage summary.
