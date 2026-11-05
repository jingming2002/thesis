# Combined AVG figure: panels A and B

A: Functional-module protein counts from the full AVG catalogue. All 32 modules are shown, ordered within class by full-catalogue protein count; panels have independent horizontal ranges. Unmapped or unresolved annotations are not shown.

B: Functional prevalence and lifestyle association among 6,061 predicted virulent and 3,076 predicted temperate vOTUs. All 32 modules are displayed, ordered within class by pooled vOTU carrier counts. OR is virulent relative to temperate; intervals are exact 95% confidence intervals; q values use BH correction across all 32 modules. Filled points indicate q < 0.05. Arrows show intervals extending beyond the axis; Inf denotes an infinite odds ratio.

The datasets, mapping qualifications, TA subset and statistical values are unchanged from the source figures. Both panels display all 32 mapped modules (12 AMG, 9 APG, 11 AReG). Unmapped annotations are not represented as modules. A uses protein counts from the full catalogue; B uses vOTU-level presence in the eligible lifestyle background. The component README files retain the full interpretation and mapping notes. This combined figure does not implement the newer 172/628 TA reclassification.

Both panels are drawn in R. Panel B was re-exported at the same physical width as panel A, retaining comparable text sizes. combine_panels.py composes the vector PDFs, adds A/B labels and replaces the separate figure titles with panel headings. It does not rasterize the PDF or change plotted data. The combined SVG and 300 dpi PNG are also included.

Reproduce each component using its R script with its directory as the command-line argument, then run python combine_panels.py (requires PyMuPDF). All component scripts and source data are bundled.
