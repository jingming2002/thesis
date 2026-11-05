# Supplementary tables

All tables use UTF-8 tab-separated text with one header row. They can be opened in spreadsheet software using tab as the delimiter. Import protein IDs, contig IDs and ORF partial codes as text.

## Table S1. Curated auxiliary viral gene catalogue and associated viral metadata

`Table_S1_AVG_catalogue.tsv` contains 11,725 candidate proteins with class assignments, original functions, annotation identifiers, viral-origin confidence, ORF coordinates and completeness, genome quality, final functional-module assignment, host predictions and lifestyle information. `Table_S1_vOTU_metadata.tsv` contains the 5,384 distinct AVG-positive vOTUs. `Table_S1_lifestyle_background.tsv` contains all 9,137 eligible vOTUs used in lifestyle tests, including genomes without retained AVGs. Confidence fields from different tools must not be treated as equivalent. Full, unabridged protein annotations are retained under `results/curation/`.

## Table S2. Functional modules and associations with predicted viral lifestyle

`Table_S2_module_definitions.tsv` reports 32 modules and their recorded assignment criteria, protein counts, carrier counts and numbers of distinct Function labels. `Table_S2_annotation_mapping.tsv` provides the original Function-to-module mapping. `Table_S2_protein_assignments.tsv` includes all 11,725 proteins, including unmapped and excluded-from-module records with their reasons. `Table_S2_lifestyle_tests.tsv` gives the 2×2 table numerators and denominators, carrier percentages, conditional odds ratios, exact 95% confidence intervals, two-sided Fisher P values and Benjamini–Hochberg q values across all 32 tests. Odds ratios compare virulent with temperate viruses. Infinite bounds are represented by `Inf`/`inf`.

Module counts cover 7,979 proteins. The remaining 3,746 proteins were retained in the catalogue but not assigned to these modules. Module labels summarize annotations and do not establish complete pathways or experimentally verified activities. The historical `functional_theme_L2` field is retained for traceability and is not the current module label.

## Table S3. Sequences included in the main phylogenetic analyses

The table contains 547 sequence records across Fhs, FolD, NrdG, FolE, GlyA, IscS, IscU-like and MetK. It reports tree-tip identifiers, sequence sources, accessions or genome IDs, available taxonomy, lengths, viral-origin confidence and ORF integrity. Fhs uses five viral candidates and 49 reference sequences (54 tips). An accession may occur in more than one analysis; rows represent analysis-specific sequence inclusion, not necessarily globally unique proteins. Server paths record provenance and are not portable download links.

## Table S4. Phylogenetic alignment and model-selection summary

Eight main analyses are reported with query scope, sequence-source counts, input and trimmed alignment lengths, BIC-selected substitution models and recorded MAFFT, trimAl and IQ-TREE settings. Inference was unrooted; midpoint rooting was applied subsequently for visualization. S3 and S4 exclude the expanded52 Fhs sensitivity tree, whose unchanged records remain in `results/phylogeny/server_audit/`.

## Table S5 (proposed). Representative genomic contexts of Fhs, FolD and NrdG

Three separate tables preserve the native annotations used to draw vOTU_1650 (Fhs), vOTU_17468 (FolD) and vOTU_8476 (NrdG). Native columns differ between sources and have deliberately not been forced into a misleading common annotation schema. Coordinates refer to the original contigs; drawing scripts handle orientation and relative coordinates. Hallmark labels for Fhs use the archived PHROG-based evidence, whereas the FolD/NrdG examples use positioned VirSorter categories. The contig-wide geNomad hallmark count is not substituted for individually located hallmark genes. See the context evidence JSON files and their README in `figures/figure4/composite/sources/`.

S1 is a proposed catalogue designation consistent with the existing S2–S4 references. S5 is additional supporting evidence and is not yet cited in the supplied manuscript.
