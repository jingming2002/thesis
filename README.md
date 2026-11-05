# Auxiliary viral genes in the pig gut virome

Data, analysis outputs and figure scripts accompanying **Auxiliary Viral Genes in the Pig Gut Virome: Functional Diversity and Characterization of One-Carbon Metabolism Genes**. Compiled against the supplied manuscript on 5 November 2026.

The final catalogue contains **11,725 candidate proteins in 5,384 vOTUs**: 3,028 AMGs, 1,808 APGs and 6,889 AReGs. The lifestyle-comparison background contains 9,137 vOTUs, including 4,860 final AVG-positive vOTUs. Functional modules cover 7,979 proteins; unmapped proteins remain in the catalogue.

## Contents

| Directory | Contents |
| --- | --- |
| `figures/figure1/` | Overview plot, R script and panel input tables |
| `figures/figure2/` | All 32 modules, protein counts, carriage statistics, R scripts and Python composition script |
| `figures/figure3/` | Exact manuscript PNG and an editable R reconstruction with labels/edges |
| `figures/figure4/` | Three phylogenies, tip metadata, collapse audits, context annotations and plotting/composition scripts |
| `figures/supplementary_trees/` | Complete and collapsed views for FolE, GlyA, IscS, IscU-like and MetK |
| `figures/manuscript_embedded/` | Unmodified EMF objects extracted from the Word document for provenance |
| `supplementary/` | Tables S1–S4 and proposed S5, in UTF-8 TSV format |
| `tables/` | Table 1 extracted from the manuscript |
| `results/` | CheckAMG, curation, CheckV, geNomad, PhaTYP, PhaGCN3, DRAM-v and phylogenetic results |
| `data/` | PVD quality/host annotations and reconciled vOTU indexes |
| `scripts/` | Figure reproduction, statistical recalculation and validation |
| `docs/` | Figure/table index, tool output guide, provenance, version discrepancies and missing server files |
| `validation/` | Checks performed during packaging and figure smoke-test results |

## Reproduce the figures

The scripts run from the repository after installation of the dependencies in `environment/`. They write to a separate directory and preserve the archived figures.

```bash
python scripts/validate_repository.py
Rscript scripts/recalculate_lifestyle_tests.R
python scripts/reproduce_supplementary.py
python scripts/reproduce_figures.py --output reproduced
```

To render selected figures only:

```bash
python scripts/reproduce_figures.py --figures 1 2 --output reproduced
```

All four figure workflows were executed successfully during packaging. The exact vector PDFs embedded inside the manuscript's EMF objects were recovered as `figures/Figure1_manuscript.pdf`, `Figure2_manuscript.pdf` and `Figure4_manuscript.pdf`. These were rendered and compared against regenerated outputs; see `validation/figure_comparison.json`. Figure 3 is a conceptual diagram: its R output reproduces the pathway layout but is not the exact styled PNG embedded in the manuscript. That PNG is supplied as `figures/Figure3_manuscript.png`.

## Supplementary tables

- **S1:** final protein catalogue, AVG-positive vOTU metadata and the complete lifestyle background.
- **S2:** module definitions, annotation mapping, protein assignments and all 32 lifestyle tests.
- **S3:** 547 sequence records across the eight main phylogenetic analyses, using the five-candidate Fhs analysis.
- **S4:** model selection and alignment summaries for the same eight analyses.
- **S5 (proposed):** native gene-context annotations for vOTU_1650, vOTU_17468 and vOTU_8476. Add an S5 citation if used in the manuscript.

TSV files are the editable, version-control-friendly supplementary tables. Multiple files with the same S number represent sections/sheets of that table. See `supplementary/README.md` for captions and field definitions.

## Scope and outstanding source files

This is a figure-and-results release, **not yet a complete rerun from the original viral FASTA files**. The current Fhs tree and tip manifest are present, but its original input FASTA, alignments and IQ-TREE logs were not supplied. Seven other families have matching input sequences, alignments and logs. The full protein-level source for the initial 40,314 candidates is also missing; the strict-region CheckAMG table is a later stage and must not be substituted for that initial screen.

`docs/RECONCILIATION.md` records manuscript/data differences without silently changing the manuscript. `docs/SERVER_EXPORT_REQUEST.txt` gives the exact remaining collection task for the server. The source Word manuscript and supervisor comments are not included in this analysis repository.

## Data conventions

Join proteins by full `Protein` and genomes by full `Contig`; short vOTU labels are for display. Counts are protein records or distinct carriers, not metagenomic abundance. Gene classes and modules can co-occur in a genome. `NA` denotes missing values; a missing annotation is not proof of biological absence. Preserve ORF `partial` fields as text. High-confidence viral origin, predicted ORF completeness and high-confidence lifestyle are separate concepts.

Large annotation tables are losslessly gzip-compressed. Original software outputs are distinguished from projections of merged tables. Existing `functional_theme_L2` annotations are historical; use `Display_module` in S2 for the final Figure 2 mapping. The original Fhs expanded52 analysis is retained only under the server audit directory and is not part of S3/S4 main analyses.

## Integrity and provenance

`FILE_MANIFEST.tsv` gives byte sizes and SHA-256 hashes. `docs/source_file_provenance.tsv` records recovered source files and transformations. `docs/manuscript_source.tsv` identifies the exact manuscript version. The independent validation script checks catalogue totals, module counts, Fisher tests, tree tips, alignment sizes and file integrity. Generated files under `reproduced/` are ignored by Git.

## Upload to your Git repository

From this directory, use your chosen remote URL:

```bash
git init
git add .
git commit -m "Add pig gut AVG data, figures and analysis outputs"
git branch -M main
git remote add origin YOUR_REPOSITORY_URL
git push -u origin main
```

No remote repository was created and nothing was published during preparation. A reuse licence has not been assigned; choose one with the project owners before granting reuse permissions.
