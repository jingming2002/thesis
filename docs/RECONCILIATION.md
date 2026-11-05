# Manuscript and source-data reconciliation

The manuscript was not edited. The following items affect reproducibility or interpretation and should be reconciled before release.

| Item | Manuscript | Recovered source evidence | Action |
| --- | --- | --- | --- |
| Initial screen | 40,314 candidates from 18,220 vOTUs | Figure 1 contains the class/confidence summaries; available strict-region output has 819,843 total protein rows and 13,928 AVG-classified rows, at a later stage | Obtain initial protein-level output and the exact 18,220 input IDs |
| High-quality genomes | 5,411 selected as high-quality | Quality audit: 5,411 of 6,074 candidate-bearing vOTUs pass all four criteria. The broad 48,299-row merged table has 14,740 rows satisfying completeness >=90, length >10 kb, score >0.90 and >=1 hallmark, without imposing the 18,220-input membership or AAI-only restriction | State the selection population explicitly; 14,740 is not a proposed replacement count |
| Completeness wording | AAI-based completeness >=90% | The merged quality table also contains terminal-repeat and other completeness methods | Recover the exact input filter and clarify which completeness field/method was applied |
| BLAST search E-value | Uniform 1e-5 | Server audit: 0.001 for Fhs and five rebuilt families; 1e-5 for FolD/NrdG. Recorded hit filter: bit score >=50 and E-value <=0.001 | Reconcile with run logs; see `results/phylogeny/server_audit/Methods_verified.md` |
| NCBI reference selection | One best hit per viral query | Seven families used one best hit; current Fhs retained up to five accessions per each of five queries, producing 23 distinct NCBI tips | Describe Fhs exception |
| Fhs source version | Five viral candidates in main discussion | Main current tree: 54 tips, LG+R5; expanded52: 101 tips, separate sensitivity set. Historical 199-tip tree is not the current analysis | S3/S4 main tables include high5 only; do not substitute historical alignments |
| Full-catalogue high confidence | Abstract calls all 11,725 high-confidence candidates | CheckAMG high counts are 2,643 AMG, 1,589 APG and 6,512 AReG; medium and low records remain | Distinguish curated candidates from CheckAMG high-confidence assignments |
| Figure 2 terminology | “mild”, “bilateral”, “precise” | Categories and tests are temperate, two-sided Fisher, exact confidence intervals | Correct caption wording in manuscript |
| Plotting software | ggtree, ape and patchwork are listed | Final supplied phylogeny plot uses base R/grid and ape; contexts use matplotlib, and composites use PyMuPDF. Actual versions are archived in sessionInfo files and validation environment | Update visualization Methods to match final scripts |
| Figure 3 | Styled pathway PNG | Exact embedded PNG preserved; available R script is an editable reconstruction with different styling | Treat PNG as manuscript version; do not claim pixel-identical R output |
| RNA/TA annotation exclusions | Final 32-module catalogue | 83 records were omitted from prior modules: 76 selected toxin labels and 7 translation-repressor records. These remain in the overall catalogue | S2 preserves both final assignments and reasons; no fresh annotation decisions were made here |

Figure 2A and B use the same final protein mapping. The full catalogue has 11,725 proteins, while the eligible lifestyle subset has 10,604 proteins in 4,860 vOTUs. The statistical background has 9,137 vOTUs and includes zero-carrier genomes. These are different, valid denominators.

Archived Fisher statistics were independently recomputed from the protein mapping. R-to-R comparisons use the original `fisher.test` implementation; the SciPy cross-check permits relative tolerance 0.001 for conditional odds ratios and confidence limits because the implementations use different numerical root-finding tolerances. Counts and memberships require exact equality.
