# Main outputs by tool

| Tool or stage | Included main results | Scope and interpretation |
| --- | --- | --- |
| CheckAMG | `results/checkamg/final_results_strict_viral.tsv.gz` | 819,843 protein rows in strict viral regions; 13,928 carry metabolic, physiological or regulatory classifications. This is not the initial 40,314-candidate output. |
| ORF and candidate curation | `results/curation/AVG_final_checkv90_complete.tsv.gz`, `AVG_final_filter_audit.tsv.gz`, ORF/quality audit tables | Full 11,725-protein catalogue and 13,640-candidate audit universe; source flags and exclusions retained. |
| CheckV | `results/checkv/quality_summary_from_merged_source.tsv` | Projection of CheckV fields from the 48,299-row merged quality/host source, not a recovered native run directory. No new contamination analysis was performed. |
| geNomad | `results/genomad/virus_summary_from_merged_source.tsv` | Projection of genome scores, hallmark counts, topology and taxonomy. Missing values remain missing. Native gene/marker outputs are not available. |
| PhaTYP | `results/phatyp/phatyp_raw.csv`, `eligible_background.tsv`, `PhaTYP_background_audit.tsv` | Raw predictions for 9,927 records; 9,137 eligible backgrounds. Historical protein-level merged output is preserved separately from final S1/S2 assignments. |
| PhaGCN3 | `results/phagcn3/final_prediction.csv`, `PhaGCN3_parsed_by_rank.tsv`, `final_AVG_vOTU_taxonomy.tsv` | Original predictions, 5,982 parsed prediction rows and assignments for 5,384 final AVG-positive vOTUs. Raw prediction CSV uses variable-width taxonomic paths; do not impose a fixed five-column parser. |
| DRAM-v | `results/dramv/annotations.tsv.gz`, `amg_summary.tsv` | Original annotations and 4,858 summary rows. These are additional annotation evidence, not an alternative definition of the final CheckAMG catalogue. |
| MAFFT / trimAl | `results/phylogeny/{family}/` | Input and aligned/trimmed FASTA for seven families. The current 54-tip Fhs alignments are missing. |
| IQ-TREE / ModelFinder | `results/phylogeny/{family}/`, S4 | Current trees for all eight main families; reports and logs for seven. Fhs high5 model settings are recorded in the server audit table. |
| BLASTp / HMMER / CD-HIT | `results/phylogeny/server_audit/` | Sequence manifests and verified methods record are included; complete current raw search hits and cluster files were not recovered. |
| R statistics | S2 lifestyle tests and `scripts/recalculate_lifestyle_tests.R` | Reconstructed from final protein assignments and the eligible background. Two-sided Fisher tests, conditional ORs, exact 95% CIs and BH adjustment over all 32 modules. |
| Genomic context | `results/genome_context/`, Figure 4 source tables and S5 | 642-vOTU NrdG co-annotation screen plus detailed examples. Absence of an nrdD annotation is not a sequence-level proof of gene absence. |

`tool_file_index.tsv` provides file sizes and line counts. Counts of lines in raw outputs are not automatically counts of distinct vOTUs, unique proteins or initial-screen candidates.

The full 48,299-row source in `data/PVD_quality_host_annotations.tsv.gz` preserves the connection between CheckV, geNomad and host annotations. The derived per-tool projections retain source row order and IDs but must not be described as byte-identical original tool exports.
