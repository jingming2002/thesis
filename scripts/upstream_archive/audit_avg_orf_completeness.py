#!/usr/bin/env python3
"""Audit AVG ORF completeness without changing the candidate universe."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


CATEGORY_ORDER = {"AMG": 0, "APG": 1, "AReG": 2}
HEADER_RE = re.compile(r"^(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#\s+(.*)$")
KO_RE = re.compile(r"K\d{5}")


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_attrs(text: str):
    out = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def iter_fasta(path: Path):
    header = None
    seq = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].rstrip("\r\n")
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)


def load_candidate_proteins(path: Path, wanted: set[str]):
    records = {}
    occurrences = Counter()
    malformed = []
    for header, sequence in iter_fasta(path):
        protein = header.split(None, 1)[0]
        if protein not in wanted:
            continue
        occurrences[protein] += 1
        match = HEADER_RE.match(header)
        if not match:
            malformed.append(protein)
            continue
        pid, start, end, strand, attr_text = match.groups()
        attrs = parse_attrs(attr_text)
        clean = sequence[:-1] if sequence.endswith("*") else sequence
        internal_stops = clean.count("*")
        nonstandard = sorted(set(clean.upper()) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO*"))
        nt_len = int(end) - int(start) + 1
        aa_len = len(clean)
        expected = {nt_len // 3, max(0, nt_len // 3 - 1)}
        records[protein] = {
            "fasta_header": header,
            "orf_start": int(start),
            "orf_end": int(end),
            "strand": "+" if strand == "1" else "-",
            "strand_numeric": strand,
            "partial": attrs.get("partial", ""),
            "partial_5prime": attrs.get("partial", "")[:1] if len(attrs.get("partial", "")) == 2 else "",
            "partial_3prime": attrs.get("partial", "")[1:] if len(attrs.get("partial", "")) == 2 else "",
            "start_type": attrs.get("start_type", ""),
            "rbs_motif": attrs.get("rbs_motif", ""),
            "rbs_spacer": attrs.get("rbs_spacer", ""),
            "gc_cont": attrs.get("gc_cont", ""),
            "genetic_code_header": attrs.get("genetic_code", ""),
            "protein_length_aa": aa_len,
            "terminal_stop_present": str(sequence.endswith("*")).lower(),
            "internal_stop_count": internal_stops,
            "nonstandard_aa": "".join(nonstandard),
            "translation_length_consistent": str(aa_len in expected).lower(),
        }
    return records, occurrences, malformed


def load_target_contig_lengths(path: Path, wanted: set[str]):
    lengths = {}
    occurrences = Counter()
    for header, sequence in iter_fasta(path):
        contig = header.split(None, 1)[0]
        if contig in wanted:
            occurrences[contig] += 1
            lengths[contig] = len(sequence)
    return lengths, occurrences


def load_gene_index(path: Path, wanted: set[str]):
    records = {}
    occurrences = Counter()
    for row in read_tsv(path):
        protein = row.get("protein", "")
        if protein in wanted:
            occurrences[protein] += 1
            records[protein] = row
    return records, occurrences


def load_by_key(path: Path, key: str, wanted: set[str]):
    records = {}
    occurrences = Counter()
    for row in read_tsv(path):
        value = row.get(key, "")
        if value in wanted:
            occurrences[value] += 1
            records[value] = row
    return records, occurrences


def numeric(value, kind=float):
    try:
        return kind(float(value))
    except (TypeError, ValueError):
        return None


def bool_text(value):
    if value is None:
        return ""
    return str(bool(value)).lower()


def load_theme_maps(category_dir: Path):
    maps = {}
    for category, stem in (("AMG", "AMGs"), ("APG", "APGs"), ("AReG", "AReGs")):
        mapping = {}
        for row in read_tsv(category_dir / f"{stem}_categorized.tsv"):
            mapping[(row["db"], row["id"])] = row
        maps[category] = mapping
    return maps


def assign_theme(row, mapping):
    top_db = row.get("top_hit_db", "")
    top_id = row.get("top_hit_hmm_id", "")
    if (top_db, top_id) in mapping:
        hit = mapping[(top_db, top_id)]
        return hit, "top_hit"
    for db, field in (
        ("KEGG", "KEGG_hmm_id"), ("FOAM", "FOAM_hmm_id"),
        ("Pfam", "Pfam_hmm_id"), ("dbCAN", "dbCAN_hmm_id"),
        ("METABOLIC", "METABOLIC_hmm_id"), ("CAMPER", "CAMPER_hmm_id"),
        ("PHROG", "PHROG_hmm_id"),
    ):
        value = row.get(field, "")
        if (db, value) in mapping:
            return mapping[(db, value)], f"{db}_fallback"
    return {"category_L1": "Unmapped", "category_L2": "", "category_L3": "", "name": ""}, "unmapped"


def count_checkamg_start(files):
    memberships = defaultdict(set)
    contigs = {}
    for category, path in files:
        for row in read_tsv(path):
            protein = row["Protein"]
            memberships[protein].add(category)
            contigs[protein] = row["Contig"]
    chosen = {p: min(cats, key=CATEGORY_ORDER.get) for p, cats in memberships.items()}
    output = {}
    for category in CATEGORY_ORDER:
        proteins = {p for p, cat in chosen.items() if cat == category}
        output[category] = (len(proteins), len({contigs[p] for p in proteins}))
    output["TOTAL"] = (len(chosen), len(set(contigs.values())))
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--contig-fasta", required=True, type=Path)
    parser.add_argument("--gene-index", required=True, type=Path)
    parser.add_argument("--circular-table", required=True, type=Path)
    parser.add_argument("--quality-table", required=True, type=Path)
    parser.add_argument("--category-dir", required=True, type=Path)
    parser.add_argument("--checkamg-metabolic", required=True, type=Path)
    parser.add_argument("--checkamg-physiology", required=True, type=Path)
    parser.add_argument("--checkamg-regulation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    candidates = list(read_tsv(args.candidates))
    proteins = [row["Protein"] for row in candidates]
    contigs = [row["Contig"] for row in candidates]
    wanted_proteins, wanted_contigs = set(proteins), set(contigs)
    candidate_occ = Counter(proteins)

    faa, faa_occ, malformed = load_candidate_proteins(args.protein_fasta, wanted_proteins)
    index, index_occ = load_gene_index(args.gene_index, wanted_proteins)
    contig_lengths, contig_occ = load_target_contig_lengths(args.contig_fasta, wanted_contigs)
    circular, circular_occ = load_by_key(args.circular_table, "contig", wanted_contigs)
    quality, quality_occ = load_by_key(args.quality_table, "seq_name", wanted_contigs)
    themes = load_theme_maps(args.category_dir)

    audit_rows = []
    for source in candidates:
        protein, contig = source["Protein"], source["Contig"]
        protein_rec = faa.get(protein, {})
        idx = index.get(protein, {})
        q = quality.get(contig, {})
        circ = circular.get(contig, {})
        length = contig_lengths.get(contig)
        start = protein_rec.get("orf_start")
        end = protein_rec.get("orf_end")
        strand = protein_rec.get("strand", "")
        partial = protein_rec.get("partial", "")
        left = start - 1 if start is not None else None
        right = length - end if length is not None and end is not None else None
        near_edge = left is not None and right is not None and min(left, right) <= 100

        topology = q.get("topology", "")
        repeat_prediction = circ.get("prediction_type", "")
        repeat_signal = topology in {"DTR", "ITR"} or repeat_prediction in {"DTR", "ITR"}
        at_edge = left == 0 or right == 0
        cutpoint_case = bool(partial and partial != "00" and repeat_signal and at_edge)

        conflicts = []
        if protein not in faa:
            conflicts.append("protein_id_missing_from_fasta")
        if protein in malformed:
            conflicts.append("malformed_fasta_header")
        if faa_occ[protein] > 1:
            conflicts.append("duplicate_protein_id_in_fasta")
        if protein not in index:
            conflicts.append("protein_id_missing_from_gene_index")
        if index_occ[protein] > 1:
            conflicts.append("duplicate_protein_id_in_gene_index")
        if contig not in contig_lengths:
            conflicts.append("contig_id_missing_from_source_fasta")
        if contig_occ[contig] > 1:
            conflicts.append("duplicate_contig_id_in_source_fasta")
        if contig not in quality:
            conflicts.append("contig_id_missing_from_quality_table")
        if start is not None and end is not None and length is not None and not (1 <= start <= end <= length):
            conflicts.append("coordinates_out_of_bounds")
        if idx and protein_rec:
            if idx.get("contig") != contig:
                conflicts.append("candidate_gene_index_contig_conflict")
            if numeric(idx.get("contig_pos_start"), int) != start or numeric(idx.get("contig_pos_end"), int) != end:
                conflicts.append("fasta_gene_index_coordinate_conflict")
            idx_strand = "+" if idx.get("frame") == "1" else "-" if idx.get("frame") == "-1" else idx.get("frame", "")
            if idx_strand != strand:
                conflicts.append("fasta_gene_index_strand_conflict")
        if protein_rec and partial not in {"00", "01", "10", "11"}:
            conflicts.append("partial_attribute_missing_or_invalid")
        if numeric(protein_rec.get("internal_stop_count"), int) not in (None, 0):
            conflicts.append("internal_stop_in_protein")
        if protein_rec.get("nonstandard_aa"):
            conflicts.append("nonstandard_amino_acid")
        if protein_rec.get("translation_length_consistent") == "false":
            conflicts.append("translation_coordinate_length_conflict")

        if conflicts:
            decision = "review"
            reason = ";".join(conflicts)
        elif cutpoint_case:
            decision = "review"
            reason = "partial_at_terminal_repeat_boundary;DTR_or_ITR_not_proof_of_circularity"
        elif partial == "00":
            decision = "predicted_complete"
            reason = "partial=00;exact_source_and_coordinate_match"
        else:
            decision = "predicted_partial"
            reason = f"partial={partial}"

        theme, theme_method = assign_theme(source, themes[source["AVG_category"]])
        completeness = numeric(q.get("completeness"))
        qlen = numeric(q.get("contig_length"), int)
        virus_score = numeric(q.get("virus_score"))
        hallmarks = numeric(q.get("n_hallmarks"), int)
        cand_score = numeric(source.get("virus_score"))
        cand_hallmarks = numeric(source.get("n_hallmarks"), int)
        row = dict(source)
        row.update(protein_rec)
        row.update({
            "gene_index_contig": idx.get("contig", ""),
            "gene_index_start": idx.get("contig_pos_start", ""),
            "gene_index_end": idx.get("contig_pos_end", ""),
            "gene_index_frame": idx.get("frame", ""),
            "contig_length": length if length is not None else "",
            "distance_left_bp": left if left is not None else "",
            "distance_right_bp": right if right is not None else "",
            "within_100bp_of_edge": bool_text(near_edge),
            "checkv_topology": topology,
            "checkamg_repeat_prediction": repeat_prediction,
            "repeat_length": circ.get("repeat_length", ""),
            "terminal_repeat_signal": bool_text(repeat_signal),
            "potential_artificial_cutpoint_case": bool_text(cutpoint_case),
            "checkv_quality": q.get("checkv_quality", ""),
            "miuvig_quality": q.get("miuvig_quality", ""),
            "checkv_completeness": q.get("completeness", ""),
            "checkv_completeness_method": q.get("completeness_method", ""),
            "checkv_contamination": q.get("contamination", ""),
            "quality_table_contig_length": q.get("contig_length", ""),
            "quality_table_virus_score": q.get("virus_score", ""),
            "quality_table_n_hallmarks": q.get("n_hallmarks", ""),
            "passes_checkv_completeness_ge90": bool_text(completeness is not None and completeness >= 90),
            "passes_length_gt10kb": bool_text(length is not None and length > 10000),
            "passes_virus_score_gt0.90": bool_text(virus_score is not None and virus_score > 0.90),
            "passes_n_hallmarks_ge1": bool_text(hallmarks is not None and hallmarks >= 1),
            "candidate_quality_values_match": bool_text(
                cand_score is not None and virus_score is not None and abs(cand_score - virus_score) < 1e-12
                and cand_hallmarks == hallmarks
            ),
            "source_contig_length_matches_quality_table": bool_text(length is not None and qlen is not None and length == qlen),
            "functional_theme_L1": theme.get("category_L1", "") or "Unmapped",
            "functional_theme_L2": theme.get("category_L2", ""),
            "functional_theme_L3": theme.get("category_L3", ""),
            "functional_theme_mapping_method": theme_method,
            "orf_qc_decision": decision,
            "orf_qc_reason": reason,
            "viral_origin_conclusion": source.get("evidence_route", ""),
            "functional_identity_conclusion": source.get("review_status", "") or "annotation_not_manually_resolved",
        })
        audit_rows.append(row)

    audit_extra = [
        "fasta_header", "orf_start", "orf_end", "strand", "strand_numeric", "partial",
        "partial_5prime", "partial_3prime", "start_type", "rbs_motif", "rbs_spacer",
        "gc_cont", "genetic_code_header", "protein_length_aa", "terminal_stop_present",
        "internal_stop_count", "nonstandard_aa", "translation_length_consistent",
        "gene_index_contig", "gene_index_start", "gene_index_end", "gene_index_frame",
        "contig_length", "distance_left_bp", "distance_right_bp", "within_100bp_of_edge",
        "checkv_topology", "checkamg_repeat_prediction", "repeat_length", "terminal_repeat_signal",
        "potential_artificial_cutpoint_case", "checkv_quality", "miuvig_quality",
        "checkv_completeness", "checkv_completeness_method", "checkv_contamination",
        "quality_table_contig_length", "quality_table_virus_score", "quality_table_n_hallmarks",
        "passes_checkv_completeness_ge90", "passes_length_gt10kb", "passes_virus_score_gt0.90",
        "passes_n_hallmarks_ge1", "candidate_quality_values_match",
        "source_contig_length_matches_quality_table", "functional_theme_L1", "functional_theme_L2",
        "functional_theme_L3", "functional_theme_mapping_method", "viral_origin_conclusion",
        "orf_qc_decision", "orf_qc_reason", "functional_identity_conclusion",
    ]
    fields = list(candidates[0]) + audit_extra
    write_tsv(args.outdir / "AVG_ORF_audit.tsv", audit_rows, fields)
    complete = [r for r in audit_rows if r["orf_qc_decision"] == "predicted_complete"]
    partial = [r for r in audit_rows if r["orf_qc_decision"] == "predicted_partial"]
    review = [r for r in audit_rows if r["orf_qc_decision"] == "review"]
    write_tsv(args.outdir / "AVG_predicted_complete.tsv", complete, fields)
    write_tsv(args.outdir / "AVG_predicted_partial.tsv", partial, fields)
    write_tsv(args.outdir / "AVG_ORF_review.tsv", review, fields)

    summary_rows = []
    for category in ("AMG", "APG", "AReG", "TOTAL"):
        subset = audit_rows if category == "TOTAL" else [r for r in audit_rows if r["AVG_category"] == category]
        for decision in ("all_candidates", "predicted_complete", "predicted_partial", "review"):
            rows = subset if decision == "all_candidates" else [r for r in subset if r["orf_qc_decision"] == decision]
            summary_rows.append({"AVG_category": category, "ORF_status": decision, "protein_records": len(rows), "unique_vOTUs": len({r["Contig"] for r in rows})})
    write_tsv(args.outdir / "AVG_ORF_summary.tsv", summary_rows, ["AVG_category", "ORF_status", "protein_records", "unique_vOTUs"])

    reason_rows = []
    for category in ("AMG", "APG", "AReG", "TOTAL"):
        subset = audit_rows if category == "TOTAL" else [r for r in audit_rows if r["AVG_category"] == category]
        counts = Counter(r["orf_qc_reason"] for r in subset if r["orf_qc_decision"] != "predicted_complete")
        for reason, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
            affected = [r for r in subset if r["orf_qc_reason"] == reason]
            reason_rows.append({"AVG_category": category, "reason": reason, "protein_records": count, "unique_vOTUs": len({r["Contig"] for r in affected})})
    write_tsv(args.outdir / "AVG_ORF_reason_summary.tsv", reason_rows, ["AVG_category", "reason", "protein_records", "unique_vOTUs"])

    checkamg_start = count_checkamg_start([
        ("AMG", args.checkamg_metabolic), ("APG", args.checkamg_physiology), ("AReG", args.checkamg_regulation)
    ])
    stage_rows = []
    for category in ("AMG", "APG", "AReG", "TOTAL"):
        strict_rows = audit_rows if category == "TOTAL" else [r for r in audit_rows if r["AVG_category"] == category]
        complete_rows = [r for r in strict_rows if r["orf_qc_decision"] == "predicted_complete"]
        for stage, proteins_n, contigs_n in (
            ("CheckAMG curated candidates", *checkamg_start[category]),
            ("Strict viral candidates", len(strict_rows), len({r["Contig"] for r in strict_rows})),
            ("Predicted complete ORFs", len(complete_rows), len({r["Contig"] for r in complete_rows})),
        ):
            stage_rows.append({"stage": stage, "AVG_category": category, "protein_records": proteins_n, "unique_vOTUs": contigs_n})
    write_tsv(args.outdir / "AVG_stage_counts.tsv", stage_rows, ["stage", "AVG_category", "protein_records", "unique_vOTUs"])

    quality_rows = []
    for contig in sorted(wanted_contigs):
        q = quality.get(contig, {})
        length = contig_lengths.get(contig)
        completeness = numeric(q.get("completeness"))
        score = numeric(q.get("virus_score"))
        hallmark = numeric(q.get("n_hallmarks"), int)
        quality_rows.append({
            "Contig": contig, "source_contig_length": length if length is not None else "",
            "checkv_completeness": q.get("completeness", ""), "checkv_completeness_method": q.get("completeness_method", ""),
            "checkv_quality": q.get("checkv_quality", ""), "miuvig_quality": q.get("miuvig_quality", ""),
            "checkv_topology": q.get("topology", ""), "geNomad_virus_score": q.get("virus_score", ""),
            "geNomad_n_hallmarks": q.get("n_hallmarks", ""),
            "passes_checkv_completeness_ge90": bool_text(completeness is not None and completeness >= 90),
            "passes_length_gt10kb": bool_text(length is not None and length > 10000),
            "passes_virus_score_gt0.90": bool_text(score is not None and score > 0.90),
            "passes_n_hallmarks_ge1": bool_text(hallmark is not None and hallmark >= 1),
            "passes_all_four_quality_checks": bool_text(
                completeness is not None and completeness >= 90 and length is not None and length > 10000
                and score is not None and score > 0.90 and hallmark is not None and hallmark >= 1
            ),
        })
    q_fields = list(quality_rows[0])
    write_tsv(args.outdir / "AVG_viral_quality_audit.tsv", quality_rows, q_fields)
    threshold_rows = []
    for field in q_fields[9:]:
        count = sum(r[field] == "true" for r in quality_rows)
        threshold_rows.append({"criterion": field, "passing_vOTUs": count, "total_vOTUs": len(quality_rows), "percent": f"{100*count/len(quality_rows):.3f}"})
    write_tsv(args.outdir / "AVG_viral_quality_summary.tsv", threshold_rows, ["criterion", "passing_vOTUs", "total_vOTUs", "percent"])

    theme_stats = []
    ko_stats = []
    for scope, rows in (("all_strict_candidates", audit_rows), ("predicted_complete_ORFs", complete)):
        grouped = defaultdict(list)
        for row in rows:
            key = (row["AVG_category"], row["functional_theme_L1"], row["functional_theme_L2"], row["functional_theme_L3"])
            grouped[key].append(row)
        for key, members in grouped.items():
            theme_stats.append({"scope": scope, "AVG_category": key[0], "theme_L1": key[1], "theme_L2": key[2], "theme_L3": key[3], "protein_records": len({m["Protein"] for m in members}), "unique_vOTUs": len({m["Contig"] for m in members})})
        ko_group = defaultdict(list)
        for row in rows:
            for ko in sorted(set(KO_RE.findall(row.get("KEGG_hmm_id", "")))):
                ko_group[(row["AVG_category"], ko, row.get("KEGG_Description", ""), row["functional_theme_L1"], row["functional_theme_L2"])].append(row)
        for key, members in ko_group.items():
            ko_stats.append({"scope": scope, "AVG_category": key[0], "KO": key[1], "KO_description": key[2], "theme_L1": key[3], "theme_L2": key[4], "protein_records": len({m["Protein"] for m in members}), "unique_vOTUs": len({m["Contig"] for m in members})})
    write_tsv(args.outdir / "AVG_functional_theme_stats.tsv", theme_stats, ["scope", "AVG_category", "theme_L1", "theme_L2", "theme_L3", "protein_records", "unique_vOTUs"])
    write_tsv(args.outdir / "AVG_KO_pathway_stats.tsv", ko_stats, ["scope", "AVG_category", "KO", "KO_description", "theme_L1", "theme_L2", "protein_records", "unique_vOTUs"])

    match_rows = [
        {"check": "candidate_rows", "count": len(candidates)},
        {"check": "unique_candidate_protein_ids", "count": len(wanted_proteins)},
        {"check": "duplicate_candidate_rows_by_protein", "count": sum(v - 1 for v in candidate_occ.values() if v > 1)},
        {"check": "protein_ids_matched_to_fasta", "count": len(wanted_proteins & set(faa))},
        {"check": "protein_ids_missing_from_fasta", "count": len(wanted_proteins - set(faa))},
        {"check": "candidate_protein_ids_duplicated_in_fasta", "count": sum(v > 1 for v in faa_occ.values())},
        {"check": "protein_ids_matched_to_gene_index", "count": len(wanted_proteins & set(index))},
        {"check": "protein_ids_missing_from_gene_index", "count": len(wanted_proteins - set(index))},
        {"check": "candidate_protein_ids_duplicated_in_gene_index", "count": sum(v > 1 for v in index_occ.values())},
        {"check": "contig_ids_matched_to_source_fasta", "count": len(wanted_contigs & set(contig_lengths))},
        {"check": "contig_ids_missing_from_source_fasta", "count": len(wanted_contigs - set(contig_lengths))},
        {"check": "candidate_contig_ids_duplicated_in_source_fasta", "count": sum(v > 1 for v in contig_occ.values())},
        {"check": "contig_ids_matched_to_quality_table", "count": len(wanted_contigs & set(quality))},
        {"check": "contig_ids_missing_from_quality_table", "count": len(wanted_contigs - set(quality))},
    ]
    write_tsv(args.outdir / "AVG_source_matching_summary.tsv", match_rows, ["check", "count"])

    missing_rows = []
    for value in sorted(wanted_proteins - set(faa)):
        missing_rows.append({"record_type": "protein_missing_fasta", "identifier": value})
    for value in sorted(wanted_proteins - set(index)):
        missing_rows.append({"record_type": "protein_missing_gene_index", "identifier": value})
    for value in sorted(wanted_contigs - set(contig_lengths)):
        missing_rows.append({"record_type": "contig_missing_source_fasta", "identifier": value})
    for value in sorted(wanted_contigs - set(quality)):
        missing_rows.append({"record_type": "contig_missing_quality_table", "identifier": value})
    write_tsv(args.outdir / "AVG_missing_and_conflicting_IDs.tsv", missing_rows, ["record_type", "identifier"])

    manifest_rows = []
    for role, path in (
        ("strict_candidate_table", args.candidates), ("pyrodigal_gv_filtered_protein_fasta", args.protein_fasta),
        ("checkamg_source_contig_fasta", args.contig_fasta), ("checkamg_gene_index", args.gene_index),
        ("checkamg_terminal_repeat_table", args.circular_table), ("checkv_genomad_host_quality_table", args.quality_table),
    ):
        stat = path.stat()
        manifest_rows.append({"role": role, "path": str(path), "size_bytes": stat.st_size, "mtime": datetime.fromtimestamp(stat.st_mtime).isoformat()})
    write_tsv(args.outdir / "source_manifest.tsv", manifest_rows, ["role", "path", "size_bytes", "mtime"])

    total_by_status = Counter(r["orf_qc_decision"] for r in audit_rows)
    quality_lookup = {r["criterion"]: r for r in threshold_rows}
    report = f"""# AVG ORF completeness audit

## Scope

This audit retains the full strict candidate universe: **{len(candidates):,} protein records on {len(wanted_contigs):,} vOTUs**. No candidate was removed because its annotation was broad, domain-only, `-like`, or functionally disputed.

Input category counts are AMG **{sum(r['AVG_category']=='AMG' for r in candidates):,}**, APG **{sum(r['AVG_category']=='APG' for r in candidates):,}**, and AReG **{sum(r['AVG_category']=='AReG' for r in candidates):,}**.

## Source correspondence

- Exact full Protein-ID matches to retained pyrodigal-gv FASTA: **{len(wanted_proteins & set(faa)):,}/{len(wanted_proteins):,}**.
- Exact full Protein-ID matches to CheckAMG gene index: **{len(wanted_proteins & set(index)):,}/{len(wanted_proteins):,}**.
- Exact contig-ID matches to the CheckAMG source FASTA: **{len(wanted_contigs & set(contig_lengths)):,}/{len(wanted_contigs):,}**.
- Candidate protein IDs duplicated in the FASTA: **{sum(v > 1 for v in faa_occ.values()):,}**; in the gene index: **{sum(v > 1 for v in index_occ.values()):,}**.
- CheckAMG metadata records that this protein FASTA derives from the same `run_pyrodigal_gv` step and subsequent CDS filtering. The original temporary pyrodigal-gv directory, GFF, and CDS nucleotide FASTA are no longer present. Coordinates, strand, start type, and `partial` attributes remain in the retained protein FASTA headers.
- Predictor versions in the preserved environment: pyrodigal-gv 0.3.2 and pyrodigal 3.5.2.

## ORF result

- Predicted complete (`partial=00` with matching sources and valid coordinates): **{total_by_status['predicted_complete']:,}**.
- Explicitly partial: **{total_by_status['predicted_partial']:,}**.
- Review required (missing/conflicting evidence, sequence anomaly, or terminal-repeat boundary case): **{total_by_status['review']:,}**.

The two `partial` digits are interpreted as 5-prime and 3-prime incompleteness flags. Edge proximity (<=100 bp) is reported descriptively and was not used as an exclusion rule. DTR/ITR is treated only as terminal-repeat evidence, not proof of circularity; an explicitly partial ORF at such a boundary is placed in review rather than automatically restored.

## Viral-quality cross-check (not applied to the ORF subset)

Among {len(quality_rows):,} vOTUs:

- CheckV completeness >=90%: **{quality_lookup['passes_checkv_completeness_ge90']['passing_vOTUs']:,}**.
- Length >10 kb: **{quality_lookup['passes_length_gt10kb']['passing_vOTUs']:,}**.
- geNomad virus score >0.90: **{quality_lookup['passes_virus_score_gt0.90']['passing_vOTUs']:,}**.
- geNomad hallmark count >=1: **{quality_lookup['passes_n_hallmarks_ge1']['passing_vOTUs']:,}**.
- All four simultaneously: **{quality_lookup['passes_all_four_quality_checks']['passing_vOTUs']:,}**.

These quality checks remain separate from ORF completeness. A complete ORF is not by itself evidence of complete function or confirmed AVG identity.

## Functional statistics

Functional themes are mapped using CheckAMG's bundled, version-traceable `AMGs_categorized.tsv`, `APGs_categorized.tsv`, and `AReGs_categorized.tsv`. The original function, `-like`, family, and domain labels remain unchanged in all record-level outputs. KO statistics count unique Protein IDs per KO and report unique vOTUs separately; they are not sample-read abundance.
"""
    (args.outdir / "REPORT.md").write_text(report, encoding="utf-8")
    print(f"candidates={len(candidates)} complete={len(complete)} partial={len(partial)} review={len(review)}")


if __name__ == "__main__":
    main()
