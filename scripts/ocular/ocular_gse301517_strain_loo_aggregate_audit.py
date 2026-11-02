#!/usr/bin/env python3
"""Recheck per-strain direction and leave-one-strain-out stability in GSE301517.

All source inputs come from an already-verified public *aggregate-only* GitHub
Actions artifact (#13037335888, run #41186736618). Source study offers 3 TM
strain groups with 4/2/1 technical samples/arm; NOT 7 donors. No raw sequencing
files, hidden clinical data or human genotype tables are used.
"""
import argparse
import csv
import hashlib
import json
import math
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

STRAINS = ("TM62", "TM70", "TM72")
CORE = ("MYOC", "TGFBI", "MMP2", "ACTA2")
PRIOR_CANDIDATES = ("MYOC", "TGFBI", "MMP2")
PIN = {
    "GSE301517_core4_3study_direction.tsv":
        "926b7899639d9b759d03688b911b89becf6b776103927a226e4f301415a86180",
    "GSE301517_EMT_194_3study_signed_direction.tsv":
        "edac2a23c47e40fde419214c8c45fed6f1bd33861a72f8c37c789db2f7d88e1d",
}
REQUIRED = ("GSE301517_16_p9_panel_3strain_vs_48h.tsv",
            "GSE301517_EMT_shared_3strain_direction.tsv",
            "GSE301517_core4_3study_direction.tsv",
            "GSE301517_EMT_194_3study_signed_direction.tsv",
            "GSE301517_3strain_vs_48h_evidence_gate.json")


def signed_match(x, y):
    return x * y > 0


def numeric(x):
    a = float(x)
    if not math.isfinite(a):
        raise ValueError("NONFINITE_EFFECT")
    return a


def read_tsv(path):
    with path.open("r", encoding="utf8", newline="") as fh:
        out = list(csv.DictReader(fh, delimiter="\t"))
    if not out:
        raise ValueError("EMPTY_GENE_TABLE_" + path.name)
    return out


def write_tsv(path, records):
    if not records:
        raise ValueError("EMPTY_OUTPUT")
    fields = sorted(set().union(*(set(x) for x in records)))
    with path.open("w", encoding="utf8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(records)


def compute_panel(panel, core):
    found = {r["gene_symbol"]: r for r in panel}
    if len(panel) != 16 or len(found) != 16 or not set(CORE).issubset(found):
        raise ValueError("EXPECTED_16_UNIQUE_GENES_AND_CORE4")
    c = {r["gene_symbol"]: r for r in core}
    if len(core) != 4 or set(c) != set(CORE):
        raise ValueError("EXPECTED_CORE4_VALUES")
    output = []
    for g, row in sorted(found.items()):
        acute = numeric(row["primary_500nM_48h_log2FC"])
        seven_day = numeric(row["third_500nM_7d_log2_CPM_ratio"])
        per = [numeric(row[f"third_{strain}_log2_CPM_ratio"]) for strain in STRAINS]
        reported = int(row["three_strain_sign_match_n"])
        actual = sum(signed_match(acute, x) for x in per)
        if reported != actual:
            raise ValueError("STRAIN_SIGN_QC_DISAGREES_" + g)
        if abs(sum(per) / 3 - seven_day) > .00003:
            raise ValueError("7DAY_EQUAL_STRAIN_MEAN_DRIFT_" + g)
        leave_one = [(sum(per) - v) / 2 for v in per]
        mean_positive = signed_match(acute, seven_day)
        all_loo = all(signed_match(acute, v) for v in leave_one)
        # Sign might be stable but affected by one unusually large strain.
        result = {
            "gene_symbol": g, "gene_id": row["gene_id"],
            "prior_candidate_core_three": g in PRIOR_CANDIDATES,
            "acute_500nM_48h_edgeR_log2FC": acute,
            "seven_day_3strain_CPM_pseudocount_log2ratio": seven_day,
            "three_of_three_strains_same_direction_as_acute": actual == 3,
            "n_strains_same_direction_as_acute": actual,
            "seven_day_mean_matches_acute": mean_positive,
            "seven_day_all_three_leave_one_strain_out_same_direction_as_acute": all_loo,
            "two_of_three_7day_strains_consensus": actual >= 2,
            "7day_log2ratio_estimand_not_edgeR_log2FC": True,
            "independent_7day_p_or_fdr_available": False,
        }
        for idx, strain in enumerate(STRAINS):
            result[f"{strain}_7day_log2ratio"] = per[idx]
            result[f"exclude_{strain}_other_two_equal_mean_log2ratio"] = leave_one[idx]
            result[f"exclude_{strain}_same_direction_as_acute"] = signed_match(acute, leave_one[idx])
        if g in c:
            ext = numeric(c[g]["published_100nM_48h_log2FC"])
            same = signed_match(acute, ext) and signed_match(ext, seven_day)
            if str(c[g]["all_three_same_direction"]).lower() != str(same).lower():
                raise ValueError("CORE_THREE_STUDY_SIGN_MISMATCH_" + g)
            result["mehrotra_100nM_48h_published_log2FC"] = ext
            result["all_three_study_mean_direction_same"] = same
            result["all_three_7day_strains_and_two_48h_studies_same"] = (
                signed_match(acute, ext) and all(signed_match(ext, v) for v in per))
        else:
            result["mehrotra_100nM_48h_published_log2FC"] = None
            result["all_three_study_mean_direction_same"] = None
            result["all_three_7day_strains_and_two_48h_studies_same"] = None
        output.append(result)
    return output


def compute_emt(emt, three):
    m = {row["gene_id"].split(".")[0]: row for row in emt}
    t = {row["gene_id"].split(".")[0]: row for row in three}
    if len(m) != 194 or len(t) != 194 or set(m) != set(t):
        raise ValueError("EXPECTED_194_SAME_ENSEMBL_EMT_IDS")
    rows = []
    for gene in sorted(m):
        row, cross = m[gene], t[gene]
        acute = numeric(cross["acute_500nM_48h_log2FC"])
        ext = numeric(cross["external_100nM_48h_log2FC"])
        seven = numeric(cross["third_500nM_7d_CPM_log2ratio"])
        if abs(acute - numeric(row["primary_48h_log2FC"])) > .000002:
            raise ValueError("EMT_ACUTE_SOURCE_MISMATCH")
        if abs(seven - numeric(row["third_7d_CPM_log2ratio"])) > .000002:
            raise ValueError("EMT_7DAY_SOURCE_MISMATCH")
        count = int(row["strains_3_matching_primary"])
        if not (0 <= count <= 3):
            raise ValueError("OUT_OF_RANGE_STRAIN_COUNT")
        all3 = signed_match(acute, ext) and signed_match(ext, seven)
        if str(cross["three_studies_all_same_direction"]).lower() != str(all3).lower():
            raise ValueError("EMT_STUDY_SIGN_MISMATCH")
        rows.append({
            "gene_id": gene, "gene_symbol": row["gene_symbol"],
            "three_studies_same_signed": all3,
            "acute_mehrotra_48h_same_signed": signed_match(acute, ext),
            "7day_equal_strain_mean_vs_acute_same_signed": signed_match(acute, seven),
            "7day_strains_same_as_acute_n": count,
            "7day_three_strain_unanimous": count == 3,
            "three_study_same_and_7day_3of3_strains":
                all3 and count == 3,
            "single_study_7day_values_not_in_this_emt_file":
                True,
        })
    return rows


def selftest():
    per = [4., 3., -1.]
    acute = 1.
    assert sum(signed_match(acute, x) for x in per) == 2
    assert all(signed_match(acute, (sum(per) - x) / 2) for x in per)
    assert not all(signed_match(acute, x) for x in per)
    assert not signed_match(0, 1) and not signed_match(1, 0)
    # One extreme strong strain can flip a three-strain mean, exposed by LOO.
    per = [-1., -1., 5.]
    assert signed_match(1., sum(per) / 3)
    assert not all(signed_match(1., (sum(per) - x) / 2) for x in per)
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "test.tsv"
        p.write_text("gene_id\tgene_symbol\nENSG00000000001\tMYOC\n", encoding="utf8")
        assert len(read_tsv(p)) == 1
    print("OCULAR_STRAIN_LOO_SELFTEST_PASS", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path)
    p.add_argument("--out", type=Path)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        selftest()
        return
    if a.source is None or a.out is None:
        p.error("--source and --out required")
    src = {}
    for name in REQUIRED:
        files = list(a.source.rglob(name))
        if len(files) != 1:
            raise ValueError("EXPECTED_SINGLE_ARCHIVED_FILE " + name)
        src[name] = files[0]
    hashes = {name: hashlib.sha256(f.read_bytes()).hexdigest() for name, f in src.items()}
    for name, sha in PIN.items():
        if hashes[name] != sha:
            raise ValueError("PINNED_2026_THREE_STUDY_SOURCE_HASH_DRIFT " + name)
    gate = json.loads(src["GSE301517_3strain_vs_48h_evidence_gate.json"].read_text(encoding="utf8"))
    if gate["status"] != "PASS_SOURCE_INTEGRITY_AND_DESCRIPTIVE_CROSS_STUDY_CPM":
        raise ValueError("UNVERIFIED_ORIGINAL_THREE_STUDY_SOURCE")
    if gate["GSE301517_strains"] != list(STRAINS) or gate["GSE301517_technical_replicates_per_strain"] != {"TM62":4,"TM70":2,"TM72":1}:
        raise ValueError("THREE_STRAIN_UNIT_OR_REPLICATE_DRIFT")
    panel = compute_panel(read_tsv(src["GSE301517_16_p9_panel_3strain_vs_48h.tsv"]),
                          read_tsv(src["GSE301517_core4_3study_direction.tsv"]))
    emt = compute_emt(read_tsv(src["GSE301517_EMT_shared_3strain_direction.tsv"]),
                      read_tsv(src["GSE301517_EMT_194_3study_signed_direction.tsv"]))
    if sum(x["three_studies_same_signed"] for x in emt) != 124:
        raise ValueError("FROZEN_THREE_STUDY_EMT124_DRIFT")
    result = {
        "status": "PASS_DESCRIPTIVE_7DAY_THREE_STRAIN_LEAVE_ONE_OUT_AUDIT",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_run_id": 41186736618, "source_artifact_id": 13037335888,
        "source_sha256": hashes,
        "source_study": "GSE301517",
        "seven_day_independent_unit": "three named TM strains (not seven technical replicates; donor independence unverified)",
        "technical_libraries_per_condition": {"TM62":4,"TM70":2,"TM72":1},
        "different_studies_have_different_drug_doses_or_exposure_times": True,
        "seven_day_log2ratio_not_same_estimand_as_edgeR_log2FC": True,
        "not_a_per_strain_statistical_significance_test": True,
        "core_four": [p for p in panel if p["gene_symbol"] in CORE],
        "panel_16": {
            "genes":len(panel),
            "three_strains_unanimous_n":sum(p["three_of_three_strains_same_direction_as_acute"] for p in panel),
            "7day_mean_matches_acute_n":sum(p["seven_day_mean_matches_acute"] for p in panel),
            "all_leave_one_strain_out_mean_match_acute_n":
                sum(p["seven_day_all_three_leave_one_strain_out_same_direction_as_acute"] for p in panel),
        },
        "emt_194": {
            "genes":len(emt),
            "all_three_study_signed":sum(x["three_studies_same_signed"] for x in emt),
            "7day_3of3_strains_same_as_acute":sum(x["7day_three_strain_unanimous"] for x in emt),
            "three_study_and_7day_3of3":sum(x["three_study_same_and_7day_3of3_strains"] for x in emt),
            "distribution_n_7day_strains_same_as_acute":
                dict(sorted(Counter(str(x["7day_strains_same_as_acute_n"]) for x in emt).items())),
        },
        "not_original_author_80_gene_comparison": True,
        "no_author_outreach_needed": True,
        "not_causal_dose_or_time_effect": True,
        "not_clinical_steroid_iop": True,
    }
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "GSE301517_16gene_three_strain_LOO_source_audit.tsv", panel)
    write_tsv(a.out / "GSE301517_EMT194_per_strain_agreement_categories.tsv", emt)
    (a.out / "GSE301517_strain_LOO_audit_manifest.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    lines = ["# Cross-study TM 7-day strain contribution audit",
             "Data source: previously verified public gene-level aggregate only.",
             "This is descriptive context-sensitivity, not proof of donor independence, causal treatment effect or a dose/time interaction.",
             "",
             "## Candidate core under per-strain scrutiny",
             "| Gene | Acute 500nM 48h | Mehrotra 100nM 48h | GSE301517 7d 3-strain mean | Strains same as acute | 3/3 LOO means same as acute |",
             "|---|---:|---:|---:|---:|---|"]
    for x in result["core_four"]:
        lines.append("| {gene_symbol} | {acute_500nM_48h_edgeR_log2FC:+.3f} | {mehrotra_100nM_48h_published_log2FC:+.3f} | {seven_day_3strain_CPM_pseudocount_log2ratio:+.3f} | {n_strains_same_direction_as_acute}/3 | {seven_day_all_three_leave_one_strain_out_same_direction_as_acute} |".format(**x))
    lines.extend(["", "## EMT gene set", json.dumps(result["emt_194"], indent=2),
                  "", "No author DEG80 list was available or needed. Published 7-day CPM ratios are not gene-wise inferential log2FC."])
    (a.out / "GSE301517_strain_LOO_report.md").write_text("\n".join(lines) + "\n", encoding="utf8")
    print("OCULAR_REAL_STRAIN_LOO_PASS", json.dumps({
        "core_four":result["core_four"], "panel":result["panel_16"],
        "emt":result["emt_194"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
