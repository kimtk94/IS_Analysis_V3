#!/usr/bin/env python3
"""Read-only source audit for the ischemic-stroke candidate and molecular evidence universe.

Writes new audit artifacts only inside --out. It does not rerun coloc,
change canonical inputs, or equate provisional distance clusters with loci.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT_FILES = {
    "abf": "results/is/stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
    "expanded_regions": "results/is/stage5_functional/broad_discovery_v2_ancestry/IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv",
    "expanded_genes": "results/is/stage5_functional/broad_discovery_v2_ancestry/IS_ANCESTRY_GENE_EVIDENCE_V3.tsv",
    "expanded_genes_summary": "results/is/stage5_functional/broad_discovery_v2_ancestry/IS_ANCESTRY_GENE_EVIDENCE_V3_SUMMARY.json",
    "eas_summary": "results/is/stage5_functional/broad_discovery_v1/IS_ALL_GENE_UNIVERSE_SUMMARY.json",
    "expanded_summary": "results/is/stage5_functional/broad_discovery_v2_ancestry/IS_ANCESTRY_EXPANSION_DISCOVERY_SUMMARY.json",
    "human_presence": "results/is/stage5_functional/phase9f_e_human_r3_4_import/HUMAN_TARGET_GENE_PRESENCE.tsv",
    "human_assay": "results/is/stage5_functional/phase9f_e_human_r3_4_import/HUMAN_ASSAY_AUDIT.tsv",
    "human_metadata": "results/is/stage5_functional/phase9f_e_human_r3_4_import/HUMAN_RUN_MANIFEST.json",
    "human_celltypes": "results/is/stage5_functional/phase9f_e_human_r3_4_import/HUMAN_AUTHOR_ANNOTATION_COUNTS.tsv",
    "phase11": "results/is/audits/IS_PHASE11_COLOC_EVIDENCE_AUDIT_20261009.json",
}
CORE = {"FGF5", "ALDH2", "SH3PXD2A", "COL4A2"}
SECONDARY = {"CALHM2", "NEURL1", "C4orf22", "INA", "COL4A1"}

def rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")
        return list(reader)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def number(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None

def analyze_abf(data):
    by_locus = Counter()
    by_status = Counter()
    by_gene = defaultdict(list)
    datasets = set()
    keys = set()
    n_05 = n_075 = n_08 = 0
    pass_count = 0
    duplicates = 0
    for row in data:
        by_locus[row.get("locus", "")] += 1
        status = row.get("status", "")
        by_status[status] += 1
        key = (row.get("locus"), row.get("dataset_key"), row.get("gene_base"))
        if key in keys:
            duplicates += 1
        keys.add(key)
        datasets.add(row.get("dataset_key", ""))
        if status != "PASS":
            continue
        pass_count += 1
        h4 = number(row.get("PP.H4"))
        if h4 is None or not 0 <= h4 <= 1:
            raise ValueError(f"Invalid H4 for {key}")
        n_05 += h4 >= 0.5
        n_075 += h4 >= 0.75
        n_08 += h4 >= 0.8
        gene = row.get("gene_symbol") or row.get("gene_base") or ""
        by_gene[gene].append(row)
    best = {}
    for gene, gene_rows in by_gene.items():
        winner = max(gene_rows, key=lambda r: float(r["PP.H4"]))
        h3 = number(winner.get("PP.H3"))
        h4 = number(winner.get("PP.H4"))
        best[gene] = {
            "locus": winner["locus"], "best_tissue": winner["dataset_key"],
            "h4": h4, "h3": h3,
            "h4_h3_ratio": round(h4 / h3, 5) if h3 and h3 > 0 else None,
            "pass_test_count": len(gene_rows),
        }
    return {
        "total_test_rows": len(data), "pass_test_rows": pass_count,
        "status_counts": dict(sorted(by_status.items())),
        "counts_by_locus": dict(sorted(by_locus.items())),
        "unique_test_keys": len(keys), "duplicate_test_keys": duplicates,
        "unique_genes_tested": len(by_gene), "unique_tissue_datasets": len(datasets),
        "h4_ge_0_5": n_05, "h4_ge_0_75": n_075, "h4_ge_0_8": n_08,
        "best_by_gene": dict(sorted(best.items())),
        "multiple_testing_adjustment": "NOT_PERFORMED",
        "remark": "Threshold counts are descriptive, not family-wise/FDR-adjusted discoveries",
    }

def partition_candidate_rows(data):
    mapped, anchor_only, unresolved = [], [], []
    for row in data:
        if row.get("group_id") and row.get("gene_id_stable"):
            mapped.append(row)
        elif (row.get("source") == "LEGACY_ANCHOR_ONLY"
              and not row.get("group_id") and not row.get("gene_id_stable")):
            anchor_only.append(row)
        else:
            unresolved.append(row)
    if unresolved:
        raise ValueError(f"Unexpected unmapped candidate rows: {len(unresolved)}")
    return mapped, anchor_only

def audit(root):
    files = {k: root / v for k, v in ROOT_FILES.items()}
    missing = [k for k, p in files.items() if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing required source files: " + ", ".join(missing))
    digests = {k: {"path": str(v), "sha256": sha256(v)} for k, v in files.items()}
    abf = analyze_abf(rows(files["abf"]))
    regions = rows(files["expanded_regions"])
    candidate_rows = rows(files["expanded_genes"])
    genes, anchor_only = partition_candidate_rows(candidate_rows)
    genes_summary = json.loads(files["expanded_genes_summary"].read_text(encoding="utf-8"))
    eas = json.loads(files["eas_summary"].read_text(encoding="utf-8"))
    expansion = json.loads(files["expanded_summary"].read_text(encoding="utf-8"))
    cells = rows(files["human_celltypes"])
    presence = rows(files["human_presence"])
    human_assay = rows(files["human_assay"])
    human = json.loads(files["human_metadata"].read_text(encoding="utf-8"))
    phase11 = json.loads(files["phase11"].read_text(encoding="utf-8"))
    groups = {r["group_id"] for r in regions}
    pairs = {(r["group_id"], r["gene_id_stable"]) for r in genes}
    unique_ids = {r["gene_id_stable"] for r in genes}
    region_count = Counter(r.get("ancestry", "UNKNOWN") for r in regions)
    gene_rows_per_group = Counter(r["group_id"] for r in genes)
    mismatches = sorted(groups - set(gene_rows_per_group))
    expected = {"gene_region_pairs": len(genes), "legacy_anchor_only_rows": len(anchor_only), "candidate_evidence_rows": len(candidate_rows), "unique_gene_ids": len(unique_ids)}
    drift = {key: {"source_summary": genes_summary.get(key), "measured": actual} for key, actual in expected.items() if genes_summary.get(key) != actual}
    if drift:
        raise ValueError(f"Source summary/count mismatch: {drift}")
    selected = []
    for gene in sorted(CORE | SECONDARY):
        rr = [r for r in genes if r.get("gene_symbol") == gene]
        selected.append({
            "gene": gene, "legacy_tier": "CORE" if gene in CORE else "SECONDARY",
            "n_expanded_region_rows": len(rr),
            "legacy_anchor_only_rows": sum(r.get("gene_symbol") == gene for r in anchor_only),
            "expanded_regions": sorted({r["group_id"] for r in rr}),
            "max_abf_h4": abf["best_by_gene"].get(gene, {}).get("h4"),
            "abf_tissue": abf["best_by_gene"].get(gene, {}).get("best_tissue"),
            "causal_gene_verified": "NO",
        })
    human_by_gene = {r["target_gene"]: r for r in presence}
    fgf = human_by_gene.get("FGF5", {})
    rna_assay = next((r for r in human_assay if r.get("assay") == "RNA"), {})
    microglia_count = sum(int(r["n_cells"]) for r in cells if r["celltype"] == "Microglia and Macrophages")
    n_cells = sum(int(r["n_cells"]) for r in cells)
    fgf_gap = {
        "gene": "FGF5", "present_in_processed_count_matrix": fgf.get("present"),
        "resolution_method": fgf.get("resolution_method"),
        "matched_features": fgf.get("matched_features"),
        "assay": fgf.get("assay"), "layer": fgf.get("layer"),
        "human_reference": human.get("reference_type"),
        "matrix_gene_count": int(rna_assay["n_features"]) if rna_assay.get("n_features") else None,
        "matrix_cell_count": int(rna_assay["n_cells"]) if rna_assay.get("n_cells") else None,
        "source_rds": human.get("input_rds"),
        "source_rds_sha256": human.get("input_sha256"),
        "technical_root_cause": "UNRESOLVED_WITHOUT_SOURCE_FEATURE_QC",
        "biological_absence_inferred": False,
        "next_checks": [
            "Inspect original Seurat RNA counts rownames and annotations: FGF5 / ENSG00000138675",
            "Audit feature-name mappings and genes removed before GSE256493 distribution",
            "Use independent cerebellum and relevant vessel expression reference; record TPM/detection and cell count",
            "Do not call missing feature a zero-expression result",
        ],
        "microglia_macrophage_cells": microglia_count,
        "all_author_annotated_cells": n_cells,
    }
    discovery = {
        "expanded_distance_components": len(regions),
        "ancestry_label_components": dict(sorted(region_count.items())),
        "provisional_gws_components": sum(r.get("has_gws") == "1" for r in regions),
        "provisional_not_gws_components": sum(r.get("has_gws") != "1" for r in regions),
        "region_gene_rows": len(genes), "candidate_evidence_total_rows": len(candidate_rows),
        "legacy_anchor_only_rows": len(anchor_only), "unique_gene_ids": len(unique_ids),
        "unique_region_gene_pairs": len(pairs), "duplicate_region_gene_rows": len(genes) - len(pairs),
        "missing_gene_region_groups": mismatches,
        "legacy_eas_regions": eas.get("regions"),
        "eur_candidate_components": expansion.get("eur_distance_clusters"),
        "eur_gws_components": expansion.get("eur_gws_clusters"),
        "eur_overlapping_eas_components": expansion.get("eur_overlapping_original_components"),
        "analysis_limitations": [
            "Distance-based components are not independent GWAS loci",
            "Region–gene overlap is not causal prioritization",
            "EAS/EUR cohort and phenotype cannot be pooled without ancestry-aware harmonization",
            "Expanded 2,225 positional gene set has NOT undergone 646 bulk coloc tests",
        ],
    }
    phase11_status = {
        "source_status": "EXISTING_OUTPUT_AUDITED_NOT_INDEPENDENTLY_RECOMPUTED",
        "independent_rerun": bool(phase11.get("independent_rerun", False)),
        "paper_grade_multisignal_coloc": "PENDING",
        "reason": phase11.get("reason_no_independent_rerun"),
    }
    return {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "files": digests,
            "discovery": discovery, "abf": abf, "selected_gene_snapshot": selected,
            "fgf5_feature_gap": fgf_gap, "phase11_gate": phase11_status}

def write_results(payload, out):
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Audit output directory must be empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    (out / "IS_P0_AUDIT.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (out / "IS_P0_SELECTED_GENE_SNAPSHOT.tsv").open("w", newline="", encoding="utf-8") as f:
        columns = ["gene", "legacy_tier", "n_expanded_region_rows", "legacy_anchor_only_rows", "expanded_regions",
                   "max_abf_h4", "abf_tissue", "causal_gene_verified"]
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for r in payload["selected_gene_snapshot"]:
            w.writerow({**r, "expanded_regions": ";".join(r["expanded_regions"])})
    d, a = payload["discovery"], payload["abf"]
    lines = [
        "# IS P0 — Observed source audit",
        "",
        f"- Expanded **distance components**: {d['expanded_distance_components']}; not independent loci.",
        f"- Expanded gene–region rows: {d['region_gene_rows']}; unique positional gene IDs: {d['unique_gene_ids']}.",
        f"- Additional anchor-only rows: {d['legacy_anchor_only_rows']}; candidate table total: {d['candidate_evidence_total_rows']}.",
        f"- GWS-tagged provisional components: {d['provisional_gws_components']}.",
        f"- Prior 4-locus ABF test rows: {a['total_test_rows']} (PASS: {a['pass_test_rows']}).",
        f"- H4 >= 0.5: {a['h4_ge_0_5']}, >= 0.75: {a['h4_ge_0_75']}, >= 0.8: {a['h4_ge_0_8']}.",
        f"- Tested molecular genes: {a['unique_genes_tested']}; tissue datasets: {a['unique_tissue_datasets']}.",
        "- ABF counts are NOT multiple-test-corrected and do NOT cover the expanded gene universe.",
        "- Phase 11B remains exploratory; matched-cohort LD and independent evaluation are pending.",
        "- FGF5: processed RNA counts feature is UNRESOLVED, not established biological nonexpression.",
        "- Core genes are provisional legacy priorities, not demonstrated causal genes.",
        "",
        "## Next gates",
        "1. Freeze full source-input and test-denominator ledger and check missingness.",
        "2. Check FGF5 row names in source Seurat RNA counts and independent reference expression.",
        "3. Run p12/H3/H4 sensitivity, allele and LD provenance QC.",
        "4. Define tier promotions using orthogonal molecular and cell-specific evidence.",
        "5. Reconsider widening molecular-QTL computation only after QC/resource budgeting.",
    ]
    (out / "IS_P0_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="IS data workspace root")
    parser.add_argument("--out", type=Path, required=True, help="New, empty output directory")
    args = parser.parse_args(argv)
    data = audit(args.root)
    write_results(data, args.out)
    print(json.dumps({"out": str(args.out), "discovery": data["discovery"],
                      "abf": {k: v for k,v in data["abf"].items() if k != "best_by_gene"},
                      "fgf5_feature_gap": data["fgf5_feature_gap"]}, indent=2))

if __name__ == "__main__":
    main()
