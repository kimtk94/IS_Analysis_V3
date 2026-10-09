#!/usr/bin/env python3
"""Non-destructive audit of IS LD matrices and actual species-tagged cell evidence."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

LD = "stage3_finemap/japan/bbj/BBJ_LD_MATRIX_AUDIT_V2.tsv"
CELL = "stage5_functional/phase9f_c_target_expression/CELLTYPE_TARGET_EVIDENCE_MASTER.tsv"
GENES = {"FGF5", "ALDH2", "SH3PXD2A", "COL4A2"}

def read_tsv(p):
    if not p.is_file():
        return None
    with p.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def audit(root):
    ld = read_tsv(root / LD)
    cell = read_tsv(root / CELL)
    out = {"version": 1, "read_only": True,
           "ld_source": str(root / LD), "cell_source": str(root / CELL),
           "ld_matrix_audit": {}, "genes": {}, "warnings": []}
    if ld is None:
        out["warnings"].append("LD_AUDIT_FILE_MISSING")
    else:
        for row in ld:
            locus = row.get("locus", "")
            status = str(row.get("status", "UNKNOWN")).upper()
            out["ld_matrix_audit"][locus] = {
                "reported_status": status, "n_variants": row.get("n_variants"),
                "ld_gate": "REVIEW_REQUIRED" if status != "PASS" else "STRUCTURAL_PASS_ONLY",
                "cohort_matched_ld_verified": False}
        out["warnings"].append("STRUCTURAL_LD_PASS_DOES_NOT_PROVE_COHORT_MATCH")
    if cell is None:
        out["warnings"].append("CELL_EVIDENCE_FILE_MISSING")
    for gene in sorted(GENES):
        rows = [r for r in (cell or []) if str(r.get("human_gene", "")).upper() == gene]
        species = Counter(r.get("species", "UNKNOWN") for r in rows)
        datasets = sorted({r.get("dataset", "") for r in rows if r.get("dataset")})
        out["genes"][gene] = {
            "evidence_rows": len(rows), "species_counts": dict(species),
            "datasets": datasets, "top_celltype_labels": sorted({r.get("top_celltype", "") for r in rows if r.get("top_celltype")}),
            "human_validation": "DESCRIPTIVE_ONLY" if any(s.lower() in ("human", "homo sapiens") for s in species) else "NOT_ESTABLISHED",
            "mouse_evidence": "DESCRIPTIVE_ONLY" if any(s.lower() in ("mouse", "mus musculus") for s in species) else "NOT_ESTABLISHED",
            "causal_cell_mechanism": "NOT_ESTABLISHED"}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("/srv/is-analysis/results/is"))
    p.add_argument("--output", type=Path, required=True)
    a=p.parse_args()
    result=audit(a.root)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({"ld_loci": len(result["ld_matrix_audit"]), "genes": result["genes"], "warnings": result["warnings"]}, ensure_ascii=False))
if __name__ == "__main__":
    main()
