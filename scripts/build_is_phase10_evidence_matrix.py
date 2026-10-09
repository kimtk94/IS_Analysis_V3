#!/usr/bin/env python3
"""Build conservative IS phase 10 evidence matrix from source status, no inferred validation."""
import argparse
import csv
import json
from pathlib import Path
from scripts.audit_is_phase10_readiness import audit

GENES = [("FGF5", "BBJ_IS_L001", "regulatory_protein"),
         ("ALDH2", "BBJ_IS_L003", "coding_metabolic"),
         ("SH3PXD2A", "BBJ_IS_L002", "immune_vascular"),
         ("COL4A2", "BBJ_IS_L004", "mural_vascular")]
HANDOFF = "stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_COLAB_R3_HANDOFF.tsv"
DRIVE = "stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_DRIVE_HANDOFF.tsv"

def pairs(path):
    if not path.is_file(): return {}
    with path.open(newline="", encoding="utf-8") as f:
        return {r["component"]: r["status"] for r in csv.DictReader(f, delimiter="\t")}

def build(root):
    a = audit(root)
    h = pairs(root / HANDOFF)
    d = pairs(root / DRIVE)
    run_status = h.get("STATUS", "MISSING")
    human_gate = "AWAITING_COLAB_OUTPUT_VALIDATION"
    if run_status != "READY_TO_RUN_NOT_EXECUTED":
        human_gate = "REVIEW_HANDOFF_STATUS"
    result = []
    for gene, locus, mechanism in GENES:
        evidence = a["genes"][gene]
        ld = a["ld_matrix_audit"].get(locus, {})
        result.append({
            "gene": gene, "locus": locus, "mechanism": mechanism,
            "mouse_rows": evidence["species_counts"].get("Mouse", 0),
            "mouse_celltypes_descriptive": ";".join(evidence["top_celltype_labels"]),
            "mouse_evidence": evidence["mouse_evidence"],
            "human_rds_local": d.get("LOCAL_RDS", "UNKNOWN"),
            "human_rds_drive": d.get("REMOTE_RDS", "UNKNOWN"),
            "human_colab_status": run_status,
            "human_celltype_validation": "NOT_ESTABLISHED",
            "ld_matrix_status": ld.get("reported_status", "NOT_AUDITED"),
            "cohort_matched_ld_verified": "NO",
            "phase10_gate": human_gate,
            "final_causal_claim": "NOT_ESTABLISHED",
        })
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("/srv/is-analysis/results/is"))
    p.add_argument("--output", type=Path, required=True)
    a=p.parse_args()
    rows=build(a.root)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader();w.writerows(rows)
    print(json.dumps({"genes":len(rows), "human_validated":sum(r["human_celltype_validation"]=="VALIDATED" for r in rows), "ld_missing":sum(r["ld_matrix_status"]=="MISSING" for r in rows), "output":str(a.output)}))
if __name__=="__main__": main()
