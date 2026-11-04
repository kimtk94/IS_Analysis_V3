#!/usr/bin/env python3
"""STRICTLY READ-ONLY research server output inventory.

Intentionally only reports allowlisted aggregate publication artifacts (not
participant data), file existence/count, SHA256, and field completeness.
No GWAS/MR/coloc computations; no modifications or network access.
JSON printed to stdout; missing files remain MISSING, never PASS.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/srv/is-analysis")
REPO = ROOT / "IS_Analysis_V3"
SOURCES = {
    "IS_GENES_2225": ("results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2225_GENE_EVIDENCE_MATRIX.tsv", 2225, ("gene_id_stable",)),
    "IS_PAIRS_2425": ("results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2425_GENE_REGION_EVIDENCE.tsv", 2425, ("gene_id_stable",)),
    "IS_UNIVERSE_MANIFEST": ("results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_GENE_EVIDENCE_MATRIX_MANIFEST.json", None, ()),
    "CKD_9GENE_MATRIX": ("results/ckd/integrated_evidence/CKD_9GENE_EVIDENCE_MATRIX.tsv", 9, ("gene_symbol", "eur_abf_h4", "eas_egfr_abf_h4")),
    "CKD_9GENE_QA": ("results/ckd/integrated_evidence/CKD_9GENE_EVIDENCE_QA.json", None, ()),
    "CKD_STAGE2C_SUSIE": ("results/ckd/stage2c_susie/susie_results/STAGE2C_SUSIE_DEFAULT.tsv", None, ("gene_symbol",)),
    "CKD_STAGE3A_KIDNEY": ("results/ckd/stage3a_kidney/STAGE3A_KIDNEY_EVIDENCE.tsv", None, ("gene_symbol",)),
    "CKD_STAGE3B_CELL": ("results/ckd/stage3b_celltype/STAGE3B_INTEGRATED_EVIDENCE.tsv", None, ("gene_symbol",)),
    "CKD_EAS_ANCHOR": ("results/ckd/stage4_eas/anchor_mr/EAS_ANCHOR_WALD_MR.tsv", None, ("gene_symbol",)),
    "CKD_EAS_EGFR_COLOC": ("results/ckd/stage4_eas/regional/eGFR/coloc_results/EAS_COLOC_DEFAULT.tsv", None, ("gene_symbol",)),
}
GENES = {"ACP1", "CPVL", "F12", "GSTA1", "GSTA3", "HLA-E", "INHBC", "SDCCAG8", "UMOD"}
MISSING_VALUES = {"", "NA", "N/A", "NULL", "NONE", "NAN", "MISSING"}


def sha_file(p):
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def read_aggregates(path, key_columns):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        header = reader.fieldnames or []
        missing = sorted(set(key_columns) - set(header))
        if missing:
            return {"schema": "MISSING_FIELDS", "missing_required_fields": missing}
        count = 0
        found = []
        summary = Counter()
        id_unique = set()
        fields = {
            "eur_abf_h4", "eas_egfr_abf_h4", "eas_bun_abf_h4",
            "eur_susie_max_h4", "top_cell_type", "koges_rsid",
        }
        nonempty = {k: 0 for k in fields if k in header}
        for row in reader:
            count += 1
            if None in row or any(v is None for v in row.values()):
                return {"schema": "MALFORMED_ROW", "row": count}
            gene = (row.get("gene_symbol") or "").strip().upper()
            if gene in GENES:
                found.append(gene)
            ens = (row.get("gene_id_stable") or "").strip()
            if ens:
                id_unique.add(ens)
            for key in nonempty:
                if str(row.get(key, "")).strip().upper() not in MISSING_VALUES:
                    nonempty[key] += 1
            status = (row.get("technical_status") or "").upper()
            if "UNRESOLVED" in status or "FAIL" in status:
                summary["technical_unresolved_records"] += 1
        ret = {
            "schema": "OK", "rows": count,
            "required_fields_present": True,
            "number_unique_ensembl_ids": len(id_unique) if "gene_id_stable" in header else None,
            "number_expected_ckd_gene_names_found": len(set(found)),
            "ckd_9_gene_set_exact": set(found) == GENES and len(found) == 9 if "gene_symbol" in header and count == 9 else None,
            "field_present_counts": nonempty,
            "technical_unresolved_records": summary["technical_unresolved_records"],
        }
        return ret


def audit_one(name, spec):
    rel, expected_rows, required_fields = spec
    path = ROOT / rel
    result = {"file": rel, "required_rows_if_frozen": expected_rows,
              "status": "MISSING"}
    if path.is_symlink():
        result["status"] = "SYMLINK_DISALLOWED"
        return result
    if not path.is_file():
        return result
    st = path.stat()
    result.update({
        "status": "FOUND",
        "bytes": st.st_size,
        "last_modified_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
        "sha256": sha_file(path),
    })
    if path.suffix == ".json":
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            result["json_schema"] = "OBJECT" if isinstance(data, dict) else "OTHER"
            # Never copy unreviewed JSON data or values from source to public logs.
            if isinstance(data, dict):
                result["json_keys"] = sorted(str(k) for k in data.keys())[:30]
        except (OSError, ValueError, UnicodeError):
            result["status"] = "INVALID_JSON"
    else:
        result["aggregate"] = read_aggregates(path, required_fields)
        if result["aggregate"].get("schema") != "OK":
            result["status"] = "SCHEMA_ERROR"
        elif expected_rows is not None and result["aggregate"]["rows"] != expected_rows:
            result["status"] = "COUNT_MISMATCH"
    return result


def git_head():
    if not (REPO / ".git").exists():
        return {"repo": "MISSING"}
    ret = {}
    for key, cmd in [
        ("sha", ["git", "-C", str(REPO), "rev-parse", "HEAD"]),
        ("branch", ["git", "-C", str(REPO), "branch", "--show-current"]),
        ("dirty_count", ["git", "-C", str(REPO), "status", "--porcelain"]),
    ]:
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=8, check=False)
            if p.returncode:
                ret[key] = "UNKNOWN"
            elif key == "dirty_count":
                ret[key] = len([x for x in p.stdout.splitlines() if x])
            else:
                ret[key] = p.stdout.strip()[:100]
        except (OSError, subprocess.TimeoutExpired):
            ret[key] = "UNKNOWN"
    return ret


def main():
    sources = {name: audit_one(name, spec) for name, spec in SOURCES.items()}
    statuses = Counter(z["status"] for z in sources.values())
    result = {
        "kind": "READ_ONLY_AGGREGATE_FILE_INVENTORY",
        "audit_utc": datetime.now(timezone.utc).isoformat(),
        "science_execution_performed": False,
        "scientific_validation_complete": False,
        "no_participant_level_data_collected": True,
        "server_git": git_head(),
        "file_status_counts": dict(statuses),
        "sources": sources,
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
