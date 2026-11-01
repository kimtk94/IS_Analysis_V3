#!/usr/bin/env python3
"""P3: complete, public GSE313736 matrix integrity by donor pool.

Download to transient runner only; upload hashes and aggregate counts, never
expression matrices, private genotypes, or donor identifiers.
"""
import argparse
import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from ocular_p1_integrity import BASE, download, metadata_summary, paired_probe

BASELINE = {"A": 8, "B": 7, "C": 7, "D": 7, "E": 7, "F": 6,
            "G": 7, "H": 6, "I": 4}
P1_C = {
    "vehicle_barcodes": 16980,
    "dex_barcodes": 22011,
    "vehicle_nnz": 92732638,
    "dex_nnz": 100908670,
}


def audit_pool(pool):
    idx = ord(pool) - ord("A")
    out = {
        "pool": pool,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "GSE313736",
        "status": "HOLD",
        "scope": "in-vitro DEX response; NOT clinical delta-IOP",
    }
    with tempfile.TemporaryDirectory(prefix="ocular-p3-") as td:
        tmp = Path(td)
        mfile = tmp / "cell_metadata.tsv.gz"
        try:
            out["metadata_transfer"] = download(
                BASE + "series/GSE313nnn/GSE313736/suppl/"
                "GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
                mfile,
                12 * 1024 * 1024,
            )
            meta = metadata_summary(mfile)
        except Exception as exc:
            out["error"] = f"METADATA: {type(exc).__name__}: {str(exc)[:160]}"
            return out
        donors = {}
        out["samples"] = {}
        for condition, gsm, suffix in (
            ("vehicle", f"GSM{9373276 + idx}", f"pool_{pool}"),
            ("dex", f"GSM{9373285 + idx}", f"pool_dex{pool}"),
        ):
            entry, donors[condition] = paired_probe(
                condition, gsm, suffix, meta, tmp
            )
            out["samples"][condition] = entry
        out["aggregate"] = {
            "vehicle_donors": len(donors.get("vehicle", set())),
            "dex_donors": len(donors.get("dex", set())),
            "paired_donors": len(
                donors.get("vehicle", set()) & donors.get("dex", set())
            ),
            "expected_paired_from_p2": BASELINE[pool],
        }
    valid = (
        all(
            x["status"] == "PASS_MATRIX_AND_BARCODE_JOIN"
            for x in out["samples"].values()
        )
        and out["aggregate"]["paired_donors"] == BASELINE[pool]
    )
    out["status"] = "PASS_FULL_TWO_MATRIX_INTEGRITY" if valid else "HOLD"
    return out


def aggregate(paths):
    reports = [json.loads(Path(p).read_text()) for p in paths]
    pools = [d.get("pool") for d in reports]
    expected = [p for p in BASELINE if p != "C"]
    errors = []
    if len(pools) != len(expected) or set(pools) != set(expected):
        errors.append("EXPECTED_EXACTLY_EIGHT_UNIQUE_NON_C_POOLS")
    for d in reports:
        p = d.get("pool")
        if d.get("status") != "PASS_FULL_TWO_MATRIX_INTEGRITY":
            errors.append(f"{p}_NOT_VALIDATED")
        if p in BASELINE and p != "C":
            x = d.get("aggregate", {})
            if x.get("paired_donors") != BASELINE[p]:
                errors.append(f"{p}_DONOR_COUNT_MISMATCH")
    combined = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS_18_MATRICES_INTEGRITY" if not errors else "HOLD",
        "scope": "16 new matrices plus 2 previously validated C matrices",
        "new_pools": sorted(pools),
        "new_matrices_passed": sum(
            x.get("status") == "PASS_MATRIX_AND_BARCODE_JOIN"
            for d in reports for x in d.get("samples", {}).values()
        ),
        "previous_pool_C_matrix_passed": 2,
        "new_paired_donor_pool_counts_sum": sum(
            d.get("aggregate", {}).get("paired_donors", 0) for d in reports
        ),
        "combined_expected_paired_donors_from_P2": sum(BASELINE.values()),
        "errors": errors,
        "notes": [
            "Donor IDs are not written to artifacts.",
            "Total unique donor count (59) comes from P2, not the sum alone.",
            "Matrix integrity is NOT donor-by-celltype pseudobulk readiness.",
            "Remaining: cell-state pseudobulk QC, independent replication, novelty.",
        ],
    }
    return combined


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", choices=list(BASELINE))
    ap.add_argument("--aggregate", nargs="*", metavar="REPORT")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if bool(args.pool) == (args.aggregate is not None):
        ap.error("specify either --pool or --aggregate")
    out = audit_pool(args.pool) if args.pool else aggregate(args.aggregate)
    p = Path(args.output)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    print(
        "OCULAR_P3_RESULT",
        out["status"],
        "pool" if args.pool else "aggregate",
        args.pool or len(args.aggregate),
    )
    if args.pool:
        for name, row in out.get("samples", {}).items():
            print(
                "MATRIX", name, row.get("status"),
                row.get("matrix", {}), "error", row.get("reason"),
            )
    else:
        print("COUNTS", out["new_matrices_passed"],
              out["new_paired_donor_pool_counts_sum"])
        print("ISSUES", out["errors"])
    if out["status"] not in (
        "PASS_FULL_TWO_MATRIX_INTEGRITY",
        "PASS_18_MATRICES_INTEGRITY",
    ):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
