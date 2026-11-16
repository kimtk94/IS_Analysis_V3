#!/usr/bin/env python3
"""GIANT BMI: reconcile exact duplicate variants and exported P=0 values.

Two bounded, SHA256-verified passes over the corrected original publisher gzip.
Only 9 duplicate candidate groups and aggregate Z diagnostics enter JSON.
No raw GWAS is retained and LDSC inputs are not silently approved.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.request import urlopen

from scripts.ckm_stroke.audit_bmi_full_stream import (
    SOURCE, EXPECTED, MAX_BYTES, BoundedHashReader
)
from scripts.ckm_stroke.audit_bmi_hm3_overlap import get_url, EXPECTED_BMI_SHA256

MAX_ROWS = 12_000_000
MAX_DUP_ROWS = 1000
MAX_ZERO_P = 500
NUCS = set("ACGT")
VALID_CHROMOSOMES = {str(i) for i in range(1, 23)}
L10 = math.log(10)


def signed_beta_se_z(beta: str, se: str) -> float:
    b, s = float(beta), float(se)
    if not (math.isfinite(b) and math.isfinite(s) and s > 0):
        raise ValueError("invalid BETA or SE")
    z = b / s
    if not math.isfinite(z):
        raise ValueError("non-finite BETA/SE Z")
    return z


def minus_log10_two_sided_p_from_z(z: float) -> float:
    """Stable normal-Wald -log10(P) for both small and extreme |Z|.

    Do NOT replace original P=0 with any arbitrary finite P. For large Z,
    use the Gaussian Mills ratio asymptotic; logarithmic output only.
    """
    x = abs(z)
    if not math.isfinite(x):
        raise ValueError("non-finite Z")
    if x < 35:
        p = math.erfc(x / math.sqrt(2))
        return -math.log10(p)
    q = 1.0 / (x * x)
    mills = 1-q+3*q*q-15*q**3+105*q**4-945*q**5
    if mills <= 0:
        raise ValueError("nonpositive Mills ratio approximation")
    log_p = math.log(2)-x*x/2 - math.log(2*math.pi)/2 - math.log(x) + math.log(mills)
    return -log_p/L10


def iter_original_rows(opener=urlopen):
    """Streaming iterator whose generator cleanup closes gzip and HTTP response."""
    with get_url(SOURCE, opener) as response:
        content_length = response.headers.get("Content-Length")
        if not content_length or not content_length.isdecimal() or int(content_length) > MAX_BYTES:
            raise ValueError("full original Content-Length absent or exceeds cap")
        reader = BoundedHashReader(response, MAX_BYTES)
        with gzip.GzipFile(fileobj=reader, mode="rb") as gz:
            heading = gz.readline(8192).decode("utf-8-sig", "strict").rstrip("\r\n").split("\t")
            if tuple(heading) != EXPECTED:
                raise ValueError("original BMI column mapping changed")
            for i, raw in enumerate(gz, start=1):
                if i > MAX_ROWS:
                    raise ValueError("full original exceeds row count cap")
                if len(raw) > 8192:
                    raise ValueError("oversized BMI variant row")
                row = raw.decode("utf-8", "strict").rstrip("\r\n").split("\t")
                if len(row) != len(EXPECTED):
                    raise ValueError("unexpected variant field count")
                yield i, row
        if reader.count != int(content_length):
            raise ValueError("original stream incomplete")
        if reader.sha256.hexdigest() != EXPECTED_BMI_SHA256:
            raise ValueError("original compressed SHA256 mismatch")
        yield 0, ["__SOURCE_COMPLETE_SHA256_PASS__", str(reader.count)]


def analyze(opener=urlopen):
    seen_ids: set[str] = set()
    seen_coords: set[tuple[int, int]] = set()
    dup_ids: set[str] = set()
    dup_coords: set[tuple[int, int]] = set()
    zero_ps: list[dict] = []
    compare_p: list[float] = []
    nums = Counter()
    pass_bytes = []

    for i, row in iter_original_rows(opener):
        if i == 0:
            pass_bytes.append(int(row[1]))
            break
        ch, bp, snp, tested, other, hrs_freq, beta, se, p, n = row
        coord = (int(ch), int(bp))
        if ch not in VALID_CHROMOSOMES or tested not in NUCS or other not in NUCS:
            raise ValueError("unexpected original autosomal SNP schema")
        if not snp or not bp.isdigit() or tested == other:
            raise ValueError("invalid original rsID, BP, or allele")
        if snp in seen_ids:
            dup_ids.add(snp)
            nums["duplicate_rsid_extra_rows"] += 1
        else:
            seen_ids.add(snp)
        if coord in seen_coords:
            dup_coords.add(coord)
            nums["duplicate_coordinate_extra_rows"] += 1
        else:
            seen_coords.add(coord)
        nums["total_rows"] += 1
        z = signed_beta_se_z(beta, se)
        pval = float(p)
        if not math.isfinite(pval) or pval < 0 or pval > 1:
            raise ValueError("P outside 0..1")
        if pval == 0:
            if len(zero_ps) >= MAX_ZERO_P:
                raise ValueError("P=0 record limit exceeded")
            zero_ps.append({
                "SNP": snp, "CHR": coord[0], "POS": coord[1],
                "tested_allele": tested, "other_allele": other,
                "beta": float(beta), "se": float(se), "N": float(n),
                "signed_z_from_beta_se": z,
                "neglog10p_normal_wald_from_rounded_beta_se": minus_log10_two_sided_p_from_z(z),
                "published_P": "0", "derived_P_is_original": False,
            })
        elif i % 1000 == 0:
            # Scaled P can differ from beta/SE if model or rounding differs.
            compare_p.append(
                abs(-math.log10(pval) - minus_log10_two_sided_p_from_z(z))
            )
    if len(pass_bytes) != 1:
        raise ValueError("first original stream verification did not complete")
    if not dup_ids and not dup_coords:
        pairs = []
    else:
        # Second bounded complete original read captures only candidate records,
        # allowing actual BETA/SE/P/N comparisons without retaining all 2.3M rows.
        pairs = []
        for i, row in iter_original_rows(opener):
            if i == 0:
                pass_bytes.append(int(row[1]))
                break
            key = (int(row[0]), int(row[1]))
            if row[2] in dup_ids or key in dup_coords:
                pairs.append({"row_number": i, **dict(zip(EXPECTED, row))})
                if len(pairs) > MAX_DUP_ROWS:
                    raise ValueError("too many duplicate candidate records")
        if len(pass_bytes) != 2:
            raise ValueError("second original stream verification did not complete")
    by_id = defaultdict(list)
    by_coord = defaultdict(list)
    for r in pairs:
        by_id[r["SNP"]].append(r)
        by_coord[(r["CHR"], r["POS"])].append(r)
    # Only groups actually seen more than once are reconciled.
    groups = []
    for snp in sorted(dup_ids):
        entries = by_id[snp]
        if len(entries) < 2:
            raise ValueError("duplicate rsID not recovered in pass two")
        vals = [tuple(e[field] for field in EXPECTED) for e in entries]
        identical = all(x == vals[0] for x in vals[1:])
        numeric_keys = ["BETA", "SE", "P", "N"]
        stat_equal = all(
            all(e[x] == entries[0][x] for x in numeric_keys)
            for e in entries[1:]
        )
        coords_equal = all((e["CHR"], e["POS"]) ==
                           (entries[0]["CHR"], entries[0]["POS"]) for e in entries)
        alleles_equal = all(
            (e["Tested_Allele"], e["Other_Allele"]) ==
            (entries[0]["Tested_Allele"], entries[0]["Other_Allele"])
            for e in entries
        )
        groups.append({
            "snp": snp, "rows": entries,
            "identical_all_fields": identical,
            "identical_beta_se_p_n": stat_equal,
            "identical_chrpos": coords_equal,
            "identical_ordered_alleles": alleles_equal,
            "dedup_decision": "SAFE_TO_DROP_IDENTICAL_REPEAT" if identical else "HOLD_DISCORDANT_REPEAT",
        })
    # A coordinate can be duplicated under distinct rsIDs, too.
    cross_id_coords = [
        {"chr": a, "pos": b, "rsids": sorted({e["SNP"] for e in group})}
        for (a, b), group in sorted(by_coord.items())
        if len(group) > 1 and len({e["SNP"] for e in group}) > 1
    ]
    stats = {
        "samples_checked": len(compare_p),
        "abs_log10p_error_median": (sorted(compare_p)[len(compare_p)//2] if compare_p else None),
        "abs_log10p_error_max": max(compare_p) if compare_p else None,
        "note": "P-based and beta/SE Wald diagnostics on deterministic ~0.1% subset; rounding/model differences possible",
    }
    nums["distinct_duplicate_rsids"] = len(dup_ids)
    nums["distinct_duplicate_coordinates"] = len(dup_coords)
    nums["p_underflow_zero_rows"] = len(zero_ps)
    nums["two_complete_sha256_matched_passes"] = len(pass_bytes)
    discordant = [g["snp"] for g in groups if g["dedup_decision"] != "SAFE_TO_DROP_IDENTICAL_REPEAT"]
    neglogs = [r["neglog10p_normal_wald_from_rounded_beta_se"] for r in zero_ps]
    return {
        "schema": "ckm-giant-bmi-duplicates-zero-p-v1",
        "at_utc": datetime.now(timezone.utc).isoformat(),
        "source": SOURCE,
        "source_sha256": EXPECTED_BMI_SHA256,
        "verification": {"pass_bytes": pass_bytes, "source_gzip_retained": False},
        "counts": dict(nums),
        "duplicate_groups": groups,
        "duplicate_different_rsids_at_same_coordinate": cross_id_coords,
        "discordant_duplicate_rsids": discordant,
        "zero_p_records": zero_ps,
        "zero_p_derived_neglog10p_range": [
            min(neglogs), max(neglogs)
        ] if neglogs else None,
        "nonzero_p_wald_comparison": stats,
        "dedup_ready_for_downstream": not discordant and not cross_id_coords,
        "zero_P_policy": "RETAIN_ORIGINAL_ZERO; INDEPENDENT_SIGNED_BETA_OVER_SE_WALD_LOGP_FOR_DIAGNOSTICS_ONLY; DO_NOT_IMPUTE_P",
        "final_sumstats_written": False,
        "ldsc_h2_rg_performed": False,
        "scientific_status": "HOLD",
        "ready_for_ldsc": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "technical_status": "HOLD_SOURCE_OR_VARIANT_QC",
        "scientific_status": "HOLD", "ready_for_ldsc": False
    }
    try:
        report = analyze()
        report["technical_status"] = (
            "PASS_DUPLICATE_EXACT_AND_ZERO_P_DIAGNOSTICS"
            if report["dedup_ready_for_downstream"]
            else "HOLD_DISCORDANT_DUPLICATE_OR_COORDINATE_COLLISION"
        )
    except Exception as exc:
        report["error_type"] = type(exc).__name__
        report["error_detail"] = str(exc)[:200]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "technical_status": report["technical_status"],
        "counts": report.get("counts"),
        "discordant_ids": report.get("discordant_duplicate_rsids"),
        "zero_p_neglog10_range": report.get("zero_p_derived_neglog10p_range"),
        "scientific_status": report["scientific_status"],
    }))
    if not report["technical_status"].startswith("PASS_DIAGNOSTICS"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
