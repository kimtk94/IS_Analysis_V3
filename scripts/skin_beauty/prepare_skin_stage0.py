#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, json, math, re
from pathlib import Path
from typing import Optional

OUTCOMES = {
    "ebi-a-GCST90094903": "facial_wrinkles_under_eye",
    "ebi-a-GCST90094904": "facial_wrinkles_crows_feet",
}

CANONICAL_ST16 = Path("/srv/is-analysis/data/ckd/analysis_ready/UKBPPP_ST16_cis_independent.tsv.gz")

PATTERNS = [
    "*UKBPPP*ST16*", "*UKB*PPP*ST16*", "*UKBPPP*cis*", "*UKB*PPP*cis*",
    "*pQTL*cis*", "*pqtl*cis*", "*UKBPPP*ST9*", "*UKB*PPP*ST9*",
]

ALIASES = {
    "protein": ["protein_id", "protein", "aptamer", "assay", "gene_symbol", "gene"],
    "variant": ["rsid", "rs_id", "snp", "variant", "variant_id", "variant_id_hg19", "markername", "id"],
    "chr": ["chr_hg19", "chromosome", "chr", "chrom"],
    "pos": ["pos_hg19", "base_pair_location", "position", "pos", "bp"],
    "ea": ["effect_allele", "ea", "alt", "a1", "tested_allele"],
    "oa": ["other_allele", "oa", "ref", "a2", "non_effect_allele"],
    "beta": ["beta_exposure", "beta", "effect", "b"],
    "se": ["se_exposure", "standard_error", "se", "stderr"],
    "p": ["p_exposure", "p_value", "pval", "p", "pvalue"],
    "mlogp": ["minus_log10p_exposure", "minus_log10_p", "mlog10p", "log10p"],
    "eaf": ["eaf", "effect_allele_frequency", "af", "freq", "a1freq"],
    "n": ["n", "n_total", "samplesize", "sample_size"],
    "f": ["f_stat", "fstat", "f_statistic"],
}

def norm(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", x.strip().lower()).strip("_")

def open_text(p: Path):
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.suffix == ".gz" else p.open("rt", encoding="utf-8", errors="replace")

def detect(cols):
    m = {norm(c): c for c in cols}
    out = {}
    for role, aliases in ALIASES.items():
        out[role] = next((m[a] for a in aliases if a in m), None)
    return out

def safe_float(x):
    try:
        y = float(x)
        return y if math.isfinite(y) else None
    except Exception:
        return None

def exposure_score(p: Path) -> int:
    rp = str(p.resolve()) if p.exists() else str(p)
    s = (str(p) + " " + rp).lower()
    score = 0
    if rp == str(CANONICAL_ST16): score += 100000
    if "st16" in s: score += 10000
    if "cis_independent" in s or "cis-independent" in s: score += 5000
    if "harmonization_safe" in s: score += 600
    if "primary_non_mhc" in s: score += 500
    if "primary" in s: score += 300
    if "canonical" in s: score += 250
    if "st9" in s: score -= 2000
    if "strongest_cis_per_protein" in s: score -= 1000
    if "marginal" in s: score -= 500
    return score

def collect_exposures(search_root: Path):
    raw = []
    if CANONICAL_ST16.exists(): raw.append(CANONICAL_ST16)
    for pat in PATTERNS:
        try:
            raw.extend(p for p in search_root.rglob(pat) if p.is_file())
        except PermissionError:
            pass
    # Include the project symlink explicitly even if rglob does not follow it.
    link = Path("/srv/is-analysis/IS_Analysis_V3/data/skin_beauty/stage0_gwas/UKBPPP_CIS_INSTRUMENTS.tsv.gz")
    if link.exists(): raw.append(link)

    # De-duplicate by displayed path, but retain symlink and target as separate audit rows.
    seen, out = set(), []
    for p in raw:
        sp = str(p)
        if sp not in seen:
            seen.add(sp); out.append(p)
    return sorted(out, key=lambda p: (-exposure_score(p), str(p)))

def audit_table(path: Path, max_rows=5_000_000):
    resolved = path.resolve()
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            return {"path": str(path), "resolved_path": str(resolved), "error": "empty_file"}
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        cols = first.rstrip("\r\n").split(delim)
        cmap = detect(cols)
        idx = {k: (cols.index(v) if v in cols else None) for k,v in cmap.items()}
        nrow = invalid_eaf = invalid_allele = dup = 0
        pmin = None; eaf_min = None; eaf_max = None
        seen_variants = set()
        for line in fh:
            if nrow >= max_rows: break
            vals = line.rstrip("\r\n").split(delim); nrow += 1
            p = None
            if idx["p"] is not None and idx["p"] < len(vals): p = safe_float(vals[idx["p"]])
            elif idx["mlogp"] is not None and idx["mlogp"] < len(vals):
                z = safe_float(vals[idx["mlogp"]]); p = 10**(-z) if z is not None and z < 310 else (0.0 if z is not None else None)
            if p is not None and 0 <= p <= 1: pmin = p if pmin is None else min(pmin, p)
            if idx["eaf"] is not None and idx["eaf"] < len(vals):
                e = safe_float(vals[idx["eaf"]])
                if e is None or e < 0 or e > 1: invalid_eaf += 1
                else:
                    eaf_min = e if eaf_min is None else min(eaf_min, e)
                    eaf_max = e if eaf_max is None else max(eaf_max, e)
            valid = {"A","C","G","T"}
            for role in ("ea","oa"):
                if idx[role] is not None and idx[role] < len(vals):
                    a = vals[idx[role]].upper()
                    if a and a not in valid: invalid_allele += 1
            if idx["variant"] is not None and idx["variant"] < len(vals) and len(seen_variants) < 1_000_000:
                v = vals[idx["variant"]]
                if v in seen_variants: dup += 1
                else: seen_variants.add(v)
    return {
        "path": str(path), "resolved_path": str(resolved), "is_symlink": path.is_symlink(),
        "size_bytes": resolved.stat().st_size if resolved.exists() else path.stat().st_size,
        "score": exposure_score(path), "columns": cols, "canonical_columns": cmap,
        "rows_scanned": nrow, "scan_truncated": nrow >= max_rows, "min_p": pmin,
        "eaf_min": eaf_min, "eaf_max": eaf_max, "invalid_eaf": invalid_eaf,
        "non_acgt_allele_cells": invalid_allele, "duplicate_variant_keys_first_million": dup,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/srv/is-analysis/IS_Analysis_V3")
    ap.add_argument("--search-root", default="/srv/is-analysis")
    ap.add_argument("--data-dir", default="/srv/is-analysis/IS_Analysis_V3/data/skin_beauty/stage0_gwas")
    ap.add_argument("--out-dir", default="/srv/is-analysis/IS_Analysis_V3/results/skin_beauty/stage0_audit")
    args = ap.parse_args()
    search_root = Path(args.search_root); data_dir = Path(args.data_dir); out_dir = Path(args.out_dir)
    data_dir.mkdir(parents=True, exist_ok=True); out_dir.mkdir(parents=True, exist_ok=True)

    candidates = collect_exposures(search_root)
    preferred = candidates[0] if candidates else None
    audits = []
    for p in candidates[:25]:
        try: audits.append(audit_table(p, max_rows=1_000_000))
        except Exception as e: audits.append({"path": str(p), "error": repr(e), "score": exposure_score(p)})

    report = {
        "project":"skin_beauty_mr", "stage":"stage0_audit", "outcomes":OUTCOMES,
        "search_roots":[str(search_root)],
        "source_priority":["OpenGWAS ebi-a mirror", "NHGRI-EBI accession metadata"],
        "canonical_st16": str(CANONICAL_ST16), "canonical_st16_exists": CANONICAL_ST16.exists(),
        "exposure_candidates":[str(p) for p in candidates],
        "preferred_exposure": str(preferred) if preferred else None,
        "preferred_exposure_audit": audits[0] if audits else None,
        "exposure_audits": audits,
        "outcome_files": [],
    }
    (out_dir/"STAGE0_SKIN_GWAS_AUDIT.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    with (out_dir/"STAGE0_EXPOSURE_FILE_CANDIDATES.tsv").open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f, delimiter="\t"); w.writerow(["rank","score","path","resolved_path"])
        for i,p in enumerate(candidates,1): w.writerow([i, exposure_score(p), str(p), str(p.resolve())])

    print("===== SKIN BEAUTY STAGE 0 AUDIT =====")
    print("root      =", args.root); print("search    =", search_root); print("data_dir  =", data_dir); print("out_dir   =", out_dir)
    print("canonical ST16 exists =", CANONICAL_ST16.exists())
    print("exposure candidates =", len(candidates)); print("preferred exposure   =", preferred)
    if audits:
        print("preferred rows scanned=", audits[0].get("rows_scanned")); print("preferred columns     =", ",".join(audits[0].get("columns",[])[:20]))
    return 0

if __name__ == "__main__": raise SystemExit(main())
