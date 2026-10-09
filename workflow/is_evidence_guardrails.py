#!/usr/bin/env python3
"""Fail-closed IS variant, ancestry-LD and colocalization evidence guardrails.

Additive audit layer: does not modify or reclassify frozen source results.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

DNA = {"A", "C", "G", "T"}
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
FIELDS = ("build", "chr", "pos", "ref", "alt", "effect_allele",
          "other_allele", "beta", "se", "p", "eaf", "ancestry")

def canonical_variant(row):
    """Canonical key is genomic REF/ALT, never signed effect-allele orientation."""
    problems = []
    build = str(row.get("build", "")).upper().replace("HG", "GRCH")
    chrom = str(row.get("chr", "")).upper().removeprefix("CHR")
    ref = str(row.get("ref", "")).upper()
    alt = str(row.get("alt", "")).upper()
    ea = str(row.get("effect_allele", "")).upper()
    oa = str(row.get("other_allele", "")).upper()
    if build not in {"GRCH37", "GRCH38"}: problems.append("UNKNOWN_BUILD")
    if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}: problems.append("INVALID_CHR")
    try:
        pos = int(str(row.get("pos", "")))
        if pos <= 0: raise ValueError()
    except (ValueError, TypeError):
        pos = None
        problems.append("INVALID_POS")
    if not ref or not alt or any(x not in DNA for x in ref + alt) or ref == alt:
        problems.append("INVALID_REF_ALT")
    if ea not in DNA or oa not in DNA or ea == oa:
        problems.append("INVALID_EFFECT_ALLELES")
    if len(ref) == len(alt) == 1 and {ea, oa} != {ref, alt}:
        problems.append("ALLELE_MISMATCH")
    if len(ref) != 1 or len(alt) != 1:
        problems.append("INDEL_REQUIRES_NORMALIZATION")
    elif frozenset((ref, alt)) in PALINDROMIC:
        problems.append("PALINDROMIC_REVIEW")
    try:
        beta = float(row.get("beta", ""))
        se = float(row.get("se", ""))
        pval = float(row.get("p", ""))
        eaf = float(row.get("eaf", ""))
        if not all(map(math.isfinite, (beta, se, pval, eaf))) or se <= 0 or not 0 <= pval <= 1 or not 0 <= eaf <= 1:
            raise ValueError()
    except (ValueError, TypeError):
        problems.append("INVALID_SUMMARY_STATS")
    ancestry = str(row.get("ancestry", "")).upper()
    if ancestry not in {"EAS", "EUR"}: problems.append("UNKNOWN_ANCESTRY")
    key = f"{build}:{chrom}:{pos}:{ref}:{alt}" if pos is not None else ""
    blocking = {"UNKNOWN_BUILD", "INVALID_CHR", "INVALID_POS", "INVALID_REF_ALT",
                "INVALID_EFFECT_ALLELES", "ALLELE_MISMATCH", "INDEL_REQUIRES_NORMALIZATION",
                "PALINDROMIC_REVIEW", "INVALID_SUMMARY_STATS", "UNKNOWN_ANCESTRY"}
    return {"variant_key": key, "qc": "BLOCK" if set(problems) & blocking else "PASS",
            "qc_flags": ";".join(problems) if problems else "NONE",
            "effect_orientation": "REF" if ea == ref else "ALT" if ea == alt else "UNKNOWN"}

def coloc_class(row):
    """Only matched credible-signal coloc may earn strong status."""
    if row.get("finemap_status") == "FAILED": return "FINE_MAPPING_FAILED"
    if row.get("ld_match") != "MATCHED": return "LD_UNRESOLVED"
    if row.get("harmonization_qc") != "PASS": return "HARMONIZATION_BLOCKED"
    try:
        h4 = float(row.get("pp_h4", ""))
        h3 = float(row.get("pp_h3", ""))
    except (ValueError, TypeError): return "UNRESOLVED"
    if not 0 <= h4 <= 1 or not 0 <= h3 <= 1: return "UNRESOLVED"
    if row.get("method") != "SUSIE_SIGNAL": return "ABF_ONLY"
    if row.get("signal_pair_valid") != "YES": return "MULTISIGNAL_UNRESOLVED"
    if h4 >= .8 and h4 > h3: return "STRONG_COLOC"
    if h4 >= .5 and h4 > h3: return "PROBABLE_COLOC"
    return "NO_STRONG_SHARED_SIGNAL"

def ld_audit(row):
    trait = str(row.get("gwas_ancestry", "")).upper()
    ld = str(row.get("ld_ancestry", "")).upper()
    cohort = str(row.get("ld_cohort", "")).strip()
    if not cohort or trait not in {"EAS", "EUR"} or ld != trait:
        return "BLOCK"
    if str(row.get("ld_file_verified", "")).upper() != "YES":
        return "BLOCK"
    return "MATCHED"

def functional_priority(row):
    """Transparent evidence tiers, not a fabricated numeric causal probability."""
    required = ("gene", "mechanism_branch", "gwas_qc", "coloc_class",
                "ld_match", "functional_context", "functional_support")
    if any(not str(row.get(k, "")).strip() for k in required):
        return "INCOMPLETE"
    if row["gwas_qc"] != "PASS" or row["ld_match"] != "MATCHED":
        return "QC_UNRESOLVED"
    if row["functional_support"] not in {"YES", "NO", "PENDING"}:
        return "INCOMPLETE"
    if row["coloc_class"] == "STRONG_COLOC" and row["functional_support"] == "YES":
        return "TIER_1_CONVERGENT"
    if row["coloc_class"] in {"STRONG_COLOC", "PROBABLE_COLOC", "ABF_ONLY"}:
        return "TIER_2_FOLLOWUP"
    return "EXPLORATORY"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["variant", "coloc", "ld", "functional"], required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.input.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows: parser.error("input must contain at least one record")
    if args.mode == "variant":
        missing = set(FIELDS) - set(rows[0])
        if missing: parser.error("missing canonical columns: " + ",".join(sorted(missing)))
        audit = [dict(r, **canonical_variant(r)) for r in rows]
    elif args.mode == "coloc":
        audit = [dict(r, evidence_class=coloc_class(r)) for r in rows]
    elif args.mode == "ld":
        audit = [dict(r, ld_match=ld_audit(r)) for r in rows]
    else:
        audit = [dict(r, evidence_tier=functional_priority(r)) for r in rows]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(audit[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(audit)
    print(json.dumps({"mode": args.mode, "rows": len(audit), "output": str(args.output)}, ensure_ascii=False))

if __name__ == "__main__":
    main()
