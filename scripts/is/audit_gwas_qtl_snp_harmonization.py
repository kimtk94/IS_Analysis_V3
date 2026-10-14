#!/usr/bin/env python3
"""Fail-closed SNP allele harmonization audit; no LD/coloc inference.

Inputs are explicit normalized variant tables (TSV) with columns:
chrom, pos, ref, alt, effect_allele, other_allele, beta, pvalue.
Both inputs MUST be the same documented genome build. This script never liftover.
"""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

COLUMNS = ("chrom", "pos", "ref", "alt", "effect_allele", "other_allele", "beta", "pvalue")
COMP = str.maketrans("ACGT", "TGCA")

def read_variants(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or any(c not in reader.fieldnames for c in COLUMNS):
            raise ValueError(f"Required columns missing: {path}")
        records = {}
        for line, r in enumerate(reader, 2):
            if None in r or any(r[c] is None for c in COLUMNS):
                raise ValueError(f"Malformed row {line} in {path}")
            chrom = r["chrom"].removeprefix("chr").upper()
            pos = int(r["pos"])
            ref, alt = r["ref"].upper(), r["alt"].upper()
            ea, oa = r["effect_allele"].upper(), r["other_allele"].upper()
            if pos <= 0 or chrom not in [str(i) for i in range(1, 23)] + ["X", "Y"] or ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt:
                raise ValueError(f"Unsupported or invalid SNP at row {line}")
            if {ea, oa} != {ref, alt}:
                raise ValueError(f"Effect allele inconsistency at row {line}")
            beta, pvalue = float(r["beta"]), float(r["pvalue"])
            if not math.isfinite(beta) or not math.isfinite(pvalue) or not 0 < pvalue <= 1:
                raise ValueError(f"Invalid beta/pvalue at row {line}")
            key = (chrom, pos)
            if key in records:
                raise ValueError(f"Duplicate position {key} in {path}")
            records[key] = (ref, alt, ea, oa, beta, pvalue)
        return records

def classify(g, q):
    gr, ga, ge, _, _, _ = g
    qr, qa, qe, _, _, _ = q
    if {gr, ga} == {qr, qa}:
        if {gr, ga} in ({"A", "T"}, {"C", "G"}):
            return "PALINDROMIC_UNRESOLVED", None
        return "MATCHED" if ge == qe else "SWAPPED", (1 if ge == qe else -1)
    if {gr, ga} == {qr.translate(COMP), qa.translate(COMP)}:
        # Strand flips are flagged for external reference verification, not automatically resolved.
        return "STRAND_UNRESOLVED", None
    return "MISMATCH", None

def audit(gwas, qtl, gwas_build, qtl_build, out):
    if gwas_build != qtl_build or gwas_build not in ("GRCh37", "GRCh38"):
        raise ValueError("Genome builds must be identical and explicit")
    g = read_variants(gwas)
    q = read_variants(qtl)
    if not g or not q:
        raise ValueError("GWAS and QTL must each contain at least one SNP")
    counts = Counter()
    rows = []
    for chrom, pos in sorted(g.keys() & q.keys()):
        status, multiplier = classify(g[(chrom, pos)], q[(chrom, pos)])
        counts[status] += 1
        rows.append({"chrom": chrom, "pos": pos, "status": status,
                     "gwas_effect_allele": g[(chrom, pos)][2],
                     "qtl_effect_allele": q[(chrom, pos)][2],
                     "qtl_beta_aligned_to_gwas": (
                         f"{q[(chrom, pos)][4] * multiplier:.12g}" if multiplier is not None else "")})
    counts["GWAS_ONLY"] = len(g.keys() - q.keys())
    counts["QTL_ONLY"] = len(q.keys() - g.keys())
    output = Path(out)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "ALLELE_HARMONIZATION_AUDIT.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chrom", "pos", "status", "gwas_effect_allele",
                         "qtl_effect_allele", "qtl_beta_aligned_to_gwas"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    report = {"status": "ALLELE_AUDIT_ONLY", "genome_build": gwas_build,
              "n_gwas": len(g), "n_qtl": len(q), "n_overlapping_positions": len(rows),
              "n_usable_nonpalindromic": counts["MATCHED"] + counts["SWAPPED"],
              "counts": dict(counts), "ld_verified": False, "coloc_recomputed": False}
    (output / "ALLELE_HARMONIZATION_AUDIT.json").write_text(json.dumps(report, indent=2) + "\n")
    return report

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--gwas", type=Path, required=True)
    p.add_argument("--qtl", type=Path, required=True)
    p.add_argument("--gwas-build", required=True)
    p.add_argument("--qtl-build", required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(audit(a.gwas, a.qtl, a.gwas_build, a.qtl_build, a.out), indent=2))
