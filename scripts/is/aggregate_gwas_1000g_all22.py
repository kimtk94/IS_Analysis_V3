#!/usr/bin/env python3
"""Cross-ancestry 22 autosome GWAS vs population-matched 1000G BIM QC.

Only compares technical source/reference overlap, NOT EUR/EAS biological
stroke differences. Requires both original per-population GitHub artifacts
with original publisher-MD5 check, 22/22 regression, and all fractions.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys

EXPECTED_GWAS_ROWS={"EUR":7482032,"EAS":6725216}
CATEGORIES=(
    "UNORDERED_ALLELE_SET_EXACT",
    "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION",
    "STRAND_COMPLEMENT_UNRESOLVED",
    "COORDINATE_NOT_IN_LD_BIM",
    "MULTI_ALLELIC_REFERENCE_POSITION_REVIEW",
    "COORDINATE_PRESENT_ALLELES_DISCORDANT",
    "UNSUPPORTED_NON_SNP_OR_ALLELE",
)


def valid_source(path,ancestry):
    if not path.is_file():
        raise ValueError("Missing true-GWAS source evidence artifact: "+str(path))
    d=json.loads(path.read_text(encoding="utf-8"))
    if (d.get("ancestry")!=ancestry or
            d.get("status")!="ALL_22_AUTOSOME_ORIGINAL_MD5_AND_BIM_OVERLAP_PASS" or
            d.get("all_22_complete") is not True or
            d.get("previous_chr1_2_regression",{}).get("status")!="PASS" or
            d.get("GWAS_checksum",{}).get("publisher_MD5_verified") is not True or
            d.get("LD_checksum",{}).get("publisher_MD5_verified") is not True or
            d.get("coloc_ready") is not False or
            d.get("allele_effect_direction_harmonized") is not False or
            d.get("LD_matrix_computed") is not False):
        raise ValueError("Unverified or overstated 22 chr source QC: "+ancestry)
    by=d.get("by_chr",{})
    if set(by)!={str(x) for x in range(1,23)}:
        raise ValueError("Incomplete chromosome universe: "+ancestry)
    rows=[]
    aggregate=Counter()
    for chrstr in sorted(by,key=int):
        r=by[chrstr]
        n=int(r["gwas_rows"])
        counts={key:int(r.get(key,0)) for key in CATEGORIES}
        if n<100 or sum(counts.values())!=n or any(v<0 for v in counts.values()):
            raise ValueError("Unreconciled variant categories on chromosome "+chrstr)
        coord=sum(v for k,v in counts.items() if k not in
                  ("COORDINATE_NOT_IN_LD_BIM","UNSUPPORTED_NON_SNP_OR_ALLELE"))
        if coord!=int(r["chr_pos_overlap"]):
            raise ValueError("Coordinate-level denominator discrepancy "+chrstr)
        if abs(100*coord/n-float(r["chr_pos_overlap_pct"]))>0.000051:
            raise ValueError("Claimed reference overlap differs from counts "+chrstr)
        rows.append({
            "ancestry":ancestry,
            "chromosome":int(chrstr),
            "gwas_rows":n,
            "reference_bim_rows":int(r["bim_rows"]),
            "coordinate_overlap_n":coord,
            "coordinate_overlap_pct":round(100*coord/n,4),
            "nonpalindromic_allele_match_n":counts["UNORDERED_ALLELE_SET_EXACT"],
            "nonpalindromic_allele_match_pct":round(100*counts["UNORDERED_ALLELE_SET_EXACT"]/n,4),
            "palindromic_ambiguous_n":counts["PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"],
            "coordinate_missing_n":counts["COORDINATE_NOT_IN_LD_BIM"],
            "coordinate_missing_pct":round(100*counts["COORDINATE_NOT_IN_LD_BIM"]/n,4),
            "strand_complement_n":counts["STRAND_COMPLEMENT_UNRESOLVED"],
            "allele_discordant_n":counts["COORDINATE_PRESENT_ALLELES_DISCORDANT"],
            "unsupported_SNP_n":counts["UNSUPPORTED_NON_SNP_OR_ALLELE"],
            "source_status":"PUBLISHER_MD5_AND_22CHR_COORDINATE_VERIFIED",
        })
        aggregate["rows"]+=n
        aggregate["coordinate_found"]+=coord
        aggregate["allele_nonpalindromic"]+=counts["UNORDERED_ALLELE_SET_EXACT"]
        aggregate["palindromic"]+=counts["PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"]
        aggregate["reference_absent"]+=counts["COORDINATE_NOT_IN_LD_BIM"]
    if aggregate["rows"]!=EXPECTED_GWAS_ROWS[ancestry]:
        raise ValueError("Full original GWAS autosomal row totals disagree with publisher QC")
    if int(d["overall"]["all_22_autosome_gwas_rows"])!=aggregate["rows"] or int(d["overall"]["all_22_chr_pos_overlap"])!=aggregate["coordinate_found"]:
        raise ValueError("Source report aggregate differs from independently summed chromosomes")
    return rows,dict(aggregate),d


def summarize(eur,eas):
    indexed={(x["ancestry"],x["chromosome"]):x for x in eur+eas}
    comparison=[]
    for c in range(1,23):
        e=indexed[("EUR",c)]
        a=indexed[("EAS",c)]
        comparison.append({
            "chromosome":c,
            "EUR_GWAS_rows":e["gwas_rows"],
            "EAS_GWAS_rows":a["gwas_rows"],
            "EUR_BIM_coordinate_coverage_pct":e["coordinate_overlap_pct"],
            "EAS_BIM_coordinate_coverage_pct":a["coordinate_overlap_pct"],
            "EAS_minus_EUR_coverage_percentage_points":round(
                a["coordinate_overlap_pct"]-e["coordinate_overlap_pct"],4),
            "EUR_absent_coordinate_pct":e["coordinate_missing_pct"],
            "EAS_absent_coordinate_pct":a["coordinate_missing_pct"],
            "phenotype_ancestry_differences_inferred":False,
        })
    return comparison


def write_csv(path,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(root,out):
    eur,ec,ei=valid_source(root/"is-GWAS-1000G-all22-EUR"/"all22_audit.json","EUR")
    eas,ac,ai=valid_source(root/"is-GWAS-1000G-all22-EAS"/"all22_audit.json","EAS")
    comparisons=summarize(eur,eas)
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/"all22_ancestry_by_chromosome.tsv",eur+eas)
    write_csv(out/"all22_EAS_vs_EUR_reference_coverage.tsv",comparisons)
    status={
        "as_of_utc":datetime.now(timezone.utc).isoformat(),
        "source_GitHub_run":43031787049,
        "EUR_publisher_and_LD_MD5":"PASS",
        "EAS_publisher_and_LD_MD5":"PASS",
        "all_22_autosome_source_QC":"PASS",
        "prior_chr1_chr2_replication":"PASS_BOTH",
        "EUR_original_rows":ec["rows"],
        "EAS_original_rows":ac["rows"],
        "EUR_coordinate_overlap_pct":round(100*ec["coordinate_found"]/ec["rows"],4),
        "EAS_coordinate_overlap_pct":round(100*ac["coordinate_found"]/ac["rows"],4),
        "EUR_palindromic_unresolved":ec["palindromic"],
        "EAS_palindromic_unresolved":ac["palindromic"],
        "reference_allele_frequency_alignment":"NOT_RUN",
        "association_effect_orientation":"NOT_HARMONIZED",
        "LD_correlation_matrices":"NOT_RUN",
        "genome_build_independent_validation":"NOT_VERIFIED",
        "per_variant_study_sample_sizes":"NOT_AVAILABLE",
        "EAS_full_cis_eQTL":"HOLD",
        "MAGMA":"NOT_RUN",
        "LDSC_SEG":"NOT_RUN",
        "coloc":"NOT_RUN",
        "status":"REFERENCE_POSITION_AND_ALLELE_CLASSIFICATION_PASS_NOT_EVIDENCE_OF_DISEASE_MECHANISM",
    }
    (out/"technical_qc_decision.json").write_text(json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=[
        "# EUR/EAS IS original GWAS vs population-matched 1000G Phase3: 22-chromosome source summary","",
        "**All data are original GWAS/reference-source technical QC, not variant-based etiologic inference.**","",
        "| Source audit | EUR | EAS |",
        "|---|---:|---:|",
        f"| GWAS rows | {ec['rows']:,} | {ac['rows']:,} |",
        f"| LD reference coordinate overlap | {status['EUR_coordinate_overlap_pct']}% | {status['EAS_coordinate_overlap_pct']}% |",
        f"| Palindromic unresolved | {ec['palindromic']:,} | {ac['palindromic']:,} |",
        f"| LD reference position absent | {ec['reference_absent']:,} | {ac['reference_absent']:,} |",
        "",
        "## Highest per-chromosome missing-coordinate rates (within each ancestry)",
    ]
    for anc,rows in (("EUR",eur),("EAS",eas)):
        top=sorted(rows,key=lambda x:-x["coordinate_missing_pct"])[:4]
        lines.extend([f"- {anc}: "+", ".join(f"chr{x['chromosome']} {x['coordinate_missing_pct']:.2f}% missing" for x in top)])
    lines.extend(["",
        "Different LD panel variant coverage does not equal different stroke biology.",
        "Effect allele direction, palindromic frequency, independent genome build, per-variant N, EAS full cis-eQTL and causality remain unverified.",
    ])
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(md)
    return status


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",default="all22_source")
    p.add_argument("--out-dir",default="artifacts/is_gwas_ld_all22_synthesis")
    a=p.parse_args()
    run(Path(a.source),Path(a.out_dir))


if __name__=="__main__":
    try:
        main()
    except (ValueError,KeyError,OSError,TypeError) as e:
        print("FAIL_CLOSED: "+str(e),file=sys.stderr)
        raise
