#!/usr/bin/env python3
"""Real publisher-MD5 EAS GIGASTROKE outcome EAF coverage and MAF QC.

Study allele-frequency distribution is outcome GWAS (GCST90104545), NOT
ImmuNexUT genotype/celltype eQTL MAF and NOT per-cell eQTL sample N. Does
not perform colocalization, candidate-gene selection or clinical inference.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import gzip
import json
import math
import os
from pathlib import Path
import tempfile

from audit_gwas_1000g_bim_overlap import GWAS_FTP,SOURCES,canonical_chrom
from audit_gwas_1000g_bim_all22 import resumable_source_download

ANCESTRY="EAS"
EXPECTED_ROWS=6725216
COLUMNS=("chromosome","base_pair_location","effect_allele",
         "other_allele","effect_allele_frequency","p_value")
PALINDROMIC={frozenset(("A","T")),frozenset(("C","G"))}
BINS=(("lt_0_001",0.001),("0_001_to_0_01",0.01),
      ("0_01_to_0_05",0.05),("0_05_to_0_1",0.1),
      ("0_1_to_0_4",0.4),("0_4_to_0_5",0.50000001))


def classify_freq(raw):
    if raw.strip().lower() in ("na","nan","null","none",".","","-"):
        return None
    try:
        value=float(raw)
    except ValueError as e:
        raise ValueError("Malformed original EAS GWAS allele frequency") from e
    if not math.isfinite(value) or not 0<=value<=1:
        raise ValueError("Out-of-bounds original EAS GWAS allele frequency")
    return value


def maf_bin(maf):
    if not 0<=maf<=0.5:
        raise ValueError("Invalid minor allele frequency")
    for name,limit in BINS:
        if maf<limit:
            return name
    raise ValueError("No MAF bin")


def audit_original_eaf(gzip_file):
    counts=Counter()
    by_chr=defaultdict(Counter)
    with gzip.open(gzip_file,"rt",encoding="utf-8",newline="") as source:
        header=source.readline().rstrip("\r\n").split("\t")
        if any(v not in header for v in COLUMNS):
            raise ValueError("Original GIGASTROKE EAS source column missing")
        ix={v:header.index(v) for v in COLUMNS}
        for row in source:
            if not row.strip():continue
            fields=row.rstrip("\r\n").split("\t")
            if len(fields)!=len(header):
                raise ValueError("Original EAS GWAS row column count mismatch")
            chrom=canonical_chrom(fields[ix["chromosome"]])
            if chrom is None:
                raise ValueError("Source GWAS row not a valid autosome")
            counts["original_autosomal_rows"]+=1
            by_chr[str(chrom)]["gwas_rows"]+=1
            eaf=classify_freq(fields[ix["effect_allele_frequency"]])
            if eaf is None:
                counts["eaf_missing"]+=1
                by_chr[str(chrom)]["eaf_missing"]+=1
                continue
            counts["eaf_valid"]+=1
            maf=min(eaf,1-eaf)
            name=maf_bin(maf)
            counts["maf_"+name]+=1
            by_chr[str(chrom)]["maf_"+name]+=1
            alleles=frozenset((fields[ix["effect_allele"]].upper(),
                              fields[ix["other_allele"]].upper()))
            if alleles in PALINDROMIC:
                counts["palindromic_eaf_valid"]+=1
                if maf>=0.4:
                    counts["palindromic_maf_ge_0_4_frequency_ambiguous"]+=1
            try:
                p=float(fields[ix["p_value"]])
            except ValueError as exc:
                raise ValueError("Non-numeric original EAS GWAS p-value") from exc
            if not math.isfinite(p) or not 0<=p<=1:
                raise ValueError("Invalid original EAS GWAS p-value")
            if p<5e-8:
                counts["gwas_significant_eaf_valid_rows"]+=1
    if counts["original_autosomal_rows"]!=EXPECTED_ROWS or len(by_chr)!=22:
        raise ValueError("GWAS original 22-autosome source rows differ from publisher-verified baseline")
    if counts["eaf_valid"]+counts["eaf_missing"]!=EXPECTED_ROWS:
        raise ValueError("EAS EAF counting does not conserve all source rows")
    if sum(counts["maf_"+b] for b,_ in BINS)!=counts["eaf_valid"]:
        raise ValueError("EAS MAF frequency-bin denominator mismatch")
    summary={
        "rows":dict(counts),
        "per_chromosome":{x:dict(v) for x,v in sorted(by_chr.items(),key=lambda v:int(v[0]))},
        "frequency_coverage_pct":round(100*counts["eaf_valid"]/EXPECTED_ROWS,4),
        "eaf_missing_pct":round(100*counts["eaf_missing"]/EXPECTED_ROWS,4),
        "n_ancestry_specific_sample_frequency_source":"OUTCOME_GWAS_EAS_NOT_EQTL",
        "maf_is_from_immuNexUT_donors":False,
        "maf_reference_for_coloc":"NOT_YET_VALIDATED",
        "GWAS_to_eQTL_effect_direction_harmonized":False,
        "palindromic_frequency_resolved":False,
        "all_22_autosomes_checked":True,
    }
    return summary


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_EAS_GWAS_eaf_QC")
    args=p.parse_args()
    out=Path(args.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    report={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "outcome_source":"GIGASTROKE GCST90104545 EAS ischemic stroke GRCh37",
        "status":"SOURCE_OR_QC_HOLD",
        "dataset_kind":"OUTCOME_GWAS_NOT_EXPOSURE_EQTL",
        "coloc_GWAS_to_eQTL":"NOT_RUN",
        "per_celltype_N":"NOT_VERIFIED",
    }
    try:
        s=SOURCES[ANCESTRY]
        accession=s["gwas_accession"]
        url=GWAS_FTP+accession+"/"+accession+"_buildGRCh37.tsv.gz"
        with tempfile.TemporaryDirectory(prefix="is_outcome_eaf_") as tmp:
            path=Path(tmp)/"original_EAS_GWAS.gz"
            report["original_publisher_md5_proof"]=resumable_source_download(
                url,path,s["gwas_publisher_md5"],s["gwas_publisher_size"],attempts=4)
            report["all_22_autosome_allele_frequency_evidence"]=audit_original_eaf(path)
        report["status"]="ORIGINAL_GWAS_EAS_EAF_22_AUTOSOME_QC_PASS"
    except (OSError,ValueError,TimeoutError,KeyError) as exc:
        report["issue"]=f"{type(exc).__name__}: {str(exc)[:220]}"
    (out/"EAS_GWAS_original_eaf_counts.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=["# Original EAS GIGASTROKE IS GWAS EAF / MAF technical QC",
           "", "**"+report["status"]+"**","",
           "Outcome GWAS effect allele frequencies are NOT eQTL subject frequencies; no coloc or joint biology inference."]
    if agg:=report.get("all_22_autosome_allele_frequency_evidence"):
        counts=agg["rows"]
        lines.extend([
            f"- Whole EAS GWAS autosomal rows: {counts['original_autosomal_rows']:,}",
            f"- Valid outcome-GWAS EAF: {counts['eaf_valid']:,} ({agg['frequency_coverage_pct']}%)",
            f"- Missing outcome-GWAS EAF: {counts['eaf_missing']:,} ({agg['eaf_missing_pct']}%)",
            f"- Palindromic A/T C/G with EAF: {counts['palindromic_eaf_valid']:,}",
            f"- Palindromic with near-0.5 EAF (MAF>=0.4): {counts['palindromic_maf_ge_0_4_frequency_ambiguous']:,}",
            "",
            "| GWAS outcome MAF band | Rows |","|---|---:|",
        ])
        for band,_ in BINS:
            lines.append(f"| {band} | {counts['maf_'+band]:,} |")
    if err:=report.get("issue"):
        lines.append("- QC HOLD: "+err)
    lines.extend(["","**NEVER use these outcome EAFs as if measured in the ImmuNexUT CD16+ monocyte eQTL subset.**",
                  "Actual cell-type N, EAS reference frequency, liftover, REF/ALT and candidate stroke loci still require independent checks."])
    summary="\n".join(lines)+"\n"
    (out/"summary.md").write_text(summary,encoding="utf-8")
    print(summary)
    if path:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(path,"a",encoding="utf-8") as f:
            f.write(summary)
    if report["status"]!="ORIGINAL_GWAS_EAS_EAF_22_AUTOSOME_QC_PASS":
        raise SystemExit(1)


if __name__=="__main__":
    main()
