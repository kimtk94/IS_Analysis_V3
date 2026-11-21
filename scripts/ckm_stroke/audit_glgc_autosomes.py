#!/usr/bin/env python3
"""Fail-closed technical coverage matrix for GLGC EUR autosome SNP-only QC.

Reads only source-linked per-chromosome JSON evidence created by the
per-trait and per-chromosome CI. Can show a truthful partial snapshot.
Never interprets this as whole-source hash, LDSC or MR evidence.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path

AUTOSOMES=tuple(str(i) for i in range(1,23))
TRAITS=("LDL_C","HDL_C","TG")
SCOPE="ONE_FULL_CHROMOSOME_ONLY_NOT_WHOLE_GWAS"
MUST_CHECK={
    "observed_header","source_region_only","shard_complete",
    "row_accounting","variant_qc_rows","no_duplicate_variant_keys",
    "no_p_diagnostic_mismatches","source_identity","technical_gate_passed"
}


def assess(directory):
    directory=Path(directory)
    evidence={}
    errors=[]
    ignored=[]
    for path in sorted(directory.glob("*.json")):
        try:
            row=json.loads(path.read_text(encoding="utf-8"))
        except (OSError,ValueError) as exc:
            errors.append(path.name+": invalid JSON "+str(exc)[:200])
            continue
        if row.get("scope")!=SCOPE or row.get("chromosome") is None:
            ignored.append(path.name)
            continue
        chrom=str(row["chromosome"])
        if chrom not in AUTOSOMES:
            errors.append(path.name+": invalid chromosome "+chrom)
            continue
        if chrom in evidence:
            errors.append("duplicate chromosome aggregate: chr"+chrom)
            continue
        evidence[chrom]=row

    coverage=[]
    summary={"source_rows":0,"retained_snp_rows":0,"multibase_non_snp":0,
             "symbolic_or_ambiguous_non_snp":0,"invalid_other":0,
             "duplicate_snp_keys":0,"palindromic_ambiguous_af":0,
             "underflow_p_zero_rows":0,
             "underflow_unreported_legacy_trait_chromosome_results":0}
    code_commits=set()
    for chrom in AUTOSOMES:
        raw=evidence.get(chrom)
        if raw is None:
            coverage.append({"chromosome":chrom,"status":"NOT_RUN",
                             "traits":{t:"NOT_RUN" for t in TRAITS}})
            continue
        seen=set(raw.get("traits",{}))
        values={}
        if raw.get("scope")!=SCOPE or raw.get("chromosome")!=chrom:
            errors.append("chr"+chrom+": invalid scope")
        if raw.get("missing_traits") or raw.get("errors"):
            errors.append("chr"+chrom+": missing/errors in chromosome evidence")
        if seen!=set(TRAITS):
            errors.append("chr"+chrom+": exactly three traits required")
        if raw.get("status")!="PASS" or raw.get("scientific_status")!="HOLD":
            errors.append("chr"+chrom+": technical gate not PASS")
        source_commit=raw.get("source_commit")
        if not isinstance(source_commit,str) or len(source_commit)!=40 or any(c not in "0123456789abcdef" for c in source_commit):
            errors.append("chr"+chrom+": source commit not a 40-character lowercase SHA")
        else:
            code_commits.add(source_commit)
        valid_total=0
        for t in TRAITS:
            record=raw.get("traits",{}).get(t)
            if not isinstance(record,dict):
                values[t]="NOT_RUN"
                continue
            checks=record.get("checks",{})
            good=(record.get("technical_status")=="PASS"
                  and isinstance(checks,dict)
                  and MUST_CHECK.issubset(checks)
                  and all(checks.get(k) is True for k in MUST_CHECK)
                  and record.get("duplicate_keys")==0
                  and isinstance(record.get("original_rows"),int) and record["original_rows"]>0
                  and isinstance(record.get("valid_snp_rows"),int) and record["valid_snp_rows"]>0
                  and record.get("source_commit")==source_commit)
            values[t]="PASS" if good else "HOLD"
            if good:
                valid_total+=record["valid_snp_rows"]
                summary["source_rows"]+=record["original_rows"]
                summary["retained_snp_rows"]+=record["valid_snp_rows"]
                summary["multibase_non_snp"]+=record.get("excluded_non_snp",0)
                summary["symbolic_or_ambiguous_non_snp"]+=record.get("excluded_symbolic_or_ambiguous_non_snp",0)
                summary["invalid_other"]+=record.get("excluded_other",0)
                summary["duplicate_snp_keys"]+=record.get("duplicate_keys",0)
                summary["palindromic_ambiguous_af"]+=record.get("palindromic_ambiguous_af",0)
                if "p_underflow_zero" in record:
                    summary["underflow_p_zero_rows"]+=record["p_underflow_zero"]
                else:
                    # Earlier approved aggregates lack this column; never
                    # mistake unreported underflow counts for actual zero.
                    summary["underflow_unreported_legacy_trait_chromosome_results"]+=1
        if valid_total!=raw.get("snp_rows_total"):
            errors.append("chr"+chrom+": SNP row total mismatches verified traits")
        status="PASS" if all(v=="PASS" for v in values.values()) and len(values)==3 and raw.get("status")=="PASS" else "HOLD"
        coverage.append({"chromosome":chrom,"status":status,"traits":values})
    observed=sum(x["status"]=="PASS" for x in coverage)
    no_evidence=sum(x["status"]=="NOT_RUN" for x in coverage)
    holds=sum(x["status"]=="HOLD" for x in coverage)
    if errors or holds:
        technical="HOLD"
    elif no_evidence:
        technical="PARTIAL"
    else:
        technical="PASS"
    if technical=="PASS" and observed!=22:
        raise AssertionError("coverage gate: full autosome count inconsistent")
    return {
        "status":technical,"scientific_status":"HOLD",
        "scope":"GLGC_EUR_AUTOSOME_SNP_ONLY_TECHNICAL_COVERAGE",
        "observed_chromosome_pass":observed,
        "missing_chromosome_count":no_evidence,
        "hold_chromosome_count":holds,
        "expected_chromosomes":22,
        "expected_trait_chromosome_jobs":66,
        "trait_chromosome_pass":sum(v=="PASS" for x in coverage for v in x["traits"].values()),
        "source_commits":sorted(code_commits),
        "mixed_code_commits_review_required":len(code_commits)>1,
        "coverage":coverage,
        "totals_verified_pass_only":summary,
        "ignored_files":ignored,
        "errors":errors,
        "whole_source_checksum_verified":False,
        "ldsc_ready":"NOT_ASSESSED",
        "mr_ready":"NOT_ASSESSED",
        "source_snp_export_scientific_ready":False,
        "underflow_recovery_required":(
            summary["underflow_p_zero_rows"]>0
            or summary["underflow_unreported_legacy_trait_chromosome_results"]>0),
        "note":"This checks SNP-only per-chromosome metadata. A 22/22 PASS is not original full-file source integrity or statistical model readiness."
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    report=assess(args.input_dir)
    report["timestamp_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    report["summary_generated_from_commit"]=os.getenv("GITHUB_SHA","")
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k!="coverage"},ensure_ascii=False))
    return int(report["status"]=="HOLD")


if __name__=="__main__":
    raise SystemExit(main())
