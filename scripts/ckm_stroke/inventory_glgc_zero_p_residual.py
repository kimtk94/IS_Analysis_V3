#!/usr/bin/env python3
"""Fail-closed inventory of residual GLGC other/zero-P candidates.

Reads the 22 original chromosome aggregate JSONs. Does not assert that
legacy 'other' rows are P=0 or scientifically eligible.
"""
import argparse
import csv
import json
from pathlib import Path

TRAITS=("LDL_C","HDL_C","TG")
VERIFIED={("19","LDL_C"):332,("8","HDL_C"):286,
          ("8","TG"):341,("2","LDL_C"):378,("2","TG"):199}
TOTAL_CANDIDATES=2632
VERIFIED_TOTAL=1536

def inventory(folder):
    data={}
    for path in sorted(Path(folder).rglob("*.json")):
        d=json.loads(path.read_text(encoding="utf-8"))
        if d.get("scope")!="ONE_FULL_CHROMOSOME_ONLY_NOT_WHOLE_GWAS":
            continue
        chrom=str(d.get("chromosome",""))
        if chrom in data or chrom not in {str(i) for i in range(1,23)}:
            raise ValueError("duplicate or invalid chromosome: "+chrom)
        if d.get("status")!="PASS" or d.get("scientific_status")!="HOLD":
            raise ValueError("unverified chromosome aggregate: "+chrom)
        data[chrom]=d
    if len(data)!=22:
        raise ValueError(f"expected 22 original chromosome aggregates, got {len(data)}")
    rows=[]
    for chrom in map(str,range(1,23)):
        d=data[chrom]
        for trait in TRAITS:
            item=d.get("traits",{}).get(trait)
            if not isinstance(item,dict) or item.get("technical_status")!="PASS":
                raise ValueError("missing original technical PASS "+chrom+" "+trait)
            other=item.get("excluded_other")
            if not isinstance(other,int) or other<0:
                raise ValueError("invalid other candidate count")
            reported=item.get("p_underflow_zero")
            if reported is not None and (not isinstance(reported,int) or reported<0 or reported>other):
                raise ValueError("invalid original underflow accounting")
            prior=VERIFIED.get((chrom,trait),0)
            if prior>other or (reported is not None and prior>reported):
                raise ValueError("verified count exceeds original candidates")
            residual=other-prior
            rows.append({
                "chromosome":chrom,"trait":trait,"original_other_candidates":other,
                "reported_source_zero_p":reported if reported is not None else "",
                "underflow_schema":"REPORTED" if reported is not None else "LEGACY_UNREPORTED",
                "previously_validated_zero_p":prior,
                "unverified_other_candidates":residual,
                "known_zero_p_to_validate":max(0,reported-prior) if reported is not None else "",
                "unclassified_legacy_candidates":residual if reported is None else "",
                "source_commit":d.get("source_commit",""),
                "scientific_status":"HOLD"
            })
    total=sum(x["original_other_candidates"] for x in rows)
    passed=sum(x["previously_validated_zero_p"] for x in rows)
    remaining=sum(x["unverified_other_candidates"] for x in rows)
    if total!=TOTAL_CANDIDATES or passed!=VERIFIED_TOTAL or remaining!=1096:
        raise ValueError(f"baseline mismatch: candidates={total}, verified={passed}, residual={remaining}")
    return {"status":"PASS","scientific_status":"HOLD",
            "scope":"GLGC_ORIGINAL_OTHER_CANDIDATE_INVENTORY_NOT_CAUSAL",
            "candidate_total":total,"prior_verified_zero_p":passed,
            "remaining_unverified_candidates":remaining,
            "known_zero_p_remaining":sum(x["known_zero_p_to_validate"] for x in rows if isinstance(x["known_zero_p_to_validate"],int)),
            "legacy_unclassified_remaining":sum(x["unclassified_legacy_candidates"] for x in rows if isinstance(x["unclassified_legacy_candidates"],int)),
            "rows":rows,
            "warning":"Legacy other cannot be called P=0 without re-reading original source. Counts do not establish MR or LDSC readiness."}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input-dir",required=True)
    p.add_argument("--output-prefix",required=True)
    a=p.parse_args()
    report=inventory(a.input_dir)
    dest=Path(a.output_prefix)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    with dest.with_suffix(".tsv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(report["rows"][0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(report["rows"])
    print(json.dumps({k:v for k,v in report.items() if k!="rows"}))
if __name__=="__main__":
    main()
