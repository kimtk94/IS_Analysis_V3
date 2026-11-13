#!/usr/bin/env python3
"""Summarize GLGC chromosome SNP-QC evidence, never infer full-GWAS PASS."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path

TRAITS=("LDL_C","HDL_C","TG")
EXCLUSION_FIELDS=("non_snp","malformed_columns","wrong_chromosome",
                  "invalid_coordinate","invalid_alleles","invalid_effect_or_se",
                  "p_underflow_zero","invalid_p","invalid_allele_frequency")
SCOPE="ONE_FULL_CHROMOSOME_ONLY_NOT_WHOLE_GWAS"


def summarize(folder,chromosome):
    if chromosome not in {str(i) for i in range(1,23)}:
        raise ValueError("invalid chromosome")
    entries={}
    errors=[]
    for path in sorted(Path(folder).glob("*.json")):
        try:
            row=json.loads(path.read_text(encoding="utf-8"))
        except (OSError,ValueError) as exc:
            errors.append(str(path)+": invalid JSON "+str(exc)[:200])
            continue
        if row.get("scope")!=SCOPE: continue
        if str(row.get("chromosome"))!=chromosome: continue
        trait=row.get("trait")
        if trait not in TRAITS:
            errors.append(str(path)+": unexpected trait "+str(trait))
            continue
        if trait in entries:
            errors.append("duplicate artifact for "+trait+" chr"+chromosome)
            continue
        counts=row.get("counts",{})
        variant=row.get("variant_qc",{})
        excluded=sum(counts.get(key,0) for key in EXCLUSION_FIELDS)
        original=counts.get("original_rows",0)
        retained=counts.get("retained_autosomal_snp",0)
        checks={
            "observed_header":row.get("header_evidence",{}).get("status")=="PASS",
            "source_region_only":row.get("scope")==SCOPE,
            "shard_complete":not row.get("row_cap_reached",True),
            "row_accounting":original==retained+excluded and original>0,
            "variant_qc_rows":variant.get("rows")==retained and retained>0,
            "no_duplicate_variant_keys":variant.get("duplicate_excess_rows")==0,
            "no_p_diagnostic_mismatches":variant.get("p_consistency_mismatches")==0,
            "source_identity":row.get("source")=="GLGC_Graham2021" and row.get("ancestry")=="EUR",
            "technical_gate_passed":row.get("status")=="PASS" and variant.get("status")=="PASS",
        }
        entries[trait]={
            "source_commit":row.get("source_commit",""),
            "technical_status":"PASS" if all(checks.values()) else "HOLD",
            "checks":checks,
            "original_rows":original,
            "valid_snp_rows":retained,
            "excluded_non_snp":counts.get("non_snp",0),
            "excluded_other":excluded-counts.get("non_snp",0),
            "palindromic_ambiguous_af":counts.get("palindromic_ambiguous_af_0_42_to_0_58",0),
            "duplicate_keys":variant.get("duplicate_excess_rows"),
            "p_consistency_sampled":variant.get("p_consistency_sampled"),
            "p_consistency_mismatches":variant.get("p_consistency_mismatches"),
        }
    missing=[t for t in TRAITS if t not in entries]
    commits={v["source_commit"] for v in entries.values()}
    if len(commits)>1 or "" in commits and len(entries)>0:
        errors.append("missing or inconsistent source commits")
    status="PASS" if not missing and not errors and all(v["technical_status"]=="PASS" for v in entries.values()) else "HOLD"
    summary={
        "status":status,
        "scientific_status":"HOLD",
        "scope":"ONE_FULL_CHROMOSOME_ONLY_NOT_WHOLE_GWAS",
        "chromosome":chromosome,
        "expected_traits":list(TRAITS),
        "missing_traits":missing,
        "errors":errors,
        "traits":entries,
        "snp_rows_total":sum(e["valid_snp_rows"] for e in entries.values()),
        "source_gwas_whole_genome_qc":"NOT_RUN",
        "allele_harmonization":"NOT_RUN",
        "downstream_ldsc_or_mr":"NOT_RUN",
        "note":"PASS covers technical SNP-only QC on one selected chromosome across 3 traits, not all 22 autosomes or source-file checksums."
    }
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input-dir",default="artifacts/ckm_stroke")
    parser.add_argument("--chromosome",default="22")
    parser.add_argument("--output",default="artifacts/ckm_stroke_lipid_chr22_summary.json")
    args=parser.parse_args()
    summary=summarize(args.input_dir,args.chromosome)
    summary["timestamp_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    summary["source_commit"]=os.getenv("GITHUB_SHA","")
    target=Path(args.output)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False))
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as output:
            output.write("\n## GLGC chr"+args.chromosome+" summary\n")
            output.write("\n| Trait | Original rows | SNP rows | Technical QC |\n|---|---:|---:|---|\n")
            for trait in TRAITS:
                item=summary["traits"].get(trait,{})
                output.write("| "+trait+" | "+str(item.get("original_rows","—"))+" | "+
                             str(item.get("valid_snp_rows","—"))+" | "+
                             str(item.get("technical_status","MISSING"))+" |\n")
            output.write("\n**Overall technical: "+summary["status"]+"; Scientific: HOLD.**\n")
    return int(summary["status"]!="PASS")


if __name__=="__main__":
    raise SystemExit(main())
