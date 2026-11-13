#!/usr/bin/env python3
"""Audit actual ImmuNexUT E-GEAD-420 public SDRF cohort structure and MD5.

No participant IDs, biosample names, age-by-disease profiles, or individual
metadata exported to GitHub. Aggregate age, diagnosis, sex, study phase only.
GEO SDRF records describe WGS blood sample inputs, NOT guaranteed per-celltype
eQTL effective n; label the difference and leave colocalization HOLD.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import statistics
import sys

from audit_EGEAD420_nominal_zip64 import (
    FILELIST,IDF,BASE,source_filelist,get_small,
)

SDRF=BASE+"E-GEAD-420.sdrf.txt"
EXPECTED=416
EXPECTED_HEALTHY=79
FIELD_KEYS=("Source Name","Characteristics[disease]","Characteristics[sex]",
            "Characteristics[age]","Characteristics[Phase]",
            "Derived Array Data Matrix File",
            "Comment[Derived Array Data Matrix File md5]")


def cohort_summary(sdrf_data):
    rows=list(csv.DictReader(io.StringIO(sdrf_data.decode("utf-8-sig","strict")),delimiter="\t"))
    if len(rows)!=EXPECTED:
        raise ValueError(f"Official E-GEAD-420 SDRF WGS subject rows not {EXPECTED}: {len(rows)}")
    if not rows or any(k not in rows[0] for k in FIELD_KEYS):
        raise ValueError("Missing official required SDRF cohort fields")
    n_missing=Counter()
    sex=Counter()
    phase=Counter()
    disease=Counter()
    age=[]
    ids=set()
    invalid_donor=0
    for row in rows:
        raw_id=row["Source Name"]
        if not raw_id or raw_id in ids:
            invalid_donor+=1
        ids.add(raw_id)
        if row["Derived Array Data Matrix File"]!="nominal_eQTL_2.tar.gz":
            raise ValueError("Source cohort references a different QTL matrix file")
        if row["Comment[Derived Array Data Matrix File md5]"]!="41544fb671177aeee7fd2d0634b90ae2":
            raise ValueError("Source participant metadata nominal QTL MD5 reference mismatch")
        for field in ("disease","sex","age","Phase"):
            key=f"Characteristics[{field}]"
            if not row[key] or row[key].strip().lower() in ("na","n/a","unknown"):
                n_missing[field]+=1
        sex[row["Characteristics[sex]"].strip().lower()]+=1
        phase[row["Characteristics[Phase]"].strip()]+=1
        disease[row["Characteristics[disease]"].strip()]+=1
        try:
            a=float(row["Characteristics[age]"])
        except ValueError:
            continue
        if not (0<a<110) or not math.isfinite(a):
            raise ValueError("SDRF age outside plausible range")
        age.append(a)
    if invalid_donor:
        raise ValueError("Original SDRF source names not unique")
    n_healthy=sum(n for name,n in disease.items() if name.lower().strip() in
                  ("healthy control","healthy","normal","control"))
    if n_healthy!=EXPECTED_HEALTHY:
        raise ValueError(f"Original healthy control status differs from official n={EXPECTED_HEALTHY}: {n_healthy}")
    if sum(disease.values())!=EXPECTED or sum(sex.values())!=EXPECTED or sum(phase.values())!=EXPECTED:
        raise ValueError("Cohort group values do not sum to total")
    return {
        "original_public_SDRF_wgs_sample_records":len(rows),
        "original_public_SDRF_distinct_source_names":len(ids),
        "healthy_control_WGS_records":n_healthy,
        "immune_disease_WGS_records":EXPECTED-n_healthy,
        "cohort_disease_category_counts":dict(disease),
        "overall_reported_sex_counts":dict(sex),
        "study_phase_counts":dict(phase),
        "fields_with_explicit_missing_records":dict(n_missing),
        "nonmissing_age_records":len(age),
        "age_median_years_aggregate":round(statistics.median(age),2) if age else None,
        "age_range_years_aggregate":[min(age),max(age)] if age else None,
        "all_eQTL_summaries_reference_same_original_TAR":True,
        "per_celltype_effective_eQTL_N_observed":False,
        "cell_type_specific_eQTL_genotypes_available":False,
        "no_participant_identifiers_exported":True,
        "source_study_is_stroke_cohort":False,
    }


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    report={
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "dataset":"E-GEAD-420 ImmuNexUT whole blood WGS cohort source metadata",
        "status":"SOURCE_OR_QC_HOLD",
        "nominal_all_cis_rows_examined":False,
        "total_sample_size_416_valid_per_each_cell_type":False,
        "QTL_coefficient_sampleN_verified":False,
        "causal_stroke_coloc":"NOT_RUN",
    }
    try:
        published=source_filelist(get_small(FILELIST,6000))
        content=get_small(SDRF,350000)
        expected=published["E-GEAD-420.sdrf.txt"]
        if len(content)!=expected["size_bytes"]:
            raise ValueError("SDRF original size differs from publisher")
        digest=hashlib.md5(content).hexdigest()
        if digest!=expected["publisher_MD5"]:
            raise ValueError("Original SDRF official MD5 check failed")
        report["original_SDRF_publisher_MD5_verified"]=True
        report["original_SDRF_MD5"]=digest
        report["aggregate_evidence"]=cohort_summary(content)
        report["status"]="SOURCE_SDRF_416_COHORT_METADATA_AND_PUBLISHER_MD5_PASS"
    except (ValueError,OSError,TimeoutError,UnicodeDecodeError,KeyError) as e:
        report["issue"]=f"{type(e).__name__}: {str(e)[:250]}"
    (out/"E_GEAD_420_cohort_metadata.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=["# ImmuNexUT E-GEAD-420 original public 416-row cohort metadata audit","",
           "**"+report["status"]+"**","",
           "These are WGS blood study participants, NOT proven QTL N in every 28 immune cell types.",""]
    if x:=report.get("aggregate_evidence"):
        lines.extend([
            f"- Distinct public sample records: {x['original_public_SDRF_distinct_source_names']}",
            f"- Healthy controls: {x['healthy_control_WGS_records']}",
            f"- Mixed immune-disease records: {x['immune_disease_WGS_records']}",
            f"- Reported age data (aggregate only): {x['nonmissing_age_records']}/{EXPECTED}",
            "",
            "| Disease group (source label) | WGS source records |",
            "|---|---:|",
        ])
        lines.extend(f"| {name} | {n} |" for name,n in sorted(x["cohort_disease_category_counts"].items()))
    if x:=report.get("issue"):
        lines.append("Source HOLD: "+x)
    lines.extend(["",
        "No donor identifiers or patient-level values exported.",
        "QTLtools clinical diagnosis/sex/phase adjustments are publisher-described; actual per-celltype N and variants/SE remain HOLD.",
    ])
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if path:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(path,"a",encoding="utf-8") as f:f.write(md)
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",default="artifacts/is_EGEAD420_sdrf")
    args=parser.parse_args()
    run(Path(args.out_dir))
