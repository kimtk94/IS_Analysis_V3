#!/usr/bin/env python3
"""48-patient Ng 2022 source clinical covariate balance, privacy-safe.

Requires exactly matched 50 public RDS biopsy IDs against source-SHA256-pinned
Ng Supplementary Data 1. Outputs no patient rows or IDs, no clinical raw values,
no individual age, no small (<5) category counts.
"""
import argparse
import collections
import datetime as dt
import json
import math
import pathlib
import re
import statistics

from probe_ng_gao_public import _normalized_label
from probe_ng_suter_clinical_sample_join import official_rows, categories, SAMPLE_N

HEADER_ALLOWLIST={
 "age":re.compile(r"^(?:age|ageatclassification|ageatdiagnosis|ageatbiopsy)(?:years)?$"),
 "sex":re.compile(r"^(?:sex|gender)$"),
 "bclc":re.compile(r"^bclc(?:stage)?$"),
 "edmondson":re.compile(r"^edmondson(?:grade)?$"),
 "child_pugh":re.compile(r"^(?:child|childpugh(?:score|class)?)$"),
 "meld":re.compile(r"^meld(?:score)?$"),
 "cirrhosis":re.compile(r"^cirrhosis$"),
 "treatment":re.compile(r"^(?:treatment|prior(?:systemic)?treatment|therapy|therapybeforebiopsy)$"),
 "batch":re.compile(r"^(?:batch|sequencingbatch|proteomicsbatch|swathbatch)$")
}
GROUPS=("ALD","HCV")
SAFE_MIN=5

def canonical_columns(columns):
    found={}
    for name,rule in HEADER_ALLOWLIST.items():
        hits=[index for k,index in columns.items() if rule.fullmatch(k)]
        if len(hits)==1:found[name]=hits[0]
    return found

def parse_age(s):
    s=str(s).strip()
    try:
        v=float(s)
        return v if math.isfinite(v) and 0<v<=110 else None
    except (ValueError,TypeError):return None

def qtile(nums,p):
    nums=sorted(nums)
    if not nums:return None
    pos=(len(nums)-1)*p
    lo=int(math.floor(pos)); hi=int(math.ceil(pos))
    return round(nums[lo]+(nums[hi]-nums[lo])*(pos-lo),2)

def numeric_abs_smd(a,b):
    if len(a)<2 or len(b)<2:return None
    va=statistics.variance(a);vb=statistics.variance(b)
    scale=(va+vb)/2
    return round(abs((statistics.mean(a)-statistics.mean(b))/(scale**.5)),3) if scale>0 else None

def binary_abs_smd(a,b):
    if not a or not b:return None
    pa=sum(a)/len(a);pb=sum(b)/len(b)
    scale=(pa*(1-pa)+pb*(1-pb))/2
    return round(abs((pa-pb)/(scale**.5)),3) if scale>0 else (0.0 if pa==pb else None)

def safe_counter(items):
    counts=collections.Counter(items)
    return {k:(n if n>=SAFE_MIN else None) for k,n in sorted(counts.items())}

def validate_grouping(private_labels,rows,columns):
    if len(private_labels)!=SAMPLE_N or len(set(private_labels))!=SAMPLE_N:
        raise ValueError("source RDS labels must be exactly 50 distinct")
    primary={}
    for row in rows:
        k=_normalized_label(row.get(columns["tumorbiopsyid"],""))
        if k:
            if k in primary:raise ValueError("source repeated biopsy label")
            primary[k]=row
    if any(k not in primary for k in private_labels):
        raise ValueError("RDS biopsy not found in source official table")
    by_patient={}
    for biopsy in private_labels:
        row=primary[biopsy]
        patient=_normalized_label(row.get(columns["patientid"],""))
        if not patient:raise ValueError("source patient label missing")
        by_patient.setdefault(patient,[]).append(row)
    if len(by_patient)!=48:raise ValueError("unique source patient count mismatch")
    return by_patient

def audit(private_labels,rows,columns):
    by_patient=validate_grouping(private_labels,rows,columns)
    clinical=canonical_columns(columns)
    patient_groups={};patient_covariates={};conflicts=collections.Counter()
    if "underlyingliverdisease" not in columns:
        raise ValueError("source etiology header absent")
    for patient,patient_rows in by_patient.items():
        labels={categories(r.get(columns["underlyingliverdisease"],"")) for r in patient_rows}
        # Ensure multi-biopsy patient has one agreed exclusive source category.
        if len(labels)==1:
            label=next(iter(labels))
            group=label[0] if len(label)==1 and label[0] in GROUPS else "MIXED_OTHER"
        else:
            group="CONFLICT"
        patient_groups[patient]=group
        record={}
        for field,index in clinical.items():
            values={str(r.get(index,"")).strip() for r in patient_rows if str(r.get(index,"")).strip()}
            if len(values)>1:
                conflicts[field]+=1
                record[field]=None
            elif len(values)==1:
                record[field]=next(iter(values))
            else:
                record[field]=None
        patient_covariates[patient]=record
    groups={g:[p for p in by_patient if patient_groups[p]==g] for g in GROUPS}
    if len(groups["ALD"])!=12 or len(groups["HCV"])!=11:
        raise ValueError("source strata 12/11 expected from pinned audit")
    summary={}
    for g,patients in groups.items():
        result={"unique_patients_n":len(patients)}
        if "age" in clinical:
            values=[parse_age(patient_covariates[p].get("age")) for p in patients]
            nums=[v for v in values if v is not None]
            result["age"]={"nonmissing_n":len(nums),
              "median":qtile(nums,.5) if len(nums)>=8 else None,
              "q1":qtile(nums,.25) if len(nums)>=8 else None,
              "q3":qtile(nums,.75) if len(nums)>=8 else None}
        if "sex" in clinical:
            values=[str(patient_covariates[p].get("sex") or "").lower() for p in patients]
            sex=[("male" if v in ("m","male") else
                  "female" if v in ("f","female") else "unclassified") for v in values]
            result["sex"]={"counts_suppress_lt5":safe_counter(sex)}
        if "bclc" in clinical:
            vals=[str(patient_covariates[p].get("bclc") or "").upper() for p in patients]
            stage=[("early_0_A" if x in ("0","A") else
                    "advanced_B_C_D" if x in ("B","C","D") else
                    "unclassified") for x in vals]
            result["bclc"]={"counts_suppress_lt5":safe_counter(stage)}
        if "child_pugh" in clinical:
            vals=[str(patient_covariates[p].get("child_pugh") or "").upper().strip()
                  for p in patients]
            classes=[("A" if v=="A" else "B_C" if v in ("B","C") else
                      "unclassified") for v in vals]
            result["child_pugh"]={"present_n":sum(v!="unclassified" for v in classes),
                                   "counts_suppress_lt5":safe_counter(classes)}
        if "cirrhosis" in clinical:
            vals=[str(patient_covariates[p].get("cirrhosis") or "").lower().strip()
                  for p in patients]
            status=[("yes" if v in ("yes","y","1","true") else
                     "no" if v in ("no","n","0","false") else
                     "unclassified") for v in vals]
            result["cirrhosis"]={"present_n":sum(v!="unclassified" for v in status),
                                  "counts_suppress_lt5":safe_counter(status)}
        if "meld" in clinical:
            nums=[]
            for p in patients:
                try:
                    v=float(str(patient_covariates[p].get("meld") or ""))
                    if math.isfinite(v) and 0<=v<=60: nums.append(v)
                except (TypeError,ValueError):pass
            result["meld"]={"numeric_present_n":len(nums),
                            "median":qtile(nums,.5) if len(nums)>=8 else None,
                            "q1":qtile(nums,.25) if len(nums)>=8 else None,
                            "q3":qtile(nums,.75) if len(nums)>=8 else None}
        for field in ("edmondson","treatment","batch"):
            if field in clinical:
                result[field]={"present_n":sum(
                    patient_covariates[p].get(field) is not None for p in patients)}
        summary[g]=result
    a=summary["ALD"].get("age",{})
    h=summary["HCV"].get("age",{})
    unadjusted_median_difference=(round(a["median"]-h["median"],2)
        if a.get("median") is not None and h.get("median") is not None else None)
    age_values={g:[parse_age(patient_covariates[p].get("age"))
                  for p in groups[g]] for g in GROUPS}
    age_obs={g:[v for v in age_values[g] if v is not None] for g in GROUPS}
    age_smd=numeric_abs_smd(age_obs["ALD"],age_obs["HCV"])
    sex_values={g:[str(patient_covariates[p].get("sex") or "").lower()
                   for p in groups[g]] for g in GROUPS}
    sex_male={g:[v in ("m","male") for v in sex_values[g]
                 if v in ("m","male","f","female")] for g in GROUPS}
    sex_smd=binary_abs_smd(sex_male["ALD"],sex_male["HCV"])
    stage_values={g:[str(patient_covariates[p].get("bclc") or "").upper()
                     for p in groups[g]] for g in GROUPS}
    bclc_early={g:[v in ("0","A") for v in stage_values[g]
                  if v in ("0","A","B","C","D")] for g in GROUPS}
    bclc_smd=binary_abs_smd(bclc_early["ALD"],bclc_early["HCV"])
    child_values={g:[str(patient_covariates[p].get("child_pugh") or "").upper().strip()
                     for p in groups[g]] for g in GROUPS}
    child_binary={g:[v=="A" for v in child_values[g] if v in ("A","B","C")]
                  for g in GROUPS}
    child_smd=binary_abs_smd(child_binary["ALD"],child_binary["HCV"])
    smds={"age":age_smd,"male_sex":sex_smd,"BCLC_early_0_A":bclc_smd,
          "ChildPugh_A_vs_BC":child_smd}
    complete_separation={
        "age":len(age_obs["ALD"])>=2 and len(age_obs["HCV"])>=2
            and len(set(age_obs["ALD"]+age_obs["HCV"]))>1
            and all(len(set(age_obs[g]))==1 for g in GROUPS),
        "male_sex":len(sex_male["ALD"])>0 and len(sex_male["HCV"])>0
            and len(set(sex_male["ALD"]))==len(set(sex_male["HCV"]))==1
            and sex_male["ALD"][0]!=sex_male["HCV"][0],
        "ChildPugh_A_vs_BC":len(child_binary["ALD"])>0 and len(child_binary["HCV"])>0
            and len(set(child_binary["ALD"]))==len(set(child_binary["HCV"]))==1
            and child_binary["ALD"][0]!=child_binary["HCV"][0],
        "BCLC_early_0_A":len(bclc_early["ALD"])>0 and len(bclc_early["HCV"])>0
            and len(set(bclc_early["ALD"]))==len(set(bclc_early["HCV"]))==1
            and bclc_early["ALD"][0]!=bclc_early["HCV"][0]}
    imbalance={k:(True if complete_separation[k] else
                   v>=0.2 if v is not None else None) for k,v in smds.items()}

    return {
      "status":"ORIGINAL_SOURCE_EXCLUSIVE_ETIOLOGY_CONFOUNDER_SCREEN",
      "source_biopsies_matched_n":50,
      "source_unique_patients_n":len(by_patient),
      "exclusive_group_counts":{"ALD":12,"HCV":11},
      "available_allowlisted_clinical_headers":sorted(clinical),
      "unavailable_clinical_fields":sorted(set(HEADER_ALLOWLIST)-set(clinical)),
      "multibiopsy_conflict_fields_suppress_lt5":{
          field:(n if n>=5 else None) for field,n in conflicts.items()},
      "groups":summary,
      "ALD_minus_HCV_unadjusted_age_median_years":unadjusted_median_difference,
      "absolute_standardized_mean_difference":smds,
      "complete_group_separation_no_finite_SMD":complete_separation,
      "imbalance_attention_threshold_0_2":imbalance,
      "source_child_pugh_header_detected": "child_pugh" in clinical,
      "candidate_baseline_adjustment_ready":False,
      "claim_context":"descriptive imbalance screen, not hypothesis test or clinical adjudication",
      "independent_patient_unit_used":True,
      "clinical_prognosis_coding_reviewed":False,
      "treatment_confounded_effect_disproved":False,
      "batch_confounded_effect_disproved":False,
      "PILOT_ONLY":True,
      "full_differential_omics_GO":False,
      "individual_patient_or_sample_values_in_output":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-labels",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    output=pathlib.Path(args.out)
    output.parent.mkdir(parents=True,exist_ok=True)
    report={"utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "source":"Ng 2022 official SHA256-verified Supplementary Data 1",
            "audit":None}
    try:
        inp=pathlib.Path(args.private_labels)
        if not inp.is_file() or inp.stat().st_size>20000:raise ValueError("invalid private sample labels")
        labels=[_normalized_label(line) for line in inp.read_text(encoding="utf-8").splitlines()]
        rows,columns=official_rows()
        report["audit"]=audit(labels,rows,columns)
    except Exception as e:
        report["audit"]={"status":"UNVERIFIED","error_class":type(e).__name__,
                         "PILOT_ONLY":True,"individual_patient_or_sample_values_in_output":False}
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))
    if report["audit"]["status"]=="UNVERIFIED":raise SystemExit(1)

if __name__=="__main__":main()
