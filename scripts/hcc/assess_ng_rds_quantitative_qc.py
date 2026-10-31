#!/usr/bin/env python3
"""Aggregate quantitative Ng/Suter RDS QC: verify numeric ranges and feature identity.

Reads only sanitized source audit JSON and prints allowed counts/quantiles.
No raw biological sample identifiers or individual gene expression values.
"""
import argparse,json,math,pathlib

FILES=("HCCfull.rds","HCCfullz.rds","HCCdata.rds")
OMICS=("T","P","PP")

def check(a):
    if a.get("source_commit")!="438f52271d0a86a6701baa7a96e80c516ba575c2":
        raise ValueError("unexpected unpinned source")
    objs=a.get("objects") or {}
    result={}
    for fname in FILES:
        obj=objs.get(fname) or {}
        if obj.get("status")!="PUBLIC_MULTIMODAL_RDS_SHAPES_AUDITED":
            raise ValueError("not a numeric multimodal source")
        layers=obj.get("layers") or {}
        cleaned={}
        for omic in OMICS:
            rec=layers.get(omic) or {}
            n=rec.get("feature_count"); m=rec.get("sample_count")
            if not isinstance(n,int) or not isinstance(m,int) or n<100 or m!=50:
                raise ValueError("invalid 50-sample quantitative geometry")
            total=n*m
            finite=rec.get("finite_numeric_count")
            missing=rec.get("missing_or_nonfinite_count")
            if not all(isinstance(v,int) and not isinstance(v,bool) and 0<=v<=total for v in (finite,missing)) or finite+missing!=total:
                raise ValueError("inconsistent missing and numeric counts")
            neg=rec.get("negative_numeric_count")
            zero=rec.get("zero_numeric_count")
            if not all(isinstance(v,int) and 0<=v<=finite for v in (neg,zero)) or neg+zero>finite:
                raise ValueError("invalid negative or zero value frequencies")
            q=rec.get("value_quantiles_global") or {}
            if not all(k in q and isinstance(q[k],(int,float)) and math.isfinite(q[k]) for k in ("p05","p25","p50","p75","p95")):
                raise ValueError("missing aggregate quantiles")
            if [q[k] for k in ("p05","p25","p50","p75","p95")]!=sorted(q.values()):
                raise ValueError("unordered QC percentiles")
            ids=rec.get("feature_label_style_counts") or {}
            for key,val in ids.items():
                if not isinstance(val,int) or val<0 or val>n:
                    raise ValueError("feature format counts invalid")
            cleaned[omic]={
                "features":n,"biopsies":m,
                "nonfinite_values_n":missing,
                "nonfinite_fraction":round(missing/total,5),
                "negative_values_n":neg,
                "zeros_n":zero,
                "global_value_quantiles":{k:q[k] for k in ("p05","p25","p50","p75","p95")},
                "feature_label_style_counts":ids
            }
        overlap=obj.get("transcript_protein_shared_feature_label_count")
        if not isinstance(overlap,int) or overlap<0 or overlap>min(layers["T"]["feature_count"],layers["P"]["feature_count"]):
            raise ValueError("missing exact T/P feature-name overlap")
        result[fname]={"modalities":cleaned,
                       "transcript_protein_exact_feature_label_overlap":overlap,
                       "feature_annotation_independently_reviewed":False,
                       "biological_effect_tests_permitted":False}
    return {
        "schema":"hcc-ng-numeric-qc-v1",
        "sources_from_real_rds":True,
        "files":result,
        "status":"NUMERIC_DISTRIBUTION_AND_ID_OVERLAP_VERIFIED",
        "quantitative_scale_reviewed_against_S5_appendix":False,
        "imputation_provenance_audited":False,
        "V02":"PARTIAL",
        "RNA_PROTEIN_EFFECTS_NOT_EVALUATED":True,
        "no_individual_values_or_identifiers":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",required=True)
    ap.add_argument("--out",required=True)
    z=ap.parse_args()
    path=pathlib.Path(z.source)
    if not path.is_file() or path.stat().st_size>3_000_000:raise SystemExit("missing/suspicious input")
    r=check(json.loads(path.read_text(encoding="utf-8")))
    out=pathlib.Path(z.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,indent=2))

if __name__=="__main__":main()
