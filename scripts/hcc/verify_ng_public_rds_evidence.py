#!/usr/bin/env python3
"""Validate Ng/Suter public RDS *aggregate-only* audit; no clinical claim."""
import argparse
import json
import pathlib
import re

SOURCE_COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
EXPECTED=("HCCfull.rds","HCCfullz.rds","HCCdata.rds")
MARKERS=("M","CN","T","P","PP")

def verify(doc):
    if doc.get("schema_version")!="hcc-ng-public-rds-audit-v1":
        raise ValueError("unexpected Ng provenance schema")
    if doc.get("source_commit")!=SOURCE_COMMIT:
        raise ValueError("public source revision was not pinned")
    if doc.get("raw_identifiers_or_numeric_matrices_in_artifact") is not False:
        raise ValueError("privacy guardrail absent")
    if doc.get("true_patient_by_gene_and_protein_quantitative_join_scientifically_validated") is not False:
        raise ValueError("scientific gate must not be promoted from geometry")
    if not doc.get("required_three_multimodal_objects_read"):
        raise ValueError("required public RDS objects unavailable")
    files=doc.get("files") or {}
    objects=doc.get("objects") or {}
    report={}
    for filename in EXPECTED:
        f=files.get(filename) or {}
        o=objects.get(filename) or {}
        if f.get("status")!="PUBLIC_PINNED_RDS_READ_IN_MEMORY":
            raise ValueError("required source file has not been read")
        if not re.fullmatch(r"[0-9a-f]{64}",str(f.get("sha256") or "")):
            raise ValueError("missing source SHA256")
        if o.get("status")!="PUBLIC_MULTIMODAL_RDS_SHAPES_AUDITED":
            raise ValueError("public RDS not recognized as a multilayer object")
        layers=o.get("layers") or {}
        if not all(v in layers for v in MARKERS):
            raise ValueError("expected five modalities absent")
        dimensions={}
        for k in MARKERS:
            t=layers[k]
            n, m=t.get("feature_count"),t.get("sample_count")
            if not isinstance(n,int) or not isinstance(m,int) or n<=0 or m<=0:
                raise ValueError("nonpositive omics matrix dimensions")
            if not t.get("numeric_matrix"):
                raise ValueError("layer is not numeric")
            if t.get("raw_labels_or_values_exported") is not False:
                raise ValueError("privacy guardrail absent in layer")
            expected_orientation="samples_by_features" if filename=="HCCdata.rds" else "features_by_samples"
            if t.get("orientation")!=expected_orientation:
                raise ValueError("source RDS modality axis orientation inconsistent")
            dimensions[k]={"features":n,"samples":m,"orientation":expected_orientation}
        same_n=(dimensions["T"]["samples"]==dimensions["P"]["samples"])
        names_equal=o.get("transcript_protein_identical_sample_labels") is True
        report[filename]={
            "layers":dimensions,
            "rna_protein_equal_column_counts":same_n,
            "rna_protein_shared_sample_labels_count":o.get("transcript_protein_shared_sample_labels"),
            "identical_RNA_protein_column_names":names_equal,
            "same_biopsy_quantitative_join_claimed":False,
            "status":"SOURCE_GEOMETRY_ALIGNABLE_NOT_CLINICALLY_VERIFIED" if same_n and names_equal
                     else "SOURCE_GEOMETRY_OR_LABEL_MISMATCH_REQUIRES_JOIN"
        }
    return {
        "execution":"VERIFIED_PUBLIC_QUANTITATIVE_RDS_STRUCTURE",
        "files":report,
        "V02":"BLOCKED_CLINICAL_SAMPLE_PROVENANCE_AND_QUANT_UNITS",
        "V04":"PARTIAL",
        "V12":"PARTIAL",
        "biological_claim_ready":False,
        "note":"Numeric multimodal matrix geometry is real; sample column name equality alone cannot certify patient linkage, exact biopsy origin, or independent etiology/recurrence analyses."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    args=p.parse_args()
    src=pathlib.Path(args.source)
    if not src.is_file() or src.stat().st_size>1_000_000:
        raise SystemExit("aggregate JSON missing or too large")
    report=verify(json.loads(src.read_text(encoding="utf-8")))
    print(json.dumps(report,indent=2))
if __name__=="__main__":
    main()
