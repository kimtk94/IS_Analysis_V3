#!/usr/bin/env python3
"""Ephemeral Ng/Suter sample-ordered patient group IDs for paired-omics QC.

Strict SHA256-pinned Ng Supplementary Data 1 join; writes *only* 50 integers
to the private GitHub runner directory. No patient IDs/biopsy names are output
in logs or uploaded artifacts. Numbers are within-run group IDs, not clinical IDs.
"""
import argparse
import collections
import pathlib

from probe_ng_gao_public import _normalized_label
from probe_ng_suter_clinical_sample_join import official_rows, YES


def assign_patient_indexes(source_labels, rows, columns):
    if len(source_labels)!=50 or len(set(source_labels))!=50:
        raise ValueError("expected 50 source-biopsy labels")
    by_biopsy={}
    for row in rows:
        biopsy=_normalized_label(row.get(columns["tumorbiopsyid"],""))
        patient=_normalized_label(row.get(columns["patientid"],""))
        if not biopsy:continue
        if not patient:raise ValueError("source biopsy without patient")
        flags=[str(row.get(columns[f],"")).strip().lower() in YES
               for f in ("rnaseq","proteomics")]
        if biopsy in by_biopsy:
            raise ValueError("duplicate source biopsy")
        by_biopsy[biopsy]=(patient,all(flags))
    normalized=[_normalized_label(label) for label in source_labels]
    if len(set(normalized))!=50 or any(label not in by_biopsy for label in normalized):
        raise ValueError("RDS to source mismatch")
    if any(not by_biopsy[k][1] for k in normalized):
        raise ValueError("missing original assays")
    # Consecutive integers preserve source row order but not source IDs.
    compact={}
    output=[]
    for biopsy in normalized:
        patient=by_biopsy[biopsy][0]
        if patient not in compact:compact[patient]=len(compact)+1
        output.append(compact[patient])
    counts=collections.Counter(output)
    if len(counts)!=48 or len(output)!=50 or sum(n-1 for n in counts.values())!=2:
        raise ValueError("source 50 biopsy/48 patient structure changed")
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-labels",required=True)
    ap.add_argument("--private-patient-index-out",required=True)
    args=ap.parse_args()
    path=pathlib.Path(args.source_labels)
    if not path.is_file() or path.stat().st_size>20000:
        raise ValueError("ephemeral source sample labels missing")
    labels=[x for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
    source_rows,cols=official_rows()
    index=assign_patient_indexes(labels,source_rows,cols)
    dest=pathlib.Path(args.private_patient_index_out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text("".join(str(x)+"\n" for x in index),encoding="utf-8")
    dest.chmod(0o600)
    print("STAGE24_PATIENT_ASSIGNMENT_OK samples=50 independent_patients=48 extra_biopsies=2 IDs_REDACTED")

if __name__=="__main__":
    main()
