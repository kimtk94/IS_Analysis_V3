#!/usr/bin/env python3
"""Read-only, offline audit of 646 pre-existing Drive-backed GTEx coloc inputs.

Drive filename inventory is supplied externally; this script makes no
downloads, creates only derivative QC files, and performs no causal inference.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

REQUIRED={"locus","dataset_key","gene_base","eqtl_beta","eqtl_se",
          "gwas_beta","gwas_se","match_key","qtl_n_scalar","harmonization"}

def read(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def inspect(path,expected):
    """SNP schema/key and finite-stat QC, not direct coloc ABF replay."""
    n=0
    distinct=set()
    labels=set()
    try:
        with path.open(newline="",encoding="utf-8") as f:
            it=csv.DictReader(f,delimiter="\t")
            missing=REQUIRED-set(it.fieldnames or [])
            if missing:
                return False,0,"MISSING_COLUMNS:"+",".join(sorted(missing))
            for r in it:
                n+=1
                if (r["locus"],r["dataset_key"],r["gene_base"])!=(
                    expected["locus"],expected["dataset_key"],expected["gene_base"]):
                    return False,n,"SOURCE_ASSOCIATION_KEY_MISMATCH"
                v=r["match_key"]
                if not v or v in distinct:
                    return False,n,"DUPLICATE_OR_MISSING_VARIANT"
                distinct.add(v)
                labels.add(r["harmonization"])
                for k in ("eqtl_beta","gwas_beta","eqtl_se","gwas_se"):
                    try:
                        val=float(r[k])
                        if not math.isfinite(val) or (k.endswith("_se") and val<=0):
                            raise ValueError(k)
                    except (ValueError,TypeError):
                        return False,n,"NONFINITE_OR_INVALID_"+k
        if n!=int(expected["expected_nsnps"]):
            return False,n,"SNP_N_MISMATCH"
        if labels!={"EXACT_REF_ALT_GRCH38"}:
            return False,n,"UNEXPECTED_HARMONIZATION_LABELS"
        return True,n,"SCHEMA_SNP_COUNT_STATS_QC_PASS"
    except (OSError,UnicodeError,csv.Error) as e:
        return False,n,"IO_SCHEMA_ERROR_"+type(e).__name__

def audit(master,replay,drive,local_dir,recovered_dir):
    key=lambda r:(r["locus"],r["dataset_key"],r["gene_base"])
    if len(master)!=646 or len(replay)!=646:
        raise ValueError("Expected 646 original assays and 646 replay statuses")
    if len({key(r) for r in master})!=646 or len({key(r) for r in replay})!=646:
        raise ValueError("Duplicate assay key")
    if {key(r) for r in master}!={key(r) for r in replay}:
        raise ValueError("Archived assay sets differ")
    filenames={r["filename"] for r in replay}
    if len(filenames)!=646 or len(drive)!=646 or set(drive)!=filenames:
        raise ValueError("Drive inventory must contain exact 646 historical input filenames")
    output=[]
    for r in sorted(replay,key=key):
        name=r["filename"]
        original=local_dir/name
        copied=recovered_dir/name
        chosen=original if original.is_file() else copied if copied.is_file() else None
        if chosen:
            passed,n,detail=inspect(chosen,r)
            status=("LOCAL_SOURCE_QC_PASS" if chosen==original else "DRIVE_RECOVERED_QC_PASS") if passed else "INPUT_QC_REVIEW"
            digest=sha256(chosen)
            size=chosen.stat().st_size
        else:
            n=0;detail="DRIVE_FILENAME_IDENTIFIED_NOT_YET_QC_CHECKED"
            status="DRIVE_PRESENT_NOT_YET_COPIED"
            digest="";size=0
        output.append({
            "locus":r["locus"],"dataset_key":r["dataset_key"],"gene_base":r["gene_base"],
            "filename":name,"drive_inventory_present":True,
            "previous_replay_status":r["status"],
            "source_QC_status":status,"source_QC_reason":detail,
            "source_file":str(chosen) if chosen else "",
            "actual_snps":n if chosen else "",
            "expected_snps":r["expected_nsnps"],
            "source_file_bytes":size if chosen else "",
            "sha256":digest,
            "new_SNP_ABF_replay_status":"NOT_RUN",
            "interpretation":"INPUT_QC_ONLY_NOT_CAUSAL_INFERENCE",
        })
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--master",required=True,type=Path)
    p.add_argument("--replay",required=True,type=Path)
    p.add_argument("--drive-name-list",required=True,type=Path)
    p.add_argument("--local-dir",required=True,type=Path)
    p.add_argument("--recovered-dir",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    name_list=[x.strip().rstrip("/") for x in a.drive_name_list.read_text().splitlines() if x.strip()]
    output=audit(read(a.master),read(a.replay),name_list,a.local_dir,a.recovered_dir)
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Do not overwrite previous audit")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    dest=a.out_dir/"IS_646_COLOC_EXISTING_INPUT_RECOVERY.tsv"
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
        w.writeheader();w.writerows(output)
    summary={
        "status":"BACKUP_INVENTORY_AND_LOCAL_FILE_SCHEMA_ONLY",
        "assays":len(output),"drive_filenames":len(name_list),
        "source_status":dict(Counter(r["source_QC_status"] for r in output)),
        "original_replay_status":dict(Counter(r["previous_replay_status"] for r in output)),
        "SNP_coloc_ABF_rerun_now":0,
        "source_sha256":{"master":sha256(a.master),"replay":sha256(a.replay),
                         "drive_names":sha256(a.drive_name_list)},
        "external_database_downloaded":False,
        "original_master_modified":False,
        "GTEx_molecular_cohort_matched_LD_verified":False,
        "warning":"Copying SNP-level files and schema QC are not synonymous with an independently rerun coloc analysis.",
    }
    (a.out_dir/"IS_646_COLOC_RECOVERY_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("IS_646_DRIVE_SOURCE_AUDIT_COMPLETE",json.dumps(summary["source_status"],sort_keys=True))
    print("ORIGINAL_REPLAY",json.dumps(summary["original_replay_status"],sort_keys=True))
if __name__=="__main__": main()
