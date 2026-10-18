#!/usr/bin/env python3
"""Read-only inventory of public BBJ ischemic stroke archive; no extraction."""
import argparse
import csv
import io
import json
import zipfile
from pathlib import Path

def inspect(path:Path):
    if not path.is_file():return {"status":"SOURCE_ZIP_NOT_PRESENT","entries":[]}
    if not zipfile.is_zipfile(path):raise ValueError("invalid zip")
    entries=[]
    with zipfile.ZipFile(path) as z:
        for info in z.infolist():
            if len(entries)>=250:raise ValueError("too many archive members")
            name=info.filename
            if info.is_dir():continue
            if name.startswith("/") or ".." in Path(name).parts or "\\" in name:
                raise ValueError("unsafe member name")
            record={"name":name,"compressed_bytes":info.compress_size,
                    "uncompressed_bytes":info.file_size,"file_type":"unknown"}
            ext=name.lower()
            if ext.endswith((".tsv",".txt",".csv",".tsv.gz",".txt.gz",".csv.gz")):
                record["file_type"]="possible_sumstats"
                # Member names and sizes are the only approved audit surface at
                # this stage. No raw rows or participant-level material are output.
            entries.append(record)
    return {"schema":"is-bbj-zip-schema-inventory-v1","status":"ZIP_MEMBERS_LISTED",
            "entries":entries,"scientific_identity":"NOT_VERIFIED",
            "case_control_metadata_adopted":False,
            "archive_allele_comparison_performed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    res=inspect(a.zip)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(res,ensure_ascii=False,indent=2)+"\n")
    print("BBJ_ZIP_SCHEMA_AUDIT",json.dumps(res,ensure_ascii=False))
if __name__=="__main__":main()
