#!/usr/bin/env python3
"""Bounded empirical triage of two previously underexplored Master's tracks.

Oral: retrieve the public GEO *processed* matrix and inspect actual shape.
Skin: inspect bounded authentic GWAS summary-statistics bytes, schema/numerics.
Never publish cell/sample identities, variant IDs, raw matrix or GWAS records.
Source-only success does not establish a valid research hypothesis or clinical
validation. HOLD on missing source or unknown phenotype orientation.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import pathlib
import tempfile
import urllib.request
from collections import Counter

URL_ORAL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz"
URL_SKIN="https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/perceived_age_summary_results.txt"
EXPECTED_SKIN={"SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM_INF","LOG_OR","LOG_SE","N_0","N_0.5","N_1","N_EFFECTIVE"}

def fetch_prefix(url,limit=4_000_000):
    r=urllib.request.Request(url,headers={"User-Agent":"MASTER_DEGREE_PUBLIC_RESEARCH_TRIAGE/1.0","Range":f"bytes=0-{limit-1}"})
    with urllib.request.urlopen(r,timeout=50) as response:
        first=response.read(limit)
        status=response.status
        length=response.headers.get("Content-Length")
    return first,{"http":status,"content_length":length,"observed_prefix_bytes":len(first),
                  "prefix_sha256":hashlib.sha256(first).hexdigest()}

def inspect_skin(blob):
    raw=blob.decode("utf-8",errors="replace")
    lines=raw.splitlines()
    if len(lines)<6: raise ValueError("too few complete lines")
    header=lines[0].lstrip("#").split()
    if len(set(header))!=len(header):raise ValueError("duplicate column names")
    missing=EXPECTED_SKIN-set(header)
    if missing: raise ValueError("skin required header missing "+str(sorted(missing)))
    ix={v:i for i,v in enumerate(header)}
    total=0;valid=0;pal=0;sex=Counter()
    for line in lines[1:-1]:
        if not line.strip() or line.startswith("#"):continue
        fields=line.split()
        if len(fields)!=len(header):continue
        total+=1
        try:
            a,b=fields[ix["ALLELE1"]],fields[ix["ALLELE0"]]
            freq=float(fields[ix["A1FREQ"]]);p=float(fields[ix["P_BOLT_LMM_INF"]])
            effect=float(fields[ix["LOG_OR"]]);se=float(fields[ix["LOG_SE"]])
            ok=(a.upper() in "ACGT" and b.upper() in "ACGT" and
                len(a)==1 and len(b)==1 and a!=b and 0<=freq<=1 and
                0<=p<=1 and se>0 and all(map(lambda x:abs(x)<1e100,(effect,se))))
            if ok: valid+=1
            if {a.upper(),b.upper()} in ({"A","T"},{"C","G"}):pal+=1
        except (ValueError,IndexError,OverflowError):pass
    if total<50: raise ValueError("insufficient complete GWAS records in bounded prefix")
    return {"header_sha256":hashlib.sha256(lines[0].encode()).hexdigest(),
            "columns":header,"full_source_read":False,
            "complete_rows_checked_from_prefix":total,"rows_basic_numeric_valid":valid,
            "palindromic_rows_prefix_only":pal,"fraction_valid":valid/total,
            "caution":"BOLT-LMM perceived-age LOG_OR describes coded perceived age; clinical skin-aging and proper MR phenotype validity unverified"}

def inspect_oral(path):
    n=0;first_shape=None;tot_uncompressed=0;row_columns=Counter()
    with gzip.open(path,"rb") as f:
        for line in f:
            n+=1;tot_uncompressed+=len(line)
            if tot_uncompressed>3_000_000_000: raise ValueError("uncompressed cap exceeded")
            if n<=4:
                ncols=line.count(b"\t")+1
                row_columns[str(ncols)]+=1
                if n==1:first_shape=ncols
    if n<5 or first_shape is None: raise ValueError("no data rows")
    return {"source_gzip_bytes":path.stat().st_size,"source_gzip_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "decompressed_lines":n,"first_line_field_count":first_shape,
            "first_four_line_field_counts":dict(row_columns),
            "uncompressed_bytes_scanned":tot_uncompressed,
            "gene_by_cell_orientation":"NOT_YET_AUTHENTICATED",
            "donor_mapping":"NOT_YET_AUTHENTICATED",
            "biological_interpretation":"HOLD_PATIENT_LEVEL_PSEUDOBULK_UNTIL_CELL_METADATA_JOIN"}

def oral():
    with tempfile.TemporaryDirectory() as t:
        target=pathlib.Path(t)/"gse103322.txt.gz"
        h=hashlib.sha256();n=0
        req=urllib.request.Request(URL_ORAL,headers={"User-Agent":"MASTER_DEGREE_PUBLIC_RESEARCH_TRIAGE/1.0"})
        with urllib.request.urlopen(req,timeout=60) as response, target.open("wb") as fo:
            while True:
                chunk=response.read(1024*1024)
                if not chunk:break
                n+=len(chunk)
                if n>125_000_000:raise ValueError("GEO 125 MB cap exceeded")
                h.update(chunk);fo.write(chunk)
        if n<100_000:raise ValueError("GEO transferred too few bytes")
        result=inspect_oral(target)
        result["download_sha256_matches_scanned"]=result["source_gzip_sha256"]==h.hexdigest()
        if not result["download_sha256_matches_scanned"]:raise ValueError("hash drift")
        return result

def selftest():
    sample="SNP CHR BP ALLELE1 ALLELE0 A1FREQ P_BOLT_LMM_INF LOG_OR LOG_SE N_0 N_0.5 N_1 N_EFFECTIVE\n"
    for i in range(80):sample+=f"rs{i} 1 {i+1} A G 0.1 0.5 0.01 0.2 10 0 1 11\n"
    r=inspect_skin(sample.encode())
    assert r["complete_rows_checked_from_prefix"]==79
    assert r["rows_basic_numeric_valid"]==79
    with tempfile.TemporaryDirectory() as t:
        p=pathlib.Path(t)/"x.gz"
        with gzip.open(p,"wt") as f:f.write("X\tY\n"+"\n".join(["1\t2"]*8)+"\n")
        assert inspect_oral(p)["decompressed_lines"]==9
    print("RESEARCH_DECISION_REALDATA_SELFTEST_PASS")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--track",required=False,choices=["oral","skin"])
    p.add_argument("--selftest",action="store_true")
    args=p.parse_args()
    if args.selftest: selftest();return
    if not args.track:p.error("--track required")
    out=pathlib.Path("results/research_operations/decision_gate_20261102")
    out.mkdir(parents=True,exist_ok=True)
    result={"track":args.track,"scientific_gate":"HOLD","source_access":"HOLD","primary_source":URL_ORAL if args.track=="oral" else URL_SKIN,
            "full_source_identity_or_phenotype_verified":False}
    try:
        if args.track=="oral":
            result["evidence"]=oral()
            result["source_access"]="PASS_ACTUAL_PROCESSED_GEO_GZIP_SCANNED"
        else:
            raw,metadata=fetch_prefix(URL_SKIN)
            result["evidence"]=inspect_skin(raw)
            result["http_probe"]=metadata
            result["source_access"]="PASS_BOUNDED_REAL_GWAS_NUMERIC_PREFIX"
    except Exception as exc:
        # Missing or inaccessible data is not a negative biological result.
        result["blocker"]=type(exc).__name__+":"+str(exc)[:170]
    (out/(args.track+".json")).write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("RESEARCH_DECISION_GATE_ACTUAL_SOURCE",json.dumps({"track":result["track"],
      "source_access":result["source_access"],"scientific_gate":result["scientific_gate"],
      "evidence":result.get("evidence"),"blocker":result.get("blocker")},ensure_ascii=False))

if __name__=="__main__":
    main()
