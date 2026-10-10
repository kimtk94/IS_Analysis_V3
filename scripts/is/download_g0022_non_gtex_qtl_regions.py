#!/usr/bin/env python3
"""Conservative targeted tabix retrieval from official EBI indexed QTL files.

Downloads only chr12:111650000-112370000 for six predeclared source IDs,
one by one, verifies integrity, enforces 60MB per-file and 8min total budget.
No login/private access, no change to canonical research source.
"""
import argparse,csv,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
SOURCE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alternate_qtl_sources_v1/G0022_NON_GTEX_QTL_SOURCE_CANDIDATES.tsv")
OUT=Path("/srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1/non_gtex")
REGION="12:111650000-112370000"
def download(source,out):
    rows=list(csv.DictReader(source.open(),delimiter="\t"))
    if len(rows)!=6:raise ValueError("Study whitelist length changed")
    out.mkdir(parents=True,exist_ok=True)
    audit=[]
    for i,r in enumerate(rows):
        ds=r["dataset_id"];url=r["nominal_summary_uri"]
        if not url.startswith("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"):
            raise ValueError("Untrusted QTL download URL")
        output=out/f"{ds}_chr12_111650000_112370000.tsv.gz"
        temp=out/f"{ds}_chr12_111650000_112370000.tsv.gz.part"
        record={"study":r["study_label"],"sample_group":r["sample_group"],
                "dataset_id":ds,"region":REGION,"source_url":url,
                "source_samples":r["sample_size"],"status":"","row_count":0,
                "compressed_bytes":0,"sha256":"","error":""}
        if output.exists():
            record["status"]="SOURCE_ARCHIVE_PREEXISTING_VERIFY"
        else:
            try:
                if temp.exists():raise RuntimeError("Leftover partial file; inspect manually")
                result=subprocess.run(["tabix",url,REGION],capture_output=True,timeout=65)
                if result.returncode:
                    raise RuntimeError("tabix failed "+result.stderr.decode("utf-8","replace")[-350:])
                if len(result.stdout)<100:
                    record["status"]="EMPTY_QTL_REGION"
                    raise ValueError("No QTL rows in selected source region")
                payload=gzip.compress(result.stdout,compresslevel=1,mtime=0)
                if len(payload)>60*1024*1024:
                    raise ValueError("Requested region exceeded 60MB compressed")
                temp.write_bytes(payload)
                with gzip.open(temp,"rt") as f:
                    next(f) # validate gzip and no header expected
                temp.replace(output)
                record["status"]="DOWNLOADED_TARGET_REGION"
            except Exception as e:
                record["status"]="DOWNLOAD_OR_FORMAT_FAILED"
                record["error"]=str(e)
                if temp.exists():temp.unlink()
        if output.exists():
            h=hashlib.sha256()
            n=0
            with output.open("rb") as f:
                for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
            with gzip.open(output,"rb") as f:
                for line in f:n+=1
            record.update(row_count=n,compressed_bytes=output.stat().st_size,
                          sha256=h.hexdigest())
        audit.append(record)
        print(json.dumps(record),flush=True)
        if i<len(rows)-1:time.sleep(5)
    resultfile=out/"G0022_NON_GTEX_QTL_DOWNLOAD_MANIFEST.tsv"
    with resultfile.open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(audit[0]))
        w.writeheader();w.writerows(audit)
    return audit
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--sources",type=Path,default=SOURCE)
    ap.add_argument("--out",type=Path,default=OUT)
    a=ap.parse_args()
    download(a.sources,a.out)
