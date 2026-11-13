#!/usr/bin/env python3
"""GLGC indexed EUR remote-window probe using pysam/HTSlib (fallback backend).

No full GWAS is downloaded or validated. No disease inference is run.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path

from scripts.ckm_stroke.glgc_region_probe import BASE, FILE_NAMES, validate_rows


def probe(trait, chrom="22", start=20000000, end=20050000, row_cap=20000, factory=None):
    if trait not in FILE_NAMES or chrom not in {str(i) for i in range(1, 23)}:
        raise ValueError("invalid phenotype or chromosome")
    if not (0 < start < end) or row_cap < 2:
        raise ValueError("invalid region or row cap")
    url=BASE+FILE_NAMES[trait]
    if factory is None:
        try:
            import pysam
        except ImportError as exc:
            return {"status":"TOOL_UNAVAILABLE","reason":str(exc),"scientific_status":"NOT_RUN","url":url}
        factory=pysam.TabixFile
    try:
        with factory(url) as f:
            headers=list(f.header)
            data=[]
            capped=False
            for row in f.fetch(chrom, start-1, end):
                if len(data)>=row_cap:
                    capped=True
                    break
                data.append(row)
        if not headers:
            return {"status":"HOLD","reason":"tabix source header unavailable","url":url,
                    "observed_data_rows":len(data),"scientific_status":"NOT_RUN"}
        result=validate_rows("\n".join(headers+data)+"\n",chrom,start,end)
        if capped:
            result["status"]="HOLD"
            result["reason"]="row cap reached; sampled region incomplete"
        result.update({"backend":"pysam_tabix", "scope":"INDEXED_REGION_ONLY_NOT_FULL_GWAS",
                       "source":"GLGC_Graham2021","trait":trait,"ancestry":"EUR",
                       "region":f"{chrom}:{start}-{end}","url":url,
                       "index_url":url+".tbi","scientific_status":"NOT_RUN",
                       "row_cap":row_cap,"capped":capped})
        return result
    except (OSError,ValueError,KeyError,RuntimeError) as exc:
        return {"status":"SOURCE_UNAVAILABLE","reason":str(exc)[:400],
                "backend":"pysam_tabix","scientific_status":"NOT_RUN","url":url}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--trait",required=True,choices=sorted(FILE_NAMES))
    p.add_argument("--out",required=True)
    p.add_argument("--chrom",default="22")
    p.add_argument("--start",type=int,default=20000000)
    p.add_argument("--end",type=int,default=20050000)
    args=p.parse_args()
    result=probe(args.trait,args.chrom,args.start,args.end)
    result["timestamp_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    result["source_commit"]=os.getenv("GITHUB_SHA","")
    out=Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
    # A completed CI process with SOURCE_UNAVAILABLE used to be a false PASS.
    return int(result["status"] != "PASS")


if __name__=="__main__":
    raise SystemExit(main())
