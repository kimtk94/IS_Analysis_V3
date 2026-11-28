#!/usr/bin/env python3
"""Recover nominal QTL single-position screen after remote TabixFile failures.

Read only verified official index/header and one bounded BGZF span per source.
Prioritize monocyte and blood; no bulk download or blind retries.
"""
import argparse,hashlib,importlib,io,json,re,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"indexed_helpers/scripts/is"))
import fixed_subtype_indexed as http
from qtl_original_indexed_http_range import parse_bgzf_range
prev=importlib.import_module("scripts.is.rs7203740_exact_allele_qtl")
meta=importlib.import_module("scripts.is.audit_official_p3a_qtl_catalogue_provenance")
COLS=prev.COLS
POS=prev.POS
TARGETS=("QTD000021","QTD000356")

def require(x,m):
    if not x:raise ValueError(m)

def read_source(url):
    idx,_=http.get(url+".tbi",8000000)
    t=http.parse_index(idx)
    require(t["seq"]==2 and t["begin"]==3 and t["end"] in (0,3),"QTL index coordinate fields changed")
    time.sleep(4)
    first,cr=http.get(url,65536,(0,65535))
    total=int(cr.rsplit("/",1)[1])
    h=parse_bgzf_range(first,0,0)["data_lines"]
    firstrow=h.readline().lstrip("#").rstrip("\r\n").split("\t")
    if set(COLS)<=set(firstrow):
        header=firstrow;origin="ORIGINAL_HEADER"
    else:
        require(len(firstrow)==len(COLS) and firstrow[1].removeprefix("chr").isdigit() and firstrow[2].isdigit() and firstrow[5]==f"chr{firstrow[1].removeprefix('chr')}_{firstrow[2]}_{firstrow[3]}_{firstrow[4]}","Original header absent and pinned column order not corroborated by data")
        header=COLS;origin="PINNED_COLUMN_ORDER_WITH_FIRST_ROW_CHECK"
    require(len(set(header))==len(header) and header[1:3]==["chromosome","position"],"QTL schema")
    span=http.indexed_span(t,"16",POS,POS,total)
    require(span is not None,"No indexed chunk for exact position")
    a,b,x,y=span
    # Explicit hard cap below the generic helper's 12 MB; no full-file fallback.
    require(b-a+1<=4000000,"Single-position compressed budget exceeded")
    time.sleep(4)
    body,cr=http.get(url,b-a+1,(a,b))
    decoded=parse_bgzf_range(body,x,y)
    relevant=[];n=0
    for line in decoded["data_lines"]:
        n+=1;require(n<=500000,"Decoded row cap")
        fields=line.rstrip("\r\n").split("\t")
        require(len(fields)==len(header),"Unexpected width inside indexed block")
        if fields[1].removeprefix("chr")=="16" and fields[2]==str(POS):
            relevant.append(line.rstrip("\r\n"))
    require(len(relevant)<=10000,"Exact position row cap")
    x=prev.parse(relevant,header)
    return {**x,"status":"EXACT_POSITION_READ_HTTP_RANGE","header_origin":origin,"index_sha256":http.sha(idx),"header_range_sha256":http.sha(first),"source_range_sha256":http.sha(body),"source_range":cr,"source_total_bytes":total,"decoded_rows":n,"original_full_file_downloaded":False}

def run():
    src=meta.get_official_sources();md=meta.table(src["meta"]);paths=meta.table(src["paths"])
    require("GRCh38 alternative allele (also the effect allele)" in src["columns"],"QTL orientation documentation missing")
    out={"status":"EXACT_ALLELE_HTTP_RECOVERY","GRCh38_position":"16:81775211","rsid":"rs7203740","target":"C/T","other_allele":"C/A (separate descriptive control)","datasets":{},"limits":["Recovery of two of the five previously failed preselected sources, prioritized for monocyte and blood; not five-source coverage.","Index/header plus one exact indexed range per source; source range integrity only, not publisher full-file checksum.","Same single-position post-readout QTL screen, not regional colocalization or MR. Nominal N and AN/2 are not certified regression N/sdY. No EAS ancestry-matched inference. Zero returned rows from a failed source is HOLD."]}
    for i,k in enumerate(TARGETS):
        if i:time.sleep(6)
        study,label,group,n=prev.DATASETS[k];a,b=md[k],paths[k]
        for field,value in {"study_id":study,"study_label":label,"sample_group":group,"sample_size":n,"quant_method":"ge"}.items():require(a[field]==b[field]==value,"Pinned metadata mismatch "+k+" "+field)
        canonical=f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{study}/{k}/{k}.all.tsv.gz"
        require(b["ftp_path"]==canonical,"Unexpected source path")
        url=canonical.replace("ftp://","https://",1);d={"source":url,"study":label,"sample_group":group,"nominal_dataset_n":int(n),"pmid":a["pmid"]}
        try:d.update(read_source(url))
        except (OSError,ValueError,KeyError,TypeError,IndexError) as e:d.update(status="HOLD_SOURCE_OR_SCHEMA",reason=str(e)[:220])
        out["datasets"][k]=d
    out["sources_read"]=sum(x["status"]=="EXACT_POSITION_READ_HTTP_RANGE" for x in out["datasets"].values())
    out["exact_CT_candidate_rows"]=sum(v["alt"]=="T" for d in out["datasets"].values() for v in d.get("candidate_rows",[]))
    out["exact_CA_candidate_rows"]=sum(v["alt"]=="A" for d in out["datasets"].values() for v in d.get("candidate_rows",[]))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run()
    except (OSError,ValueError,KeyError,TypeError) as e:x={"status":"HOLD_PINNED_SOURCE","reason":str(e)[:220]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
