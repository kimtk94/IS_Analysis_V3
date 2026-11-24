#!/usr/bin/env python3
"""Query official GTEx v8 collapsed exons using versioned GENCODE v26 IDs.

This only validates an expression gene-model TSS proxy. It cannot certify
the original FastQTL tested gene BED or support a causal gene/coloc claim.
"""
import argparse
import gzip
import hashlib
import json
import re
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz"
API="https://gtexportal.org/api/v2/dataset/collapsedGeneModelExon"
GENES={
 "FGF5":("ENSG00000138675","4",81266618,"T","G"),
 "NEURL":("ENSG00000107954","10",102493750,"T","C"),
 "C4orf22":("ENSG00000197826","4",None,None,None),
 "CALHM2":("ENSG00000138172","10",None,None,None),
}
MAX_GTF=85000000
MAX_API=1000000
ATTR=re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
VER=re.compile(r"ENSG\d{11}\.\d+")

def versioned_ids(gtf):
    x={}
    ids={v[0]:k for k,v in GENES.items()}
    with gzip.open(gtf,"rt",encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):continue
            cols=line.rstrip("\n").split("\t")
            if len(cols)!=9 or cols[2]!="gene":continue
            m=ATTR.search(cols[8])
            if not m or m.group(1).split(".")[0] not in ids:continue
            name=ids[m.group(1).split(".")[0]]
            if name in x or not VER.fullmatch(m.group(1)):
                raise ValueError("Invalid or duplicated GENCODE26 gene ID")
            chrom=cols[0].removeprefix("chr")
            strand=cols[6]
            if chrom!=GENES[name][1] or strand not in ("+","-"):
                raise ValueError("Invalid GENCODE gene chr/strand")
            a,b=int(cols[3]),int(cols[4])
            x[name]={"gencodeId":m.group(1),"strand":strand,"chromosome":chrom,
                     "regular_TSS":a if strand=="+" else b}
    if set(x)!=set(GENES):raise ValueError("Four versioned genes not all present")
    return x

def api_call(id,fetcher=None):
    if not VER.fullmatch(id):raise ValueError("API ID requires version")
    url=API+"?"+urllib.parse.urlencode({
        "datasetId":"gtex_v8","gencodeId":id,"page":0,"itemsPerPage":10000})
    if fetcher:
        raw=fetcher(url)
    else:
        req=urllib.request.Request(url,headers={
            "User-Agent":"IS-GTEx-v8-API-Reference-Validation/1.0",
            "Accept":"application/json"})
        with urllib.request.urlopen(req,timeout=24) as res:
            if urllib.parse.urlsplit(res.url).hostname!="gtexportal.org":
                raise ValueError("Untrusted API redirect")
            raw=res.read(MAX_API+1)
    if len(raw)>MAX_API:raise ValueError("Too many API bytes")
    return {"url":url,"sha256":hashlib.sha256(raw).hexdigest(),
            "response":json.loads(raw.decode("utf-8"))}

def extract_exons(z):
    found=[]
    def visit(v,depth):
        if depth>7:return
        if isinstance(v,list):
            for item in v:visit(item,depth+1)
        elif isinstance(v,dict):
            a=next((k for k in ("start","exonStart","chromStart") if k in v),None)
            b=next((k for k in ("end","exonEnd","chromEnd") if k in v),None)
            if a is not None and b is not None:
                try:
                    lo,hi=int(v[a]),int(v[b])
                    if lo>=1 and hi>=lo:
                        found.append({"start":lo,"end":hi,
                            "chromosome":str(v.get("chromosome",v.get("chrom",""))),
                            "strand":str(v.get("strand",""))})
                        return
                except (TypeError,ValueError):pass
            for child in v.values():
                if isinstance(child,(list,dict)):visit(child,depth+1)
    visit(z,0)
    return found

def interpret(name,gene,body):
    meta=GENES[name]
    payload=body["response"]
    entry={"gencodeId":gene["gencodeId"],"datasetId":"gtex_v8",
           "URL":body["url"],"response_sha256":body["sha256"],
           "top_level_keys":sorted(payload)[:35] if isinstance(payload,dict) else [],
           "exon_entries_found":0,"collapsed_TSS_approved":False}
    exons=extract_exons(payload)
    entry["exon_entries_found"]=len(exons)
    chrom={x["chromosome"].removeprefix("chr") for x in exons if x["chromosome"]}
    strands={x["strand"] for x in exons if x["strand"]}
    if not exons or chrom!={meta[1]} or strands!={gene["strand"]}:
        entry["gate"]="HOLD_EXON_SCHEMA_OR_CHR_STRAND"
        entry["exon_samples"]=exons[:3]
        return entry
    tss=(min(x["start"] for x in exons) if gene["strand"]=="+"
         else max(x["end"] for x in exons))
    entry["collapsed_TSS_approved"]=True
    entry["GTEx_v8_collapsed_model_TSS"]=tss
    entry["regular_GENCODE26_TSS"]=gene["regular_TSS"]
    entry["delta_TSS_bp"]=tss-gene["regular_TSS"]
    entry["collapsed_TSS_cis_lower"]=tss-1000000
    entry["collapsed_TSS_cis_upper"]=tss+1000000
    if meta[2] is not None:
        pos=meta[2]
        entry["QTD_only_variant"]=f"{meta[1]}:{pos}:{meta[3]}:{meta[4]}"
        entry["QTD_only_outside_collapsed_cis"]=not(tss-1000000<=pos<=tss+1000000)
        entry["outside_bp"]=max(tss-1000000-pos,pos-(tss+1000000),0)
    entry["gate"]="COLLAPSED_EXON_MODEL_TSS_VERIFIED_NOT_FASTQTL_GENE_BED"
    return entry

def read_bounded_gtf():
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/"gencode.v26.gtf.gz"
        request=urllib.request.Request(GTF_URL,headers={"User-Agent":"IS-GTEx-v8-GENCODE26-ID/1.0"})
        with urllib.request.urlopen(request,timeout=40) as resp,path.open("wb") as f:
            n=0
            while True:
                chunk=resp.read(1048576)
                if not chunk:break
                n+=len(chunk)
                if n>MAX_GTF:raise ValueError("GTF size cap")
                f.write(chunk)
        return versioned_ids(path)

def audit(ids=None,fetcher=None):
    out={"schema":"IS_GTExV8_API_CollapsedExon_TSS_v1",
         "GTEx_API_datasetId":"gtex_v8","GTEx_API_sequential_queries":True,
         "FastQTL_original_tested_gene_BED_TSS_verified":False,
         "source_variant_missingness_proven":False,"coloc_performed":False,
         "scientific_status":"HOLD"}
    try:
        ids=ids or read_bounded_gtf()
        out["GENCODE_v26_identifiers"]=ids
    except (OSError,ValueError,EOFError) as e:
        out.update({"api_gate":"HOLD_GENCODE26_VERSIONED_ID_UNAVAILABLE",
                    "failure_class":type(e).__name__})
        return out
    by={}
    for name in ("FGF5","NEURL","C4orf22","CALHM2"):
        try:
            response=api_call(ids[name]["gencodeId"],fetcher)
            by[name]=interpret(name,ids[name],response)
        except (OSError,ValueError,KeyError,UnicodeError) as exc:
            by[name]={"gate":"HOLD_API_REQUEST","failure_class":type(exc).__name__,
                      "reason":str(exc)[-150:]}
    out["by_gene"]=by
    valid=all(by[g].get("collapsed_TSS_approved") for g in GENES)
    outside=(all(by[g].get("QTD_only_outside_collapsed_cis") is True
                 for g in ("FGF5","NEURL")) if valid else None)
    out["all_four_gene_exon_models_parsed"]=valid
    out["two_missing_QTD_variants_outside_API_collapsed_TSS"]=outside
    out["api_gate"]=("FOUR_COLLAPSED_EXONS_VERIFIED_QTD_2_OUTSIDE_NOT_FASTQTL_CERTIFIED"
                     if valid and outside else "HOLD_API_SOURCE_INCOMPLETE_OR_DIFFERENT_BOUNDARIES")
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    result=audit()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("IS_GTEX_V8_COLLAPSED_API="+json.dumps(result,sort_keys=True,separators=(",",":")))
if __name__=="__main__":
    main()
