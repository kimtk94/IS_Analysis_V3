#!/usr/bin/env python3
"""P1 public paired 10x matrix integrity; aggregate-only output, no raw files in Git."""
import argparse,csv,gzip,hashlib,io,json,re,tempfile,urllib.parse,urllib.request
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone

SAMPLES={"vehicle":("GSM9373278","pool_C"),"dex":("GSM9373287","pool_dexC")}
BASE="https://ftp.ncbi.nlm.nih.gov/geo/"
def file_url(gsm,pool,suffix):
    return BASE+"samples/"+gsm[:-3]+"nnn/"+gsm+"/suppl/"+gsm+"_"+pool+"_"+suffix
def download(url,path,cap):
    req=urllib.request.Request(url,headers={"User-Agent":"ocular-research-p1/1.0"})
    digest=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=75) as r,open(path,"wb") as f:
        n=r.headers.get("Content-Length")
        if n and int(n)>cap:raise ValueError("CONTENT_LENGTH_OVER_CAP "+n)
        while True:
            b=r.read(2*1024*1024)
            if not b:break
            size+=len(b)
            if size>cap:raise ValueError("TRANSFER_OVER_CAP")
            digest.update(b);f.write(b)
    with open(path,"rb") as f:
        if f.read(2)!=bytes.fromhex("1f8b"):raise ValueError("BAD_GZIP_MAGIC")
    return {"bytes":size,"sha256":digest.hexdigest()}
def bar(s):
    return re.sub(r"-\d+$","",str(s).strip()).upper()
def lines(p):
    with gzip.open(p,"rt",encoding="utf8",errors="replace") as f:
        return [v.strip() for v in f if v.strip()]
def matrix_check(path,nfeat,nbar,full=True):
    with gzip.open(path,"rb") as f:
        h=f.readline().strip()
        if h!=b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("WRONG_MATRIX_FORMAT "+repr(h)[:120])
        while True:
            line=f.readline()
            if not line:raise ValueError("MISSING_MATRIX_DIMS")
            if not line.startswith(b"%"):break
        dims=tuple(map(int,line.split()))
        if len(dims)!=3:raise ValueError("WRONG_MATRIX_DIMENSIONS")
        rows,cols,nnz=dims
        if rows!=nfeat or cols!=nbar or nnz<1:
            raise ValueError("10X_DIMENSION_MISMATCH "+str(dims)+" versus "+str((nfeat,nbar)))
        readlines=0
        if full:
            # Fast buffered block scan validates sparse entry count and gzip CRC.
            carry=b""
            while True:
                data=f.read(4*1024*1024)
                if not data:break
                readlines+=data.count(b"\n")
                carry=(carry+data)[-1:]
            if carry and carry!=b"\n":readlines+=1
            if readlines!=nnz:raise ValueError("SPARSE_NNZ_MISMATCH "+str((readlines,nnz)))
        else:
            if len(f.readline().split())<3:raise ValueError("NO_SPARSE_COORDINATES")
    return {"n_features":rows,"n_barcodes":cols,"declared_nnz":nnz,
            "verified_nnz":readlines if full else None,"gzip_crc_checked":full}
def metadata_summary(path):
    grouped={}
    with gzip.open(path,"rt",encoding="utf8",errors="replace") as f:
        rd=csv.DictReader(f,delimiter="\t")
        required={"library_id","barcode","donor_id_demuxlet","DEX Treatment"}
        if not required.issubset(rd.fieldnames or []):
            raise ValueError("METADATA_REQUIRED_COLUMNS_MISSING")
        for r in rd:
            k=(r["library_id"].strip(),r["DEX Treatment"].strip().lower())
            grouped.setdefault(k,{})[bar(r["barcode"])]=r["donor_id_demuxlet"].strip()
    return grouped
def demux_barcodes(path):
    seen=set()
    with gzip.open(path,"rt",encoding="utf8",errors="replace") as f:
        rd=csv.DictReader(f,delimiter="\t")
        if "BARCODE" not in (rd.fieldnames or []):raise ValueError("DEMUX_BARCODE_COLUMN_MISSING")
        for r in rd:
            if r.get("BARCODE"):seen.add(bar(r["BARCODE"]))
    return seen
def paired_probe(name,gsm,pool,meta,tmp):
    out={"gsm":gsm,"condition":name,"files":{},"status":"HOLD"}
    try:
        f={}
        for suffix in ("barcodes.tsv.gz","features.tsv.gz","dmx.best.gz"):
            fp=tmp/(name+"_"+suffix)
            out["files"][suffix]=download(file_url(gsm,pool,suffix),fp,12*1024*1024)
            f[suffix]=fp
        barcodes=lines(f["barcodes.tsv.gz"]);features=lines(f["features.tsv.gz"])
        bset={bar(x) for x in barcodes}
        if len(bset)!=len(barcodes) or not features:raise ValueError("BARCODES_DUPLICATE_OR_EMPTY")
        out["n_barcodes"]=len(barcodes);out["n_features"]=len(features)
        dmx=demux_barcodes(f["dmx.best.gz"])
        out["demux_fraction_of_barcodes"]=round(len(bset&dmx)/len(bset),4)
        candidates=[]
        for (lib,condition),mapping in meta.items():
            if (name=="dex")!=(condition=="dex"):continue
            n=len(bset&set(mapping))
            candidates.append((n,lib,condition,len(mapping)))
        if not candidates:raise ValueError("MISSING_METADATA_CONDITION")
        candidates.sort(reverse=True)
        n,lib,condition,libsize=candidates[0]
        out["metadata_library"]=lib
        out["metadata_condition"]=condition
        out["metadata_joined_barcodes"]=n
        out["metadata_total_barcodes_in_library"]=libsize
        out["metadata_join_fraction"]=round(n/max(1,libsize),4)
        out["demux_join_fraction"]=round(len(bset&dmx)/len(bset),4)
        donors={d for bc,d in meta[(lib,condition)].items()
                if bc in bset and d and d.lower() not in ("na","unassigned","doublet")}
        out["n_joined_donors"]=len(donors)
        if out["metadata_join_fraction"]<.90:raise ValueError("METADATA_BARCODE_JOIN_LT_90pct")
        if out["demux_join_fraction"]<.90:raise ValueError("DEMUX_BARCODE_JOIN_LT_90pct")
        matrix=tmp/(name+"_matrix.mtx.gz")
        out["files"]["matrix.mtx.gz"]=download(file_url(gsm,pool,"matrix.mtx.gz"),matrix,620*1024*1024)
        out["matrix"]=matrix_check(matrix,len(features),len(barcodes),full=True)
        out["status"]="PASS_MATRIX_AND_BARCODE_JOIN"
        return out,donors
    except Exception as e:
        out["reason"]=type(e).__name__+": "+str(e)[:220]
        return out,set()
    finally:
        (tmp/(name+"_matrix.mtx.gz")).unlink(missing_ok=True)
def ena_registry():
    args={"accession":"PRJNA1243143","result":"read_run",
          "fields":"run_accession,sample_accession,study_accession,library_layout,fastq_ftp",
          "format":"tsv"}
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?"+urllib.parse.urlencode(args)
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"ocular-research-p1/1.0"})
        with urllib.request.urlopen(req,timeout=45) as r:
            data=r.read(2*1024*1024+1)
        if len(data)>2*1024*1024:raise ValueError("RESPONSE_TOO_LARGE")
        rd=csv.DictReader(io.StringIO(data.decode("utf8","replace")),delimiter="\t")
        runs=list(rd)
        if not runs:return {"status":"HOLD","reason":"EMPTY_ENA_REGISTRY"}
        if not {"run_accession","sample_accession"}.issubset(rd.fieldnames or []):
            raise ValueError("INVALID_ENA_COLUMNS")
        return {"status":"RUN_METADATA_PASS","run_count":len(runs),
                "biosample_count":len({r["sample_accession"] for r in runs}),
                "layout_counts":dict(Counter(r.get("library_layout","") for r in runs)),
                "fastq_endpoints":sum(bool(r.get("fastq_ftp")) for r in runs),
                "donor_treatment_mapping":"NOT_VALIDATED"}
    except Exception as e:
        return {"status":"HOLD","reason":type(e).__name__+": "+str(e)[:200]}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    a={"utc":datetime.now(timezone.utc).isoformat(),"dataset":"GSE313736",
       "study_scope":"in-vitro TM steroid response, NOT clinical drug-induced IOP",
       "pairs":{},"status":"HOLD"}
    with tempfile.TemporaryDirectory(prefix="ocular-p1-") as t:
        tmp=Path(t);meta={}
        try:
            fp=tmp/"cell_metadata.tsv.gz"
            a["metadata_file"]=download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",fp,12*1024*1024)
            meta=metadata_summary(fp)
            a["metadata_cells"]=sum(len(v) for v in meta.values())
        except Exception as e:
            a["metadata_error"]=type(e).__name__+": "+str(e)[:160]
        donor={}
        for name,(gsm,pool) in SAMPLES.items():
            x,donor[name]=paired_probe(name,gsm,pool,meta,tmp)
            a["pairs"][name]=x
        a["aggregate_pair_donors"]={
            "control":len(donor.get("vehicle",set())),
            "dex":len(donor.get("dex",set())),
            "paired":len(donor.get("vehicle",set())&donor.get("dex",set()))}
    a["independent_replication_PRJNA1243143"]=ena_registry()
    if all(x["status"]=="PASS_MATRIX_AND_BARCODE_JOIN" for x in a["pairs"].values()) and a["aggregate_pair_donors"]["paired"]>0:
        a["status"]="GO_POOL_C_FULL_COUNT_MATRIX_INTEGRITY"
    else:
        a["status"]="HOLD_POOL_C_INTEGRITY"
    Path(args.output).write_text(json.dumps(a,indent=2,ensure_ascii=False)+"\n")
    print("FINAL",a["status"],"donors",a["aggregate_pair_donors"])
    for name,x in a["pairs"].items():
        print("SAMPLE",name,x["status"],"join",x.get("metadata_join_fraction"),"demux",x.get("demux_join_fraction"),"reason",x.get("reason",""))
    print("REPLICATION",a["independent_replication_PRJNA1243143"]["status"])
    # Data QC failures produce HOLD evidence, without claiming the workflow itself passed QC.
if __name__=="__main__":main()
