#!/usr/bin/env python3
"""Audit two fixed chromosome-13 genes and two SNPs in public brain vascular eQTLs.

Bryois Zenodo 7276971 contains all tested cis pairs, five documented columns
(no SE). Do not infer the effect allele or reconstruct SE from nominal p.
"""
import argparse,gzip,hashlib,json,math,re
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlsplit
REC="https://zenodo.org/api/records/7276971"
FILES=("Pericytes.13.gz","Endothelial.cells.13.gz")
GENES={"COL4A1_ENSG00000187498","COL4A2_ENSG00000134871"}
SNPS={"rs9515201","rs4771674"}
def check(ok,msg):
    if not ok:raise ValueError(msg)
def get(url,cap):
    check(urlsplit(url).scheme=="https" and urlsplit(url).hostname=="zenodo.org","Unexpected source URL")
    with urlopen(Request(url,headers={"User-Agent":"IS-public-summary-research/1.0","Accept-Encoding":"identity"}),timeout=70) as f:
        check(urlsplit(f.geturl()).hostname=="zenodo.org","Unexpected redirect")
        data=f.read(cap+1)
    check(0<len(data)<=cap,"Source exceeds byte bound or empty")
    return data
def scan(path):
    counts={g:0 for g in sorted(GENES)}
    hits={g:{s:[] for s in sorted(SNPS)} for g in sorted(GENES)}
    n=0; malformed=0; candidate_seen=set()
    with gzip.open(path,"rt",encoding="utf8",errors="strict") as f:
        for l in f:
            a=l.split();n+=1
            check(n<=15000000,"Row bound")
            if len(a)!=5:malformed+=1;continue
            g,s,d,p,b=a
            if g not in GENES:continue
            counts[g]+=1
            if s not in SNPS:continue
            key=(g,s)
            check(key not in candidate_seen,"Duplicate candidate SNP-gene pair")
            candidate_seen.add(key);p,b=float(p),float(b);d=int(d)
            check(math.isfinite(p) and 0<=p<=1 and math.isfinite(b),"Invalid candidate p/beta")
            hits[g][s].append({"nominal_p":p,"source_beta_unoriented":b,"distance_to_tss":d})
    check(n>1000 and malformed==0,"Documented five-column schema not confirmed")
    return {"rows":n,"gene_pair_counts":counts,"fixed_pairs":hits}
def run(root):
    raw=get(REC,4000000);rec=json.loads(raw)
    check(str(rec["id"])=="7276971","Wrong release")
    files={v["key"]:v for v in rec["files"]}
    out={"status":"PUBLIC_CIS_PAIR_AUDIT_NOT_COLOC","record":"7276971","record_sha256":hashlib.sha256(raw).hexdigest(),"files":{},"limits":["Only exact published gene-SNP rows and tested pair coverage; no cell-type interaction test.","Source has no SE and tested allele is in separate snp_pos.txt. Beta remains unoriented; do not compare direction to GWAS.","No new MR, colocalization, causal or ancestry claim. These data are not a matched East Asian vascular QTL cohort."]}
    root.mkdir(parents=True,exist_ok=True)
    for name in FILES:
        try:
            v=files[name];size=int(v["size"]);md5=v["checksum"].removeprefix("md5:")
            check(0<size<=50000000 and re.fullmatch("[0-9a-f]{32}",md5),"Unverified manifest")
            url="https://zenodo.org/records/7276971/files/"+name+"?download=1"
            buf=get(url,size)
            check(len(buf)==size and hashlib.md5(buf).hexdigest()==md5,"Full publisher file integrity mismatch")
            p=root/name;p.write_bytes(buf)
            out["files"][name]={"status":"SCANNED","bytes":size,"md5":md5,"sha256":hashlib.sha256(buf).hexdigest(),**scan(p)}
        except (OSError,ValueError,KeyError,TypeError,EOFError) as e:out["files"][name]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:200]}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.root)
    except (OSError,ValueError,KeyError,TypeError) as e:x={"status":"HOLD_RECORD","reason":str(e)[:220]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
