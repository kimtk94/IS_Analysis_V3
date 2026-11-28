#!/usr/bin/env python3
"""Scan exact published Japanese immune-cell conditional QTL summary rows.
No donor-level records, no significant-only coloc, no inferred effect allele.
"""
import argparse, collections, gzip, hashlib, io, json, math, re, tarfile, urllib.request, zipfile
from pathlib import Path

GENES=("FGF5","ALDH2","SH3PXD2A","COL4A2","COL4A1","CALHM2","NEURL1","C4orf22","INA")
BASE="https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-398/"
ZIP_SIZE=564882709
ZIP_MD5="2fc0a61db661d179d1083bb39ec605de"
TAR_SIZE=4133191680
TAR_MD5="83236d666d3ba135b0af726932b73fb2"
MAX_SELECTED_BYTES=100_000_000

def get(url,cap=5000000):
    req=urllib.request.Request(url,headers={"Accept":"application/json, */*","Content-Type":"application/json","User-Agent":"IS-public-research/1.0"})
    with urllib.request.urlopen(req,timeout=30) as f:
        raw=f.read(cap+1)
    if len(raw)>cap:raise ValueError("Bounded response exceeded")
    return raw

def digest(p,algo):
    with p.open("rb") as f:return hashlib.file_digest(f,algo).hexdigest()

def write(out,o):
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(o,indent=2,allow_nan=False)+"\n")
    print(json.dumps(o,indent=2,allow_nan=False),flush=True)

def resolve_genes():
    result={}
    for gene in GENES:
        # NEURL1 is historical symbol; retain explicit alias, never fuzzy match.
        query="NEURL1" if gene=="NEURL1" else gene
        url="https://rest.ensembl.org/lookup/symbol/homo_sapiens/"+query+"?content-type=application/json"
        try:
            raw=get(url)
            r=json.loads(raw)
            if r.get("species")!="homo_sapiens" or r.get("assembly_name")!="GRCh38" or not re.fullmatch(r"ENSG[0-9]+",str(r.get("id"))):
                raise ValueError("Unexpected Ensembl identity")
            if r.get("display_name")!=gene and not (gene=="NEURL1" and r.get("display_name")=="NEURL"):
                raise ValueError("Symbol drift: "+str(r.get("display_name")))
            result[gene]={"id":r["id"],"symbol":r["display_name"],"chr":str(r["seq_region_name"]),
                          "url":url,"sha256":hashlib.sha256(raw).hexdigest()}
        except (OSError,ValueError,KeyError,TypeError) as e:
            result[gene]={"status":"HOLD_GENE_ID","reason":str(e)[:200]}
    if result["SH3PXD2A"].get("id")!="ENSG00000107957":
        raise ValueError("Primary SH3PXD2A identity unavailable")
    return result

def unpack(root):
    p=root/"conditional.zip"; target=root/"conditional.tar"
    if p.stat().st_size!=ZIP_SIZE or digest(p,"md5")!=ZIP_MD5:raise ValueError("DDBJ ZIP MD5/size mismatch")
    with zipfile.ZipFile(p) as z:
        m=[x for x in z.infolist() if x.filename.rsplit("/",1)[-1]=="conditional_eQTL_FDR0.05.tar"]
        if len(m)!=1 or m[0].file_size!=TAR_SIZE:raise ValueError("Unexpected published nested tar")
        h=hashlib.md5()
        with z.open(m[0]) as f,target.open("wb") as w:
            while b:=f.read(4*1024*1024):
                w.write(b);h.update(b)
        if target.stat().st_size!=TAR_SIZE or h.hexdigest()!=TAR_MD5:raise ValueError("Nested tar MD5/size mismatch")
    return target

def gene_tokens(fields, mapping):
    return sorted({mapping.get(f.split(b".")[0]) or mapping.get(f) for f in fields if f.split(b".")[0] in mapping or f in mapping})

def scan(root,out):
    gene_map=resolve_genes()
    mapping={}
    for g,d in gene_map.items():
        if "id" in d:
            for k in (d["id"],d["symbol"],g):
                if k.encode() in mapping and mapping[k.encode()]!=g:raise ValueError("Gene identity collision")
                mapping[k.encode()]=g
    # Restrict to exact token boundaries; also support GENCODE version suffix.
    pat=re.compile(rb"(?<![A-Za-z0-9_])(?:"+b"|".join(re.escape(x) for x in sorted(mapping,key=len,reverse=True))+rb")(?:\.[0-9]+)?(?![A-Za-z0-9_])")
    tar=unpack(root)
    inventory=[]; all_counts=collections.Counter(); totals=0; chosen_bytes=0
    selected=out.parent/"conditional_candidate_rows.tsv.gz"
    out.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(tar,"r:") as t,gzip.open(selected,"wb") as z:
        z.write(b"member\toriginal_summary_row\n")
        seen=set()
        for m in t:
            if not m.isfile():continue
            if m.name in seen:raise ValueError("Duplicate tar member")
            seen.add(m.name)
            if len(seen)>300:raise ValueError("Too many tar members")
            f=t.extractfile(m)
            if f is None:raise ValueError("Unreadable member")
            stream=io.BufferedReader(f)
            if stream.peek(2)[:2]==b"\x1f\x8b":stream=gzip.GzipFile(fileobj=stream)
            n=0;counts=collections.Counter();widths=collections.Counter();cols=collections.defaultdict(collections.Counter);preview=collections.defaultdict(list);first=[]
            for line in stream:
                n+=1;totals+=1
                if len(line)>100000 or b"\x00" in line:raise ValueError("Unexpected binary/non-summary row in "+m.name)
                if n<=2:first.append(line.decode("utf-8",errors="strict").strip()[:1200])
                if not pat.search(line):continue
                fs=line.split()
                genes=gene_tokens(fs,mapping)
                if not genes:continue
                widths[len(fs)]+=1
                for g in genes:
                    counts[g]+=1;all_counts[g]+=1
                    for k,fld in enumerate(fs):
                        if mapping.get(fld.split(b".")[0],mapping.get(fld))==g:cols[g][str(k+1)]+=1
                    if len(preview[g])<2:preview[g].append(line.decode("utf-8",errors="strict").strip()[:1500])
                chosen_bytes+=len(line)+len(m.name)+1
                if chosen_bytes>MAX_SELECTED_BYTES:raise ValueError("Selected rows exceed bounded storage; refine explicit schema")
                z.write(m.name.encode()+b"\t"+line.rstrip(b"\r\n")+b"\n")
            inventory.append({"member":m.name,"bytes_in_tar":m.size,"lines":n,"first_lines":first,
                "selected_rows_by_gene":dict(counts),"candidate_column_positions":{g:dict(x) for g,x in cols.items()},
                "selected_row_widths":dict(widths),"candidate_examples":dict(preview)})
            print("Scanned member",m.name,n,dict(counts),flush=True)
    if not inventory or totals<1000:raise ValueError("No meaningful source scan")
    summary={"status":"EXACT_CANDIDATE_TOKEN_SCAN_CONDITIONAL_SOURCE_ONLY",
        "source":"E-GEAD-398","zip_md5":ZIP_MD5,"tar_md5":TAR_MD5,"zip_sha256":digest(root/"conditional.zip","sha256"),
        "genes":gene_map,"total_source_lines":totals,"member_count":len(inventory),
        "candidate_rows_by_gene":dict(all_counts),"selected_bytes":chosen_bytes,"members":inventory,
        "limits":["The publisher conditional tables are significance-filtered; counts are not independent signals.",
                  "Member names and actual column identities must be verified before interpreting slopes or P values.",
                  "No tested-variant denominator, complete cis regional data, allelic orientation, or formal Coloc inferred.",
                  "Same donors across subsets; no cell-specific effect difference test or stroke disease context inferred."]}
    # Save full header and examples but keep Actions log concise.
    out.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps({**summary,"members":[{k:v for k,v in m.items() if k not in ("candidate_examples",)} for m in inventory]},indent=2,allow_nan=False),flush=True)

def metadata(out):
    records=[]
    for id in ("21296030","17190663"):
        try:
            raw=get("https://zenodo.org/api/records/"+id)
            p=json.loads(raw);m=p["metadata"]
            if str(p["id"])!=id:raise ValueError("Record identity drift")
            records.append({"record":id,"status":"METADATA","title":m.get("title"),"doi":m.get("doi"),"access_right":m.get("access_right"),
                "license":m.get("license"),"description":str(m.get("description",""))[:3000],
                "sha256":hashlib.sha256(raw).hexdigest(),
                "files":[{k:v for k,v in x.items() if k in ("key","size","checksum")} for x in p.get("files",[])]})
        except (OSError,ValueError,KeyError,TypeError) as e:
            records.append({"record":id,"status":"HOLD","reason":str(e)[:300]})
    write(out,{"status":"MAEEA_FULL_NOMINAL_SOURCE_INVENTORY_ONLY","records":records})

def selftest():
    m={b"ENSG00000107957":"SH3PXD2A",b"SH3PXD2A":"SH3PXD2A",b"COL4A2":"COL4A2"}
    assert gene_tokens(b"ENSG00000107957.9 chr10 123 456".split(),m)==["SH3PXD2A"]
    assert not gene_tokens(b"ENSG000001079570 SH3PXD2A2 COL4A20".split(),m)
    assert gene_tokens(b"COL4A2 SH3PXD2A".split(),m)==["COL4A2","SH3PXD2A"]
    print("Three exact identity tests passed")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=["selftest","scan","metadata"])
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--out",type=Path,default=Path("outputs/result.json"))
    a=p.parse_args()
    if a.mode=="selftest":selftest()
    elif a.mode=="metadata":metadata(a.out)
    else:scan(a.root,a.out)
if __name__=="__main__":main()
