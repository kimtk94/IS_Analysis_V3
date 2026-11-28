#!/usr/bin/env python3
"""Qualify the processed Gao HCC independent-source archive without inferring specimen identity."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re,sys,tarfile
from pathlib import Path,PurePosixPath
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior/scripts/publication"))
import early_track_acceleration as base
base.HOSTS.add("zenodo.org");require=base.require
RECORD=14553766;NAME="copheemap_hcc.tar.gz";MD5="6b67ac84db568b9fc6b77950d7f06a70"
PANEL=("SOAT1","HMGCR","SQLE","LDLR","FASN","ACACA","ACLY","SREBF2")
ARCH=ROOT/"docs/publication/evidence/hcc_historical_gene_model_readback_20261125.json"
def prior():
    a=json.loads(ARCH.read_text());q=a["aggregate"]
    require(a["run_id"]==45613882235 and q["patients"]==9 and q["matched_tumor_biopsies"]==22 and q["unchanged_other_genes"]==7,"archived sensitivity identity")
    return q
def safe_member(m):
    p=PurePosixPath(m.name)
    require(not p.is_absolute() and ".." not in p.parts and (m.isdir() or m.isfile()),"unsafe or nonregular tar member")
    require(len(m.name)<=240 and m.size<=500_000_000,"member size/name cap")
def table_prefix(f):
    # Bounded UTF-8 text prefix; no sample identifiers or data values emitted.
    raw=f.read(131072)
    if raw[:2]==b"\x1f\x8b":return None
    txt=raw.decode("utf-8-sig",errors="strict")
    lines=txt.splitlines()
    if len(lines)<3:return None
    opts=[]
    for d in ("\t",","):
        rr=list(csv.reader(lines[:4],delimiter=d))
        if len(rr)>=3 and len(rr[0])>=2 and all(len(v)==len(rr[0]) for v in rr[:3]):opts.append((d,rr))
    if len(opts)!=1:return None
    d,rr=opts[0];h=rr[0]
    return d,h,rr[1:]
def keyshape(v):
    if re.fullmatch(r"ENSG\d+(?:\.\d+)?",v):return "Ensembl_gene"
    if re.fullmatch(r"\d+",v):return "integer"
    if re.fullmatch(r"[A-Za-z0-9_.-]+",v):return "alphanumeric"
    return "other"
def inventory():
    prior();meta,ms=base.get("https://zenodo.org/api/records/"+str(RECORD),2_000_000);d=json.loads(meta)
    require(d.get("id")==RECORD and d.get("doi")=="10.5281/zenodo.14553766","record identity")
    files=d.get("files",[]);hits=[v for v in files if v.get("key")==NAME]
    require(len(hits)==1 and hits[0]["checksum"]=="md5:"+MD5,"published checksum")
    n=hits[0]["size"];require(isinstance(n,int) and 0<n<=150_000_000,"archive compressed cap")
    raw,src=base.get(f"https://zenodo.org/records/{RECORD}/files/{NAME}?download=1",n)
    require(len(raw)==n and hashlib.md5(raw).hexdigest()==MD5,"archive integrity")
    reports=[];headers={};total=0;names=set()
    with tarfile.open(fileobj=io.BytesIO(raw),mode="r:gz") as t:
        for m in t:
            safe_member(m);require(m.name not in names,"duplicate tar member");names.add(m.name);require(len(names)<=80,"member count cap")
            if not m.isfile():continue
            total+=m.size;require(total<=1_200_000_000,"archive expanded cap")
            label=m.name
            z={"name":label,"bytes":m.size}
            if label.lower().endswith((".tsv",".txt",".csv",".gct")):
                f=t.extractfile(m)
                try:q=table_prefix(f)
                except (UnicodeError,csv.Error):q=None
                if q:
                    delim,h,rr=q
                    terms={"gene","gene_id","geneid","gene symbol","gene_symbol","symbol","hgnc","id","name","description","sample","sample_id","sampleid","patient","patient_id","tissue","sample_type","sampletype","aliquot","aliquot_id","specimen","specimen_id","entrezgene","entrez_gene_id","ensembl_gene_id"}
                    recognized=[{"column":i+1,"header":v} for i,v in enumerate(h) if v.lower().strip() in terms]
                    # Expression matrix candidate: explicit first feature key and numeric-backed remaining columns.
                    def numeric(v):
                        try:return math.isfinite(float(v))
                        except ValueError:return v in {"NA","NaN",""}
                    candidates=[i for i in range(1,len(h)) if all(numeric(v[i]) for v in rr[:2]) and h[i] and h[i].lower().strip() not in terms]
                    keys={h[i] for i in candidates}
                    z.update(delimiter=repr(delim),header_width=len(h),recognized_public_schema=recognized,first_feature_shapes=dict(collections.Counter(keyshape(v[0]) for v in rr)),
                             numeric_backed_columns=len(candidates),unique_candidate_keys=len(keys))
                    if len(candidates)>=100 and len(keys)==len(candidates):headers[label]=keys
            reports.append(z)
    pairs=[]
    for i,(a,x) in enumerate(headers.items()):
        for b,y in list(headers.items())[i+1:]:
            pairs.append({"file_a":a,"file_b":b,"exact_column_overlap":len(x&y),"columns_a":len(x),"columns_b":len(y)})
    return {"status":"HCC_GAO_PROCESSED_ARCHIVE_INVENTORY","record_source":ms,"source":src,"md5":MD5,"archive_members":len(names),"expanded_bytes":total,"files":reports,"candidate_exact_column_intersections":pairs,
            "limits":["Published Gao et al. Cell 2019 HBV HCC cohort of 159 paired patients, processed by CoPheeMap authors; source archive may contain tumor-only or mixed tissue data. Exact column overlap alone does not prove patient or tissue identity. No arbitrary ID normalization, individual-level export, raw extraction or genotype download. Need explicit sample annotation and exact gene-level mapping before fixed-panel analyses. Prior PDC numeric parity not a proven aliquot join. No independent claim relative to other Chinese cohorts without overlap check."]}
def selftest():
    prior();assert len(PANEL)==8
    x=tarfile.TarInfo("../bad")
    try:safe_member(x)
    except ValueError:pass
    else:raise AssertionError("unsafe path accepted")
    a=table_prefix(io.BytesIO(b"gene\tS1\tS2\nACACA\t1\t2\nSQLE\t3\t4\n"))
    assert a[0]=="\t" and len(a[1])==3
    print("Prior sensitivity, safe archive and explicit table-prefix checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["selftest","inventory"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=inventory()
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError,tarfile.TarError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    z["scientific_status"]="HOLD_SAMPLE_PROVENANCE_AND_NOVELTY"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
