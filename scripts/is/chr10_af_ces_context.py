#!/usr/bin/env python3
"""Fixed chr10 AF/CES locus context. No locus-specific mediation or causal inference."""
import argparse,collections,hashlib,itertools,json,math,re,urllib.parse,urllib.request
from pathlib import Path
from statistics import NormalDist

# Freeze published stroke and AF anchors, not selected from today's GTEx outputs.
IDS=("rs12415501","rs4918058","rs55983834","rs6584555","rs35176054")
GENES={"NEURL","NEURL1","SH3PXD2A","CALHM2"}
POPS=("EAS","EUR","JPT")
TISSUES=("Heart_Atrial_Appendage","Heart_Left_Ventricle","Artery_Aorta","Artery_Coronary","Artery_Tibial","Whole_Blood")
def require(ok,msg):
    if not ok:raise ValueError(msg)
def get(url,args=None):
    if args:url+="?"+urllib.parse.urlencode(args)
    req=urllib.request.Request(url,headers={"Accept":"application/json","Content-Type":"application/json","User-Agent":"Public-locus-study/1.0"})
    with urllib.request.urlopen(req,timeout=15) as f:
        raw=f.read(6000001)
    require(len(raw)<=6000000,"Oversized API response")
    return json.loads(raw),hashlib.sha256(raw).hexdigest(),url
def safe(fn,*args):
    try:return fn(*args)
    except (OSError,ValueError,KeyError,TypeError,OverflowError) as e:return {"status":"HOLD","reason":str(e)[:400]}
def variant(rs):
    data,sha,url=get("https://rest.ensembl.org/variation/human/"+rs)
    rows={(int(x["start"]),x["allele_string"]) for x in data["mappings"] if x.get("assembly_name")=="GRCh38" and str(x.get("seq_region_name"))=="10" and x.get("start")==x.get("end") and x.get("strand")==1 and re.fullmatch("[ACGT]/[ACGT]",x.get("allele_string",""))}
    require(len(rows)==1,"No unique biallelic forward GRCh38 mapping")
    pos,alleles=next(iter(rows));ref,alt=alleles.split("/")
    require(ref!=alt and 102403919<=pos<=104460593,"Unexpected chr10 position")
    return {"status":"EXACT_REFERENCE_VARIANT","rsid":rs,"pos38":pos,"ref":ref,"alt":alt,"variant_id":f"chr10_{pos}_{ref}_{alt}_b38","source_sha256":sha,"url":url}
def ld(x,y,pop):
    name="1000GENOMES:phase_3:"+pop
    data,sha,url=get(f"https://rest.ensembl.org/ld/human/pairwise/{x}/{y}",{"population_name":name,"content-type":"application/json"})
    hits=[r for r in data if {r.get("variation1"),r.get("variation2")}=={x,y} and r.get("population_name")==name]
    require(len(hits)==1,"No unique LD record; not zero")
    r=float(hits[0]["r2"]);dp=float(hits[0]["d_prime"])
    require(math.isfinite(r) and math.isfinite(dp) and 0<=r<=1 and 0<=dp<=1,"Invalid LD")
    return {"status":"UNSIGNED_REFERENCE_LD","r2":r,"d_prime":dp,"source_sha256":sha,"url":url}
def eqtl(v):
    rows=[];sources=[]
    for page in range(5):
        data,sha,url=get("https://gtexportal.org/api/v2/association/singleTissueEqtl",{"datasetId":"gtex_v10","variantId":v["variant_id"],"itemsPerPage":1000,"page":page})
        require(isinstance(data,dict) and isinstance(data.get("data"),list),"GTEx schema")
        part=data["data"]
        for r in part:
            require(r.get("variantId")==v["variant_id"] and r.get("datasetId","gtex_v10")=="gtex_v10","GTEx identity")
            require(bool(r.get("gencodeId")) and bool(r.get("tissueSiteDetailId")),"Missing QTL identity")
            row={k:r[k] for k in ("variantId","gencodeId","geneSymbol","tissueSiteDetailId","nes","pValue") if k in r}
            rows.append(row)
        sources.append({"sha256":sha,"url":url})
        if len(part)<1000:break
    else:raise ValueError("GTEx page limit")
    require(len({(r["gencodeId"],r["tissueSiteDetailId"]) for r in rows})==len(rows),"Duplicate GTEx records")
    by_tissue={t:[r for r in rows if r["tissueSiteDetailId"]==t] for t in TISSUES}
    return {"status":"SIGNIFICANT_GTEX_V10_ONLY","rsid":v["rsid"],"n_records":len(rows),"n_genes":len({r["gencodeId"] for r in rows}),
            "n_tissues":len({r["tissueSiteDetailId"] for r in rows}),"selected_tissue_all_genes":by_tissue,
            "candidate_gene_all_tissues":[r for r in rows if r.get("geneSymbol") in GENES],
            "all_rows":rows,"source_pages":sources}
def reference(out):
    variants=[{"rsid":x,**safe(variant,x)} for x in IDS]
    # Six fixed comparisons, chosen by published AF and stroke anchor roles.
    pairs=(("rs12415501","rs4918058"),("rs12415501","rs55983834"),("rs4918058","rs55983834"),
           ("rs12415501","rs6584555"),("rs6584555","rs4918058"),("rs35176054","rs55983834"))
    results=[]
    for x,y in pairs:
        for pop in POPS:
            results.append({"rs1":x,"rs2":y,"population":pop,**safe(ld,x,y,pop)})
    eq=[{"rsid":v["rsid"],**safe(eqtl,v)} for v in variants if v["status"]=="EXACT_REFERENCE_VARIANT"]
    r={"status":"REFERENCE_LD_AND_SIGNIFICANT_QTL_NOT_MEDIATION","fixed_variants":variants,"ld":results,"gtex":eq,
       "limits":["Reference LD is unsigned and may not match discovery cohorts; absent API output is not zero LD.",
                 "GTEx significant-only output is not full nominal statistics, not Japanese eQTL or evidence of cell specificity.",
                 "Do not infer stroke or AF causal effects, mediation, or independence from r2 alone."]}
    write(out,r,{"status":r["status"],"fixed_variants":variants,"ld":results,"gtex":[{k:v for k,v in q.items() if k!="all_rows"} for q in eq],"limits":r["limits"]})
def write(out,r,log=None):
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print(json.dumps(r if log is None else log,indent=2,allow_nan=False),flush=True)
def workbook(root,out):
    from openpyxl import load_workbook
    results=[]
    for year in (2023,2025):
        p=root/f"af{year}.xlsx"
        with p.open("rb") as f:sha=hashlib.file_digest(f,"sha256").hexdigest()
        w=load_workbook(p,read_only=True,data_only=True)
        inventory=[];matches=[]
        for s in w.worksheets:
            header=[];n=0;found=0
            for i,row in enumerate(s.iter_rows(values_only=True),1):
                n+=1
                fields=["" if x is None else str(x) for x in row]
                if i<=6:header.append(fields[:65])
                # Exact fixed rsID or gene tokens. Do not assume every matching
                # row is an AF association or that an rsID is a lead.
                ids=set();genes=set()
                for f in fields:
                    tokens=set(re.split(r"[\s,;|]+",f.strip()))
                    ids.update(tokens.intersection(IDS))
                    genes.update(tokens.intersection(GENES))
                if ids or genes:
                    found+=1
                    if len(matches)>=2500:raise ValueError("Unbounded candidate rows")
                    matches.append({"sheet":s.title,"row":i,"rsids":sorted(ids),"genes":sorted(genes),"fields":fields[:100]})
            inventory.append({"sheet":s.title,"rows":n,"matched_rows":found,"header":header})
        w.close()
        results.append({"year":year,"source_sha256":sha,"bytes":p.stat().st_size,"sheets":inventory,"matching_rows":matches})
    # Save exact source fields. Log matching rows and headers of hit sheets only.
    r={"status":"PUBLISHED_AF_SOURCE_CONTEXT_NOT_NEW_GWAS","results":results,
       "limit":"Published selected tables; interpret source headers and cohort/conditional models before comparing. Not full association statistics, not independent replication or a new statistical test."}
    write(out,r,{"status":r["status"],"results":[{**s,"sheets":[q for q in s["sheets"] if q["matched_rows"]]} for s in results],"limit":r["limit"]})
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["reference","workbook"]);p.add_argument("--root",type=Path,default=Path("."));p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    if a.mode=="reference":reference(a.out)
    else:workbook(a.root,a.out)
if __name__=="__main__":main()
