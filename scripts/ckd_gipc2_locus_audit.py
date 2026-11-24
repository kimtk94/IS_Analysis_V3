#!/usr/bin/env python3
"""Public GIPC2 locus audit: exact reagent rows, reported previous index SNP, LD and GTEx significant eQTLs."""
import argparse
import hashlib
import itertools
import json
import math
import re
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

PINS={"egfr":"ad0d9181d3d1b4d598b7c9e5a3d5be11825db8681a877069a4410c96834328c5",
      "bun":"a56ee97acc8ac60231a5ae68c39bcc8357c06651805721fed3302f4454a1cdf2"}
TARGETS={("Olink","GIPC2"),("SomaScan","seq.20921.1")}
PRIOR="rs650495"
DISCOVERY="rs11162405"
POPS=("EAS","EUR","AFR")

def read(path,trait):
    with path.open("rb") as f:
        if hashlib.file_digest(f,"sha256").hexdigest()!=PINS[trait]:raise ValueError("Archive digest")
    with zipfile.ZipFile(path) as z:
        if z.testzip():raise ValueError("Archive CRC")
        n=[x for x in z.namelist() if x.rsplit("/",1)[-1]=="assay_results.json"]
        if len(n)!=1:raise ValueError("Missing derived table")
        rows=json.loads(z.read(n[0]))
    d={(x["platform"],x["assay_id"]):x for x in rows if (x["platform"],x["assay_id"]) in TARGETS}
    if set(d)!=TARGETS:raise ValueError("Missing fixed paired reagent")
    return d

def get(base,query=None):
    url=base+("?" + urllib.parse.urlencode(query) if query else "")
    req=urllib.request.Request(url,headers={"Accept":"application/json","Content-Type":"application/json","User-Agent":"Public-research-locus-audit/1.0"})
    with urllib.request.urlopen(req,timeout=25) as f:
        raw=f.read(5000001)
    if len(raw)>5000000:raise ValueError("API cap")
    return json.loads(raw),hashlib.sha256(raw).hexdigest(),url

def resolve_variant(row):
    v=row["sentinel"].split(":")
    chrom,pos,ref,alt=v[0].removeprefix("chr"),int(v[1]),v[2],v[3]
    url=f"https://rest.ensembl.org/overlap/region/human/{chrom}:{pos}-{pos}"
    p,sha,uri=get(url,{"feature":"variation","content-type":"application/json"})
    ids=sorted({x["id"] for x in p if str(x["seq_region_name"])==chrom and x["start"]==x["end"]==pos
                and x.get("assembly_name")=="GRCh38" and set((ref,alt))<=set(x.get("alleles",[]))
                and re.fullmatch("rs[0-9]+",str(x.get("id")))})
    if len(ids)!=1:raise ValueError("No unique Ensembl exact coordinate allele rsID: "+str(ids))
    return ids[0],{"url":uri,"sha256":sha,"sentinel":row["sentinel"]}

def variation(rs):
    p,sha,u=get("https://rest.ensembl.org/variation/human/"+rs,{"content-type":"application/json"})
    matches=[m for m in p.get("mappings",[]) if m.get("assembly_name")=="GRCh38" and str(m.get("seq_region_name"))=="1" and m.get("start")==m.get("end") and m.get("strand")==1]
    # Restrict to a single forward-oriented biallelic chromosome 1 SNV.
    vv={(int(m["start"]),m["allele_string"]) for m in matches if re.fullmatch("[ACGT]/[ACGT]",m.get("allele_string",""))}
    if len(vv)!=1:raise ValueError("No unique GRCh38 biallelic mapping")
    pos,alleles=next(iter(vv))
    ref,alt=alleles.split("/")
    return {"rsid":rs,"variant_id":f"chr1_{pos}_{ref}_{alt}_b38","pos":pos,"alleles":alleles,"url":u,"sha256":sha}

def ld(a,b,pop):
    name="1000GENOMES:phase_3:"+pop
    url=f"https://rest.ensembl.org/ld/human/pairwise/{a}/{b}"
    data,sha,u=get(url,{"population_name":name,"content-type":"application/json"})
    rows=[r for r in data if {r.get("variation1"),r.get("variation2")}=={a,b} and r.get("population_name")==name]
    if len(rows)!=1:raise ValueError("No unique pairwise LD record; not zero LD")
    r=rows[0];r2=float(r["r2"]);dp=float(r["d_prime"])
    if not (0<=r2<=1 and 0<=dp<=1 and math.isfinite(r2) and math.isfinite(dp)):raise ValueError("LD bounds")
    return {"status":"REFERENCE_PAIRWISE_LD","rs1":a,"rs2":b,"population":pop,"r2":r2,"d_prime":dp,"sha256":sha,"url":u}

def gtex(v):
    allrows=[];page=0;hashes=[]
    while page<5:
        data,sha,u=get("https://gtexportal.org/api/v2/association/singleTissueEqtl",
            {"datasetId":"gtex_v10","variantId":v["variant_id"],"itemsPerPage":1000,"page":page})
        if not isinstance(data,dict) or not isinstance(data.get("data"),list):raise ValueError("GTEx response schema")
        rows=data["data"]
        for r in rows:
            if r.get("variantId")!=v["variant_id"]:raise ValueError("GTEx variant identity")
            if r.get("datasetId","gtex_v10")!="gtex_v10":raise ValueError("GTEx dataset identity")
        allrows.extend(rows);hashes.append({"url":u,"sha256":sha})
        pagination=data.get("paging_info",data.get("pagingInfo",{}))
        if len(rows)<1000:break
        page+=1
    else:raise ValueError("GTEx page limit")
    ids=[(r.get("gencodeId"),r.get("tissueSiteDetailId")) for r in allrows]
    if len(ids)!=len(set(ids)):raise ValueError("Duplicate GTEx rows")
    if any(not a or not b for a,b in ids):raise ValueError("GTEx identity schema")
    summary=[]
    for r in allrows:
        keys=("variantId","gencodeId","geneSymbol","tissueSiteDetailId","nes","pValue","snpId")
        summary.append({k:r[k] for k in keys if k in r})
    return {"status":"SIGNIFICANT_EQTL_ONLY","rsid":v["rsid"],"variant_id":v["variant_id"],
            "record_count":len(summary),"tissue_count":len({r["tissueSiteDetailId"] for r in summary}),
            "gene_ids":sorted({r["gencodeId"] for r in summary}),
            "kidney_rows":[r for r in summary if r["tissueSiteDetailId"] in ("Kidney_Cortex","Kidney_Medulla")],
            "target_rows":[r for r in summary if r["gencodeId"].split(".")[0]=="ENSG00000137960"],
            "all_rows":summary,"source_pages":hashes}

def safe(fn,*args):
    try:return fn(*args)
    except (OSError,ValueError,KeyError,TypeError,OverflowError) as e:return {"status":"HOLD","reason":str(e)[:400]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    e=read(a.root/"egfr.zip","egfr");b=read(a.root/"bun.zip","bun")
    reagent=[]
    fixed={PRIOR,DISCOVERY}
    for k,r in sorted(e.items()):
        s=b[k]
        if any(r.get(f)!=s.get(f) for f in ("sentinel","beta","se","F")):raise ValueError("Exposure varies by endpoint")
        v={"platform":k[0],"assay_id":k[1],"sentinel":r.get("sentinel"),"effect_allele":r.get("effect_allele"),
           "other_allele":r.get("other_allele"),"exposure_beta":r.get("beta"),"exposure_se":r.get("se"),"F":r.get("F"),"exposure_eaf":r.get("eaf"),"source_p":r.get("source_p")}
        for trait,x in (("egfr",r),("bun",s)):
            v[trait]={f:x.get(f) for f in ("status","rsid","mr_beta","mr_ci95","mr_p","bonferroni_p")}
        if r.get("rsid"):
            fixed.add(r["rsid"])
            v["resolved_rsid"]=r["rsid"]
        else:
            try:
                rs,provenance=resolve_variant(r);fixed.add(rs);v.update(resolved_rsid=rs,resolution=provenance)
            except (OSError,ValueError,KeyError,TypeError) as exc:v["resolution_hold"]=str(exc)[:300]
        reagent.append(v)
    variants=[{"rsid":rs,**safe(variation,rs)} for rs in sorted(fixed)]
    pairs=[]
    for x,y in itertools.combinations(sorted(fixed),2):
        for pop in POPS:pairs.append({"rs1":x,"rs2":y,"population":pop,**safe(ld,x,y,pop)})
    eqtls=[]
    for v in variants:
        if "variant_id" in v:eqtls.append({"rsid":v["rsid"],**safe(gtex,v)})
    result={"status":"GIPC2_LOCUS_EVIDENCE_AUDIT_NOT_COLOC","reagents":reagent,"variants":variants,
            "pairwise_reference_ld":pairs,"gtex_v10_significant_eqtls":eqtls,
            "limitations":["Ensembl reference LD is not cohort-specific and has no signed allele covariance here.",
            "Significant-only GTEx eQTL hits are not complete regional QTL data; missing is not null.",
            "Olink/SomaScan CKB participants and renal outcome are shared; no independent exposure replication.",
            "No evidence of protein mediation or novel causal variant follows from proximity or LD alone."]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({**result,"gtex_v10_significant_eqtls":[{k:v for k,v in q.items() if k!="all_rows"} for q in eqtls]},indent=2,allow_nan=False))
if __name__=="__main__":main()
