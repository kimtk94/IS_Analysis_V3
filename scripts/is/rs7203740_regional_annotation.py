#!/usr/bin/env python3
"""Post-readout rs7203740 locus annotation, not independent confirmation.
Preserve source failures and exact population identity. No causal gene claim.
"""
import argparse,hashlib,json,math,re,time
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import urlparse,urlencode
RS="rs7203740";CHR="16";POS=81775211
BASE="https://rest.ensembl.org"
POPULATIONS=["1000GENOMES:phase_3:CHB","1000GENOMES:phase_3:JPT","1000GENOMES:phase_3:CEU"]
def require(x,m):
    if not x:raise ValueError(m)
def get(url):
    meta={"url":url}
    try:
        require(urlparse(url).hostname in {"rest.ensembl.org","www.ebi.ac.uk"},"Unapproved host")
        with urlopen(Request(url,headers={"Accept":"application/json","Content-Type":"application/json","User-Agent":"IS-research-public-locus-audit/1.0"}),timeout=35) as r:
            require(urlparse(r.geturl()).hostname==urlparse(url).hostname,"Redirect host")
            buf=r.read(8000001);require(len(buf)<=8000000,"Response cap")
            meta.update(http_status=r.status,sha256=hashlib.sha256(buf).hexdigest(),bytes=len(buf))
        return json.loads(buf),meta
    except (OSError,ValueError,TypeError) as e:
        meta.update(status="HOLD_SOURCE",reason=str(e)[:220])
        if isinstance(e,HTTPError):meta["http_status"]=e.code
        return None,meta
def request(path,params=None):
    url=BASE+path+("?"+urlencode(params) if params else "")
    time.sleep(2)
    return get(url)
def audit(manifest):
    a=json.loads(manifest.read_text());require(a["verified_count"]==24 and len(a["rows"])==24,"Frozen anchor manifest")
    hits=[x for x in a["rows"].values() if x.get("rsid")==RS]
    require(len(hits)==1,"Anchor rsID missing")
    h=hits[0];require(h["status"]=="ANCHOR_VERIFIED" and str(h["chr"])==CHR and h["pos"]==POS,"Anchor identity")
    alleles={h["ref"],h["alt"]};require(len(alleles)==2 and all(re.fullmatch("[ACGT]",x) for x in alleles),"Anchor pair")
    out={"status":"POST_READOUT_REGIONAL_ANNOTATION","rsid":RS,"GRCh38_chr":CHR,"GRCh38_pos":POS,"historical_ref":h["ref"],"historical_alt":h["alt"],"source_run":45534891070,"sources":{},"identity_verified":False,"limits":["Variant selected for follow-up after the expanded-family result, not independent discovery. No inferred causal gene, colocalization, or independent replication.","LD is population-specific 1000 Genomes Phase 3 CHB/JPT/CEU reference evidence, not matched GWAS individual-level LD; no claim that an unreturned proxy is independent.","Catalog is literature-curated, not full GWAS summary statistics. Absence of exact rsID or regional association is not evidence of novelty."]}
    v,p=request("/variation/human/"+RS);out["sources"]["variation"]=p
    if v is None:return out
    try:
        require(v["name"]==RS,"Variation ID")
        mm=[m for m in v["mappings"] if m["assembly_name"]=="GRCh38" and str(m["seq_region_name"])==CHR and m["start"]==m["end"]==POS and m["strand"]==1]
        require(len(mm)==1 and alleles<=set(mm[0]["allele_string"].split("/")),"Ensembl coordinate/allele mismatch")
        out.update(identity_verified=True,genome_reference=mm[0]["allele_string"].split("/")[0],ensembl_allele_string=mm[0]["allele_string"],most_severe_consequence=v.get("most_severe_consequence"),minor_allele=v.get("minor_allele"),MAF=v.get("MAF"))
    except (KeyError,TypeError,ValueError) as e:out["identity_reason"]=str(e);return out
    v,p=request("/vep/human/id/"+RS,{"canonical":1});out["sources"]["vep"]=p
    if isinstance(v,list):
        try:
            rows=[x for x in v if x.get("id")==RS and x.get("assembly_name")=="GRCh38" and str(x.get("seq_region_name"))==CHR and x.get("start")==x.get("end")==POS]
            require(len(rows)==1,"Unexpected VEP coordinate or multiple rows")
            x=rows[0];require(alleles<=set(x["allele_string"].split("/")),"VEP pair")
            tc=[{k:t.get(k) for k in ("gene_id","gene_symbol","transcript_id","variant_allele","consequence_terms","biotype","canonical","impact")} for t in x.get("transcript_consequences",[]) if t.get("variant_allele") in alleles]
            out["vep"]={"most_severe_consequence":x.get("most_severe_consequence"),"transcript_consequences":tc,"regulatory_consequences":x.get("regulatory_feature_consequences",[])}
        except (KeyError,ValueError,TypeError) as e:out["vep"]={"status":"HOLD_SCHEMA","reason":str(e)}
    v,p=request(f"/overlap/region/human/{CHR}:{POS-250000}..{POS+250000}",{"feature":"gene"});out["sources"]["regional_genes"]=p
    if isinstance(v,list):
        try:
            genes=[]
            for g in v:
                require(g["assembly_name"]=="GRCh38" and str(g["seq_region_name"])==CHR and g["feature_type"]=="gene","Regional gene identity")
                start,end=g["start"],g["end"];require(start<=end and start<=POS+250000 and end>=POS-250000,"Gene coordinate")
                genes.append({"id":g["id"],"symbol":g.get("external_name"),"biotype":g.get("biotype"),"start":start,"end":end,"distance_bp":max(0,start-POS,POS-end)})
            out["regional_genes"]=sorted(genes,key=lambda g:(g["distance_bp"],g["id"]))
        except (KeyError,ValueError,TypeError) as e:out["regional_genes"]={"status":"HOLD_SCHEMA","reason":str(e)}
    out["ld"]={}
    for pop in POPULATIONS:
        v,p=request("/ld/human/"+RS+"/"+pop,{"window_size":500,"r2":0.8,"attribs":1});out["sources"]["ld_"+pop]=p
        if not isinstance(v,list):out["ld"][pop]={"status":"HOLD_SOURCE_OR_SCHEMA"};continue
        try:
            proxies={}
            for q in v:
                a,b=q["variation1"],q["variation2"]
                require(RS in (a,b) and q["population_name"]==pop,"LD identity")
                other=b if a==RS else a
                rr=float(q["r2"]);dd=float(q["d_prime"])
                require(.8-1e-8<=rr<=1 and 0<=abs(dd)<=1,"LD measure")
                if other==RS:continue
                require(other not in proxies or proxies[other]["r2"]==rr,"Duplicate LD conflict")
                proxies[other]={"rsid":other,"r2":rr,"d_prime":dd}
            out["ld"][pop]={"status":"REFERENCE_LD_RETRIEVED","n":len(proxies),"proxies":sorted(proxies.values(),key=lambda x:(-x["r2"],x["rsid"]))}
        except (KeyError,ValueError,TypeError) as e:out["ld"][pop]={"status":"HOLD_SCHEMA","reason":str(e)}
    # Official Catalog V2 exact variant resource, not the retired summary statistics API.
    v,p=get("https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms/"+RS)
    out["sources"]["catalog_exact_variant"]=p
    if v is not None:out["catalog_exact_variant"]=v
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--manifest",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x=audit(a.manifest);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    brief={**x}
    brief["ld"]={k:{**{a:b for a,b in v.items() if a!="proxies"},"first_10":v.get("proxies",[])[:10]} for k,v in x.get("ld",{}).items()}
    print(json.dumps(brief,indent=2,allow_nan=False))
if __name__=="__main__":main()
