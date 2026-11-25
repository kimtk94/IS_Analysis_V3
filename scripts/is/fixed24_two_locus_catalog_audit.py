#!/usr/bin/env python3
"""Bounded GWAS Catalog v2 regional and exact high-LD proxy audit.
Never infer novelty from absence in the curated catalog.
"""
import argparse,collections,hashlib,json,math,re,time
from pathlib import Path
from urllib.parse import urlencode,urlparse
from urllib.request import Request,urlopen
BASE="https://www.ebi.ac.uk/gwas/rest/api/v2/"
WINDOWS={"rs2767060":("10",35833363,500000),"rs7203740":("16",81775211,500000)}
PAGE=100
MAX_PAGES=6
MAX_EXACT=12
POPS={"1000GENOMES:phase_3:CHB":24,"1000GENOMES:phase_3:JPT":34,"1000GENOMES:phase_3:CEU":59}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def fetch(endpoint,params=None):
    u=BASE+endpoint+("?"+urlencode(params) if params else "")
    time.sleep(2)
    req=Request(u,headers={"Accept":"application/json","User-Agent":"IS-research-regional-curated-audit/1.0"})
    with urlopen(req,timeout=40) as r:
        require(r.status==200 and urlparse(r.geturl()).hostname=="www.ebi.ac.uk","Catalog HTTP/host")
        raw=r.read(8000001);require(len(raw)<=8000000,"Catalog response cap")
    return json.loads(raw),{"url":u,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)}
def paged(endpoint,params):
    out=[];pages=[];total=None;expected_pages=None
    for i in range(MAX_PAGES):
        x,source=fetch(endpoint,{**params,"page":i,"size":PAGE})
        require(isinstance(x,dict) and isinstance(x.get("content"),list) and isinstance(x.get("page"),dict),"Unexpected documented v2 page schema")
        p=x["page"];content=x["content"]
        n,size,tp,num=[p[k] for k in ("total_elements","size","total_pages","number")]
        require(all(isinstance(v,int) and not isinstance(v,bool) for v in (n,size,tp,num)) and n>=0 and tp>=0 and 0<size<=PAGE and num==i,"Page identity")
        require(len(content)<=size and tp==math.ceil(n/size),"Page count")
        if total is None:total=n;expected_pages=tp
        require(n==total and tp==expected_pages,"Live pagination changed")
        out.extend(content);pages.append(source)
        if i+1>=tp:break
    return out,{"total_elements":total,"returned_elements":len(out),"complete":len(out)==total and len(pages)==max(1,expected_pages),"pages":pages}
def read_proxy(path):
    d=json.loads(path.read_text())
    require(d["status"]=="POST_READOUT_LD_PROXY_PUBLICATION_AUDIT" and d["rsid"]=="rs2767060" and d["reference_populations_parsed"]==3,"Prior LD artifact identity")
    require(d["workbook_sha256"]=="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4" and d["publication_matches"]==[],"Prior published workbook result")
    proxies={}
    for pop,n in POPS.items():
        a=d["populations"][pop]
        require(a["status"]=="REFERENCE_LD_PARSED" and a["two_current_reads_identical_after_canonicalization"] is True and len(a["proxies"])==n,"Complete verified proxy list")
        seen=set()
        for v in a["proxies"]:
            rs=v["rsid"];require(rs not in seen and re.fullmatch(r"rs[0-9]+",rs) and .8-1e-8<=v["r2"]<=1 and v["chr"]=="10","Proxy identity")
            seen.add(rs);proxies.setdefault(rs,{})[pop]=v["r2"]
    require(len(proxies)==d["proxy_union"]==78 and sum(len(x)==3 for x in proxies.values())==d["proxy_shared_all_three"]==15,"LD counts")
    return proxies
def audit(proxy_path):
    proxy=read_proxy(proxy_path)
    out={"status":"POST_READOUT_CURATED_CATALOG_REGIONAL_AUDIT","documentation":"https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/","documented_openapi_git_blob":"52640a5408d477853fd1d09ededd86470ec28e89","regions":{},"exact_associations":{},"limits":["Catalog v2 literature-curated associations, not complete GWAS summary statistics. Regional variant membership is not evidence of stroke association; absent entries do not establish novelty.","GRCh38 +/-500 kb windows fixed before Catalog read. Exact index/proxy matches only; no proxy set for multiallelic rs7203740 (earlier LD API empty). Reference LD not GWAS matched.","Bounded pagination and exact lookup cap; incomplete retrieval stays HOLD. No unreturned variant is interpreted as absent. No independent replication or causal gene inference."]}
    md,src=fetch("metadata");require(isinstance(md,dict),"Catalog metadata")
    build=md.get("gene_build");require(isinstance(build,str) and build.startswith("GRCh38"),"Catalog build not GRCh38")
    out["catalog_metadata"]={k:md.get(k) for k in ("version","data_release_date","api_release_date","gene_build","dbsnp_build","ensembl_version")}
    out["catalog_metadata_source"]=src
    exact=set(WINDOWS)
    for rs,(ch,pos,win) in WINDOWS.items():
        try:
            rows,pages=paged("single-nucleotide-polymorphisms",{"chromosome":ch,"bp_start":pos-win,"bp_end":pos+win,"sort":"rs_id","direction":"asc"})
            ids=set();summ=[]
            for v in rows:
                k=v["rs_id"];require(isinstance(k,str) and k not in ids,"Duplicate or invalid catalog variant ID");ids.add(k)
                loc=[q for q in v["locations"] if str(q.get("chromosome_name"))==ch and isinstance(q.get("chromosome_position"),int) and pos-win<=q["chromosome_position"]<=pos+win]
                require(bool(loc),"Returned variant outside requested GRCh38 region")
                summ.append({"rsid":k,"positions":[q["chromosome_position"] for q in loc],"mapped_genes":v.get("mapped_genes"),"alleles":v.get("alleles"),"functional_class":v.get("functional_class")})
            targets={rs}|(set(proxy) if rs=="rs2767060" else set())
            hits=sorted(ids&targets);exact.update(hits)
            out["regions"][rs]={"status":"COMPLETE_CURATED_REGION" if pages["complete"] else "HOLD_PAGE_CAP","chr":ch,"start":pos-win,"end":pos+win,"pagination":pages,"exact_index_or_proxy_hits":hits,"regional_variants":summ}
        except (OSError,ValueError,KeyError,TypeError,IndexError) as e:out["regions"][rs]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:240]}
    out["planned_exact_queries"]=sorted(exact)
    if len(exact)>MAX_EXACT:
        out["exact_query_status"]="HOLD_EXACT_QUERY_BUDGET";return out
    for rs in sorted(exact):
        try:
            rows,pages=paged("associations",{"rs_id":rs,"sort":"p_value","direction":"asc"})
            seen=set();matches=[]
            for a in rows:
                alleles=a.get("snp_allele")
                require(isinstance(alleles,list) and any(v.get("rs_id")==rs for v in alleles if isinstance(v,dict)),"Exact association rsID unverified")
                aid=a["association_id"];require(aid not in seen,"Duplicate association ID");seen.add(aid)
                p=a.get("p_value");require(isinstance(p,(int,float)) and math.isfinite(p) and 0<=p<=1,"Catalog p-value")
                matches.append({k:a.get(k) for k in ("association_id","reported_trait","p_value","accession_id","pubmed_id","first_author","mapped_genes","snp_allele","multi_snp_haplotype","snp_interaction","beta_num","beta_unit","or_per_copy_num")})
            out["exact_associations"][rs]={"status":"COMPLETE_EXACT_CURATED_ASSOCIATIONS" if pages["complete"] else "HOLD_PAGE_CAP","pagination":pages,"rows":matches}
        except (OSError,ValueError,KeyError,TypeError,IndexError) as e:out["exact_associations"][rs]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:240]}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--proxy",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=audit(a.proxy)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:240]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    brief={**x,"regions":{k:{a:b for a,b in d.items() if a!="regional_variants"} for k,d in x.get("regions",{}).items()}}
    print(json.dumps(brief,indent=2,allow_nan=False))
if __name__=="__main__":main()
