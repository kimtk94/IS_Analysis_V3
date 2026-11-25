#!/usr/bin/env python3
"""Fixed rs2767060 LAS regional heterogeneity readback, descriptive only."""
import argparse,collections,json,math,importlib
from pathlib import Path
base=importlib.import_module("scripts.is.plcg2_region_coverage_gate")
cat=importlib.import_module("scripts.is.fixed24_two_locus_catalog_audit")
START,END=35333363,36333363
POS=35833363
STUDIES=("EAS_LAS","EUR_LAS")
def req(ok,msg):
    if not ok:raise ValueError(msg)
def read(study,prior):
    old=json.loads((prior/study/(study+".json")).read_text())
    req(old["study"]==study and old["status"]=="INDEXED_FIXED_PANEL_READBACK","Prior GWAS identity")
    src,info=base.gwas.open_source(study)
    for k in ("source","metadata_sha256","index_sha256","header_range_sha256","source_total_bytes"):
        req(info[k]==old[k],"Prior GWAS source changed: "+k)
    url,t,total,h,f=src
    # Reuse validated interval reader for chr10, not chr16.
    base.CH="10"
    rows,meta=base.read_interval(url,t,h,total,START,END,f["chr"],f["pos"])
    groups=collections.defaultdict(list);c=collections.Counter()
    for d in rows:
        ea,oa=d[f["ea"]].upper(),d[f["oa"]].upper()
        if not base.pair(oa,ea):c["excluded_non_strict_allele"]+=1;continue
        groups[(int(d[f["pos"]]),"".join(sorted(ea+oa)))].append(d)
    records={}
    for k,rr in groups.items():
        if len(rr)!=1:c["duplicate_pair_keys"]+=1;continue
        z=rr[0]
        try:
            import re
            rs=z[f["rsid"]].strip().lower()
            req(re.fullmatch(r"rs[0-9]+",rs) is not None and int(z[f["code"]]) in {10,11,12,13},"SNP identity")
            beta,se,p,eaf=map(float,[z[f[x]] for x in ("beta","se","p","eaf")])
            req(all(math.isfinite(v) for v in (beta,se,p,eaf)) and se>0 and 0<=p<=1 and 0<eaf<1,"Numerical fields")
            records[k]={"rsid":rs,"pos":k[0],"ea":z[f["ea"]].upper(),"oa":z[f["oa"]].upper(),"beta":beta,"se":se,"p":p,"eaf":eaf}
        except (ValueError,KeyError,TypeError):c["invalid_rows"]+=1
    anchor=records.get((POS,"CT"));p=[x for x in old["rows"].values() if x.get("rsid")=="rs2767060" and x["status"]=="MATCHED_HARMONISED_ROW"]
    req(anchor and anchor["rsid"]=="rs2767060" and len(p)==1,"Anchor missing")
    sign=1 if anchor["ea"]=="T" else -1
    req(sign*anchor["beta"]==p[0]["beta_alt"] and anchor["se"]==p[0]["se"] and anchor["p"]==p[0]["p"],"Prior anchor drift")
    return records,{"source":info,"range":meta,"strict_pairs":len(records),"exclusions":dict(c),"anchor":anchor}
def rank(rows,field,anchor):
    v=anchor[field]
    return [1+sum(d[field]<v for d in rows),sum(d[field]<=v for d in rows)]
def summary(rows):
    if not rows:return {"n":0}
    top=lambda key:min(rows,key=lambda d:(d[key],d["rsid"]))
    return {"n":len(rows),"EAS_p_lt_5e_8":sum(d["EAS_p"]<5e-8 for d in rows),"EUR_p_lt_5e_8":sum(d["EUR_p"]<5e-8 for d in rows),"same_sign_count":sum(d["EAS_beta"]*d["EUR_beta"]>0 for d in rows),"zero_count":sum(d["EAS_beta"]*d["EUR_beta"]==0 for d in rows),"top_EAS":top("EAS_p"),"top_EUR":top("EUR_p"),"top_conditional_max_variance":top("p_max_variance"),"nominal_max_variance_p_lt_0_05":sum(d["p_max_variance"]<.05 for d in rows)}
def run(prior,proxyfile):
    proxy=cat.read_proxy(proxyfile)
    a,am=read("EAS_LAS",prior);b,bm=read("EUR_LAS",prior)
    rows=[];rsid_conflict=0
    for k in sorted(set(a)&set(b)):
        u,v=a[k],b[k]
        if u["rsid"]!=v["rsid"]:rsid_conflict+=1;continue
        # Alphabetically second nucleotide as common effect allele.
        ea=k[1][1];su=1 if u["ea"]==ea else -1;sv=1 if v["ea"]==ea else -1
        bu,bv=su*u["beta"],sv*v["beta"];delta=bu-bv
        se=u["se"]+v["se"]
        row={"rsid":u["rsid"],"pos":k[0],"effect_allele":ea,"other_allele":k[1][0],"EAS_beta":bu,"EAS_se":u["se"],"EAS_p":u["p"],"EAS_eaf":u["eaf"] if su==1 else 1-u["eaf"],"EUR_beta":bv,"EUR_se":v["se"],"EUR_p":v["p"],"EUR_eaf":v["eaf"] if sv==1 else 1-v["eaf"],"delta":delta,"p_max_variance":math.erfc(abs(delta)/se/math.sqrt(2)),"ci95_max_variance":[delta-1.959964*se,delta+1.959964*se],"reference_proxy_populations":proxy.get(u["rsid"],{})}
        rows.append(row)
    anchor=[d for d in rows if d["pos"]==POS and d["rsid"]=="rs2767060" and d["effect_allele"]=="T" and d["other_allele"]=="C"]
    req(len(anchor)==1 and math.isclose(anchor[0]["delta"],.3672,abs_tol=1e-10),"Original cross-ancestry contrast parity")
    anchor=anchor[0]
    proxy_rows=[d for d in rows if d["rsid"] in proxy];proxy_rs={d["rsid"] for d in proxy_rows}
    # Preserve all matched rows; no selection of new candidates by association P.
    return {"status":"DESCRIPTIVE_POST_READOUT_LAS_REGION","region":f"10:{START}-{END}","sources":{"EAS_LAS":am,"EUR_LAS":bm},"exact_matched":len(rows),"rsid_conflicts":rsid_conflict,"EAS_matched_fraction":len(rows)/len(a),"EUR_matched_fraction":len(rows)/len(b),"summary_all":summary(rows),"reference_proxy_union":78,"matched_proxy_rsids":len(proxy_rs),"missing_proxy_rsids":sorted(set(proxy)-proxy_rs),"summary_reference_proxies":summary(proxy_rows),"anchor":anchor,"anchor_ranks":{k:rank(rows,k,anchor) for k in ("EAS_p","EUR_p","p_max_variance")},"rows":rows,"limits":["Post-readout fixed 1 Mb region, same GIGASTROKE release; not independent confirmation or a new multiplicity family for significance claims.","No effect-based filtering; conditional max variance difference assumes comparable reported scales. Nominal counts and rank only; no newly declared significant locus.","Reference proxy membership does not certify allele-specific LD phase or GWAS-matched LD; sign across different variants not interpreted as phased support. No causal gene, fine mapping, MR or coloc."]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--prior",type=Path,required=True);p.add_argument("--proxy",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.prior,a.proxy)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:300]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in x.items() if k!="rows"},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
