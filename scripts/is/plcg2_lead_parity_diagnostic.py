#!/usr/bin/env python3
"""Bounded nominal-vs-permutation lead audit; no tolerance relaxation."""
import argparse,collections,importlib,json,math,re
from pathlib import Path
region=importlib.import_module("scripts.is.plcg2_region_coverage_gate")
ext=region.ext
GENE=region.GENE
POS=81878856
VAR=region.LEAD
def require(ok,msg):
    if not ok:raise ValueError(msg)
def run(root):
    p=json.loads((root/"is_plcg2_rosmap_egene_cs_20261125.json").read_text())["result"]["sources"]["permutation"]
    old=json.loads((root/"is_fixed24_two_locus_qtl_extension_20261125.json").read_text())["result"]["datasets"]["QTD000434"]
    require(p["status"]=="COMPLETE_OFFICIAL_PERMUTATION_READ" and p["PLCG2"]["variant"]==VAR,"Prior permutation identity")
    url=old["source"]
    require(url=="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz","Exact nominal URL")
    t,h,total,info=ext.load_index(url)
    require(all(info[k]==old[k] for k in ("index_sha256","header_range_sha256","source_total_bytes")),"Nominal index/header/size changed")
    rows,rangeinfo=region.read_interval(url,t,h,total,POS,POS,"chromosome","position")
    observed=[];variants=[]
    for d in rows:
        if d["gene_id"]!=GENE and d["molecular_trait_id"]!=GENE:continue
        require(d["gene_id"]==d["molecular_trait_id"]==GENE,"Target gene ID disagreement")
        require(d["variant"]==f"chr16_{POS}_{d['ref']}_{d['alt']}","Nominal variant key")
        vals={k:float(d[k]) for k in ("beta","se","pvalue","maf")}
        require(all(math.isfinite(v) for v in vals.values()) and vals["se"]>0 and 0<=vals["pvalue"]<=1 and 0<vals["maf"]<=.5,"Invalid nominal stats")
        item={**{k:d[k] for k in ("gene_id","molecular_trait_id","variant","rsid","type","ref","alt")},**vals,"ac":int(d["ac"]),"an":int(d["an"])}
        require(0<item["ac"]<item["an"],"Nominal allele counts")
        variants.append(item)
        if d["variant"]==VAR:observed.append(item)
    out={"status":"NOMINAL_PERMUTATION_LEAD_DIAGNOSTIC","source":url,"nominal_provenance":info,"range":rangeinfo,"permutation_lead":p["PLCG2"],"permutation_source_sha256":p["source"]["sha256"],"target_gene_position_rows":variants,"exact_lead_rows":len(observed),"parity":{},"limits":["Same QTL dataset, not independent evidence. Report raw nominal and permutation beta/p; no assumption that models, scaling, sample sets or precision are identical. Do not relax the previous gate or perform formal coloc from this diagnostic."]}
    if len(observed)==1:
        v=observed[0];b=float(p["PLCG2"]["beta"]);pv=float(p["PLCG2"]["pvalue"])
        out["parity"]={"beta_exact":v["beta"]==b,"p_exact":v["pvalue"]==pv,"beta_abs_difference":abs(v["beta"]-b),"p_abs_difference":abs(v["pvalue"]-pv),"beta_sign_equal":v["beta"]*b>0,"nominal_p_over_permutation_p":v["pvalue"]/pv if pv else None,"nominal_wald_z":v["beta"]/v["se"],"nominal_two_sided_wald_p":math.erfc(abs(v["beta"]/v["se"])/math.sqrt(2)),"original_gate_beta_pass":math.isclose(v["beta"],b,rel_tol=1e-6),"original_gate_p_pass":math.isclose(v["pvalue"],pv,rel_tol=1e-6)}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.evidence)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:300]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
