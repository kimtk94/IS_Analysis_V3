#!/usr/bin/env python3
"""Fixed GCC2 regional coverage and association-profile audit. Not a colocalization test."""
import argparse,collections,csv,gzip,hashlib,io,json,math,re
from pathlib import Path
import early_skin_molecular_context as point
base=point.base;idx=point.idx;require=base.require;ROOT=point.ROOT
ARCH=ROOT/"docs/publication/evidence/early_track_qtl_and_node_readback_20261125.json"
GENE="ENSG00000135968"
LEAD="2_108872562_C_T"
LOW=108372562;HIGH=109372562
def archive():
    d=json.loads(ARCH.read_text());require(len(d["records"])==10,"expected ten verified prior jobs")
    q=[z["aggregate"] for z in d["records"] if z["run_id"]==45612527929 and z["name"]=="combine"]
    require(len(q)==1 and q[0]["completed_queries"]==48 and q[0]["observed_nominal_tests"]==665,"prior fixed family")
    q=q[0]["strong_nominal_rows"]
    require(len(q)==1 and q[0]["gene"]==GENE and q[0]["variant"]==LEAD,"frozen GCC2 candidate")
    return q[0]
def key(pos,a,b):
    require(a!=b and re.fullmatch("[ACGT]",a) and re.fullmatch("[ACGT]",b) and {a,b} not in ({"A","T"},{"C","G"}),"nonpalindromic SNV required")
    return "2_"+str(pos)+"_"+"_".join(sorted((a,b)))
def unique(pool):
    # Exclude duplicate records and all positions with multiple allele pairs, before comparison.
    loc=collections.Counter(k.split("_")[1] for k in pool)
    return {k:v[0] for k,v in pool.items() if len(v)==1 and loc[k.split("_")[1]]==1}
def gwas():
    prior=point.leads()[LEAD]["traits"]["GCST90094903"]
    src0=json.loads(point.READBACK.read_text())["records"]
    d=next(v["aggregate"]["full_source_readout"] for v in src0 if v["aggregate"].get("accession")=="GCST90094903")
    f=d["full_source"];raw,src=base.get(f["url"],f["bytes"])
    require(len(raw)==f["bytes"] and src["sha256"]==f["sha256"] and hashlib.md5(raw).hexdigest()==f["md5"],"verified GWAS whole-file digest drift")
    pool=collections.defaultdict(list);n=0;region=0;bad=0
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
        rr=csv.DictReader(io.TextIOWrapper(stream,encoding="utf8"),delimiter="\t")
        fields={"hm_chrom","hm_pos","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value","hm_code","hm_effect_allele_frequency","hm_rsid"}
        require(fields<=set(rr.fieldnames or []) and len(rr.fieldnames)==len(set(rr.fieldnames)),"exact GWAS header")
        for row in rr:
            n+=1
            if row["hm_chrom"].removeprefix("chr")!="2":continue
            try:p=int(row["hm_pos"])
            except ValueError:continue
            if not LOW<=p<=HIGH:continue
            region+=1
            try:
                a,b=row["hm_effect_allele"],row["hm_other_allele"];k=key(p,a,b)
                require(int(row["hm_code"]) in (10,11,12,13),"unresolved harmonisation")
                beta,se,pv,eaf=[float(row[z]) for z in ("hm_beta","standard_error","p_value","hm_effect_allele_frequency")]
                require(all(math.isfinite(x) for x in (beta,se,pv,eaf)) and se>0 and 0<=pv<=1 and 0<eaf<1,"invalid statistic")
                s=1 if a==k.split("_")[-1] else -1
                pool[k].append({"beta_second":s*beta,"se":se,"p":pv,"second_eaf":eaf if s==1 else 1-eaf,"rsid":row["hm_rsid"]})
            except (ValueError,TypeError,KeyError):bad+=1
    require(n==4950621,"full GWAS row count")
    rows=unique(pool);z=rows.get(LEAD);require(z is not None and z["rsid"]=="rs4584990" and all(math.isclose(z[k],prior[v],rel_tol=1e-9,abs_tol=1e-12) for k,v in (("beta_second","beta_alt"),("se","standard_error"),("p","p_value"))),"frozen GWAS anchor mismatch")
    return {"status":"VERIFIED_GWAS_REGION","study":"GCST90094903","region":["2",LOW,HIGH],"source":src,"full_md5":f["md5"],"whole_rows":n,"regional_rows":region,"eligible_rows":len(rows),"rejected_rows":bad,"duplicate_or_multiallelic_keys":len(pool)-len(rows),"rows":rows}
def qtl(mode):
    import time
    t=point.TISSUES[mode];url=point.FILEBASE+t+".tsv.gz";prior=archive()["tissues"][mode]
    ir,_=idx.get(url+".tbi",12_000_000);tbi=idx.parse_index(ir)
    first,cr=idx.get(url,65536,(0,65535));total=int(cr.rsplit("/",1)[1])
    h=idx.parse_bgzf_range(first,0,0)["data_lines"].readline().lstrip("#").rstrip("\r\n").split("\t")
    point.validate_header(h,tbi)
    pool=collections.defaultdict(list);ranges=[];cache={};scanned=0;rejected=0;empty=0
    # Nonoverlapping query windows, split oversized index spans; no p-value-based selection.
    windows=[(s,min(s+49999,HIGH)) for s in range(LOW,HIGH+1,50000)]
    while windows:
        lo,hi=windows.pop(0)
        try:span=idx.indexed_span(tbi,"2",lo,hi,total)
        except ValueError as e:
            require("range exceeds budget" in str(e).lower() and hi-lo>=10000,"index range contract: "+str(e)[:100])
            mid=(lo+hi)//2;windows[:0]=[(lo,mid),(mid+1,hi)];continue
        if span is None:empty+=1;continue
        a,b,x,y=span
        if span not in cache:
            require(sum(len(v[0]) for v in cache.values())+(b-a+1)<=240_000_000,"regional range cap")
            time.sleep(.7);rr,rc=idx.get(url,b-a+1,(a,b));cache[span]=(rr,rc)
        rr,rc=cache[span];ranges.append({"window":[lo,hi],"range":rc,"sha256":idx.sha(rr)})
        for row in csv.reader(idx.parse_bgzf_range(rr,x,y)["data_lines"],delimiter="\t"):
            if len(row)!=len(h):continue
            d=dict(zip(h,row))
            if d["chromosome"].removeprefix("chr")!="2" or d["gene_id"].split(".")[0]!=GENE:continue
            try:p=int(d["position"])
            except ValueError:continue
            if not lo<=p<=hi:continue
            scanned+=1
            try:
                ref,alt=d["ref"],d["alt"];k=key(p,ref,alt)
                v=point.exact_row(d,"2_"+str(p)+"_"+ref+"_"+alt)
                require(v["gene"].split(".")[0]==GENE,"wrong target gene")
                s=1 if alt==k.split("_")[-1] else -1
                pool[k].append({"beta_second":s*v["beta_alt"],"se":v["se"],"p":v["p"],"maf":v["maf"],"allele_number":v["allele_number"]})
            except (ValueError,TypeError,KeyError):rejected+=1
    rows=unique(pool);z=rows.get(LEAD)
    require(z is not None and all(math.isclose(z[k],prior[v],rel_tol=1e-9,abs_tol=1e-12) for k,v in (("beta_second","beta_alt"),("se","se"),("p","p"))),"frozen QTL anchor mismatch")
    return {"status":"INDEXED_GCC2_REGION","tissue":t,"gene":GENE,"region":["2",LOW,HIGH],"source_url":url,"index_sha256":idx.sha(ir),"header_sha256":idx.sha(first),"remote_bytes":total,"query_ranges":ranges,"empty_index_windows":empty,"regional_gene_rows":scanned,"eligible_rows":len(rows),"rejected_rows":rejected,"duplicate_or_multiallelic_keys":len(pool)-len(rows),"rows":rows}
def combine(root):
    from scipy.stats import spearmanr
    a=json.loads((root/"gwas.json").read_text())
    require(a.get("status")=="VERIFIED_GWAS_REGION" and a.get("region")==["2",LOW,HIGH],"GWAS regional contract")
    g=a["rows"];result={}
    for mode,t in point.TISSUES.items():
        f=root/(mode+".json");b=json.loads(f.read_text()) if f.exists() else {}
        if b.get("status")!="INDEXED_GCC2_REGION":result[mode]={"status":"HOLD_QTL_SOURCE","reason":b.get("reason")};continue
        require(b["region"]==a["region"] and b["tissue"]==t and b["gene"]==GENE,"QTL artifact identity")
        q=b["rows"];common=set(g)&set(q);n=len(common)
        gf=n/len(g) if g else 0;qf=n/len(q) if q else 0
        z={"status":"REGIONAL_COVERAGE_ONLY","GWAS_eligible":len(g),"QTL_eligible":len(q),"exact_common":n,"GWAS_fraction":gf,"QTL_fraction":qf,
           "coverage_gate":n>=100 and gf>=.8 and qf>=.8 and LEAD in common}
        if n>=100 and LEAD in common:
            keys=sorted(common);x=[-math.log10(max(g[k]["p"],1e-300)) for k in keys];y=[-math.log10(max(q[k]["p"],1e-300)) for k in keys]
            rho=float(spearmanr(x,y).statistic)
            top=lambda d:sorted(common,key=lambda k:(d[k]["p"],k))[:10]
            z.update(status="DESCRIPTIVE_REGIONAL_PROFILE",spearman_neg_log_p=rho if math.isfinite(rho) else None,top10_overlap=len(set(top(g))&set(top(q))),
                     GWAS_top_common=top(g)[0],QTL_top_common=top(q)[0],anchor_GWAS_rank=1+sum(g[k]["p"]<g[LEAD]["p"] for k in common),anchor_QTL_rank=1+sum(q[k]["p"]<q[LEAD]["p"] for k in common))
        result[mode]=z
    return {"status":"GCC2_REGIONAL_PROFILE_NOT_COLOCALIZATION","gene":GENE,"variant":"rs4584990","region":["2",LOW,HIGH],"results":result,
            "limits":["Exploratory post-point-readout 1 Mb window fixed before regional data, exact nonpalindromic unique SNVs, no missing-as-zero. Spearman/top10 overlap descriptive; correlated SNPs are not independent, no p-value. >=100 common and >=80% of both source universes are conservative coverage gate only. Korean GWAS and mainly European GTEx, multiple nearby credible sets, LD and phenotype scale unresolved. No coloc PP.H4, SuSiE or MR claim. Future ancestry-matched LD/conditional signal analysis is necessary."]}
def selftest():
    assert archive()["variant"]==LEAD
    assert key(108872562,"T","C")==LEAD
    assert unique({"2_10_A_C":[1],"2_10_A_G":[2],"2_11_C_T":[3]})=={"2_11_C_T":3}
    try:key(10,"A","T")
    except ValueError:pass
    else:raise AssertionError("palindromic accepted")
    print("Frozen regional candidate, exact allele keys, duplicate and multiallelic exclusion passed.")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=["gwas",*point.TISSUES,"combine","selftest"]);ap.add_argument("--out",type=Path);ap.add_argument("--root",type=Path,default=Path("inputs"));a=ap.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=gwas() if a.mode=="gwas" else qtl(a.mode) if a.mode in point.TISSUES else combine(a.root)
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    z.update(mode=a.mode,scientific_status="HOLD_ANCESTRY_MATCHED_CONDITIONAL_VALIDATION")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    # Full public regional rows in artifact only; stdout compact for traceable readback.
    print(json.dumps({k:v for k,v in z.items() if k not in {"rows","query_ranges"}},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
