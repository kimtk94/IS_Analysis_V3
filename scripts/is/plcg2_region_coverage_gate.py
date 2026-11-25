#!/usr/bin/env python3
"""Post-readout PLCG2 region coverage gate; not formal coloc or fine mapping."""
import argparse,collections,csv,importlib,json,math,re,time
from pathlib import Path
ext=importlib.import_module("scripts.is.fixed24_two_locus_qtl_extension")
gwas=importlib.import_module("scripts.is.fixed_24_subtype_extension")
http=ext.http
GENE="ENSG00000197943"
CH="16"
START,END=81275211,82275211
TARGET="chr16_81775211_C_T"
LEAD="chr16_81878856_T_G"
STUDIES=("EAS_CES","EAS_SVS","EUR_CES","EUR_SVS")
EXPECTED=11248
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read_interval(url,t,h,total,start,end,chr_field,pos_field):
    span=http.indexed_span(t,CH,start,end,total)
    require(span is not None,"No indexed region chunk")
    a,b,x,y=span
    require(b-a+1<=4000000,"Compressed interval exceeds fixed 4 MB cap")
    time.sleep(5)
    raw,cr=http.get(url,b-a+1,(a,b))
    require(int(cr.rsplit("/",1)[1])==total,"Source total changed during read")
    stream=ext.base.parse_bgzf_range(raw,x,y)["data_lines"]
    rows=[];decoded=0
    for line in csv.reader(stream,delimiter="\t"):
        decoded+=1;require(decoded<=500000,"Decoded interval row cap")
        require(len(line)==len(h),"Unexpected region row width")
        d=dict(zip(h,line))
        if d[chr_field].removeprefix("chr")==CH and start<=int(d[pos_field])<=end:
            rows.append(d)
    return rows,{"region":[start,end],"range":cr,"sha256":http.sha(raw),"decoded_rows":decoded,"in_region_rows":len(rows),"compressed_bytes":len(raw)}
def pair(ref,alt):
    return len(ref)==len(alt)==1 and ref in "ACGT" and alt in "ACGT" and ref!=alt and set(ref+alt) not in (set("AT"),set("CG"))
def read_qtl(root):
    prior=json.loads((root/"is_fixed24_two_locus_qtl_extension_20261125.json").read_text())["result"]["datasets"]["QTD000434"]
    p=json.loads((root/"is_plcg2_rosmap_egene_cs_20261125.json").read_text())["result"]["sources"]["permutation"]
    require(p["status"]=="COMPLETE_OFFICIAL_PERMUTATION_READ" and int(p["PLCG2"]["n_variants"])==EXPECTED and p["PLCG2"]["variant"]==LEAD,"Pinned permutation source")
    url=prior["source"]
    require(url=="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz","Exact QTL source")
    t,h,total,info=ext.load_index(url)
    require(all(info[k]==prior[k] for k in ("index_sha256","header_range_sha256","source_total_bytes")),"QTL source changed from verified single-position read")
    records={};rejected=set();c=collections.Counter();ranges=[]
    for start in range(START,END+1,250000):
        end=min(start+249999,END)
        rows,meta=read_interval(url,t,h,total,start,end,"chromosome","position")
        ranges.append(meta)
        for d in rows:
            if d["gene_id"]!=GENE and d["molecular_trait_id"]!=GENE:continue
            require(d["gene_id"]==d["molecular_trait_id"]==GENE,"Target gene mismatch")
            ref,alt=d["ref"],d["alt"];pos=int(d["position"])
            v=d["variant"]
            require(v==f"chr16_{pos}_{ref}_{alt}","Target variant key drift")
            c["target_rows"]+=1
            try:
                beta,se,pval,maf=map(float,[d["beta"],d["se"],d["pvalue"],d["maf"]])
                ac,an=map(int,[d["ac"],d["an"]])
                require(all(math.isfinite(x) for x in (beta,se,pval,maf)) and se>0 and 0<=pval<=1 and 0<maf<=.5 and 0<ac<an,"QTL numeric eligibility")
                item={"pos":pos,"ref":ref,"alt":alt,"beta_ALT":beta,"se":se,"p":pval,"alt_eaf":ac/an,"an":an,"maf":maf,"rsids":sorted(set(re.findall(r"rs[0-9]+",d["rsid"]))),"type":d["type"]}
            except (ValueError,TypeError,KeyError):
                rejected.add(v);c["invalid_rows"]+=1;continue
            if v in records:
                a=records[v];aa={k:x for k,x in a.items() if k!="rsids"};bb={k:x for k,x in item.items() if k!="rsids"}
                if aa!=bb:rejected.add(v);c["conflicting_duplicate"]+=1
                else:a["rsids"]=sorted(set(a["rsids"]+item["rsids"]));c["identical_alias_rows"]+=1
            else:records[v]=item
    for v in rejected:records.pop(v,None)
    require(TARGET in records and LEAD in records,"Missing fixed target or official eGene lead")
    a=records[TARGET];l=records[LEAD]
    require("rs7203740" in a["rsids"] and math.isclose(a["beta_ALT"],.279006,rel_tol=1e-7) and math.isclose(a["se"],.068636,rel_tol=1e-7) and math.isclose(a["p"],.0000551034,rel_tol=1e-7),"Target nominal parity")
    require(math.isclose(l["p"],7.86674e-11,rel_tol=1e-6) and math.isclose(l["beta_ALT"],-.246995,rel_tol=1e-6),"Official eGene lead parity")
    eligible={v:d for v,d in records.items() if d["type"]=="SNP" and pair(d["ref"],d["alt"])}
    return {"status":"REGION_QTL_READ","source":url,**info,"ranges":ranges,"counts":dict(c),"unique_valid_variants_in_window":len(records),"official_cis_tested_variants":EXPECTED,"strict_nonpal_snp_count":len(eligible),"target":a,"official_lead":l,"records":eligible}
def read_gwas(study,oldroot):
    old=json.loads((oldroot/study/(study+".json")).read_text())
    require(old["study"]==study and old["status"]=="INDEXED_FIXED_PANEL_READBACK","Prior GWAS artifact identity")
    src,info=gwas.open_source(study)
    for k in ("source","metadata_sha256","index_sha256","header_range_sha256","source_total_bytes"):
        require(info[k]==old[k],"GWAS source changed: "+k)
    url,t,total,h,f=src
    rows,rangeinfo=read_interval(url,t,h,total,START,END,f["chr"],f["pos"])
    groups=collections.defaultdict(list);c=collections.Counter()
    for d in rows:
        ref,alt=d[f["oa"]].upper(),d[f["ea"]].upper()
        if not pair(ref,alt):c["excluded_non_strict_allele"]+=1;continue
        k=(int(d[f["pos"]]),"".join(sorted(ref+alt)))
        groups[k].append(d)
    valid={}
    for k,group in groups.items():
        if len(group)!=1:c["duplicate_pair_keys"]+=1;continue
        d=group[0];rs=d[f["rsid"]].lower().strip()
        try:
            require(re.fullmatch(r"rs[0-9]+",rs) is not None,"No exact rsID")
            require(int(d[f["code"]]) in {10,11,12,13},"GWAS harmonisation code")
            beta,se,p,eaf=[float(d[f[x]]) for x in ("beta","se","p","eaf")]
            require(all(math.isfinite(x) for x in (beta,se,p,eaf)) and se>0 and 0<=p<=1 and 0<eaf<1,"Invalid GWAS statistics")
            valid[k]={"pos":k[0],"ea":d[f["ea"]].upper(),"oa":d[f["oa"]].upper(),"rsid":rs,"beta":beta,"se":se,"p":p,"eaf":eaf}
        except (ValueError,KeyError,TypeError):c["invalid_rows"]+=1
    prior=[x for x in old["rows"].values() if x.get("rsid")=="rs7203740" and x.get("status")=="MATCHED_HARMONISED_ROW"]
    anchor=valid.get((81775211,"CT"))
    if prior:
        require(len(prior)==1 and anchor is not None and anchor["rsid"]=="rs7203740","Prior GWAS target missing")
        sign=1 if anchor["ea"]=="T" else -1
        a=prior[0]
        require(math.isclose(anchor["beta"]*sign,a["beta_alt"],abs_tol=1e-10) and math.isclose(anchor["se"],a["se"],abs_tol=1e-10) and math.isclose(anchor["p"],a["p"],rel_tol=1e-10,abs_tol=1e-15),"GWAS target parity")
    return {"status":"REGION_GWAS_READ","source":url,"source_provenance":{k:info[k] for k in ("metadata_sha256","index_sha256","header_range_sha256","source_total_bytes","publisher_full_file_md5")},"range":rangeinfo,"counts":dict(c),"unique_strict_keys":len(valid),"prior_target_matched":bool(prior),"records":valid}
def compare(q,g):
    matched=[];c=collections.Counter()
    for v,d in q["records"].items():
        k=(d["pos"],"".join(sorted(d["ref"]+d["alt"])))
        a=g["records"].get(k)
        if not a:c["no_unique_eligible_pair"]+=1;continue
        if a["rsid"] not in d["rsids"]:c["rsid_not_corroborated"]+=1;continue
        sign=1 if a["ea"]==d["alt"] else -1
        matched.append({"variant":v,"rsid":a["rsid"],"pos":d["pos"],"ref":d["ref"],"alt":d["alt"],"qtl_beta_ALT":d["beta_ALT"],"qtl_se":d["se"],"qtl_p":d["p"],"qtl_alt_eaf":d["alt_eaf"],"gwas_beta_ALT":sign*a["beta"],"gwas_se":a["se"],"gwas_p":a["p"],"gwas_alt_eaf":a["eaf"] if sign==1 else 1-a["eaf"]})
    n=len(matched)
    return {"status":"DESCRIPTIVE_ALLELE_MATCHED_REGION","matched":n,"qtl_denominator":q["strict_nonpal_snp_count"],"qtl_coverage":n/q["strict_nonpal_snp_count"],"exclusion_counts":dict(c),"target":[d for d in matched if d["variant"]==TARGET],"official_lead":[d for d in matched if d["variant"]==LEAD],"matched_gwas_p_lt_5e_8":sum(d["gwas_p"]<5e-8 for d in matched),"top_matched_gwas":min(matched,key=lambda d:d["gwas_p"]) if matched else None,"top_matched_qtl":min(matched,key=lambda d:d["qtl_p"]) if matched else None,"matched_rows":matched}
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--prior",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x={"status":"POST_READOUT_REGION_COVERAGE_GATE","GRCh38_region":f"16:{START}-{END}","gene":GENE,"target":TARGET,"official_eGene_lead":LEAD,"comparisons":{},"limits":["Fixed index +/-500 kb post-readout descriptive window; not a certified complete cis universe. All strict nonpalindromic SNPs retained, not P-selected.","Exact coordinate and allele pair plus corroborating rsID; no proxy or ambiguous strand substitution. Region rows and index/header bytes only; no independent full-file checksum or FASTA certificate.","ROSMAP molecular ancestry not matched to EAS; nominal N/AN not model N or sdY. No approved cohort LD, effect-scale equivalence, formal coloc, MR, fine mapping or independent replication. Coverage alone is not a pass for formal coloc."]}
    try:
        q=read_qtl(a.evidence);x["qtl"]=q
    except (OSError,ValueError,TypeError,KeyError,IndexError) as e:
        x["qtl"]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:300]}
    if x["qtl"]["status"]=="REGION_QTL_READ":
        for study in STUDIES:
            try:
                g=read_gwas(study,a.prior)
                c=compare(q,g);c["gwas"]={k:v for k,v in g.items() if k!="records"}
                x["comparisons"][study]=c
            except (OSError,ValueError,TypeError,KeyError,IndexError) as e:x["comparisons"][study]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:300]}
    else:
        for study in STUDIES:x["comparisons"][study]={"status":"HOLD_QTL_SOURCE"}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    brief={**x,"qtl":{k:v for k,v in x["qtl"].items() if k!="records"},"comparisons":{k:{a:b for a,b in v.items() if a!="matched_rows"} for k,v in x["comparisons"].items()}}
    print(json.dumps(brief,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
