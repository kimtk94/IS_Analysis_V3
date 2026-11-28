#!/usr/bin/env python3
"""Exact published Chinese skin summary tables: target coverage, not automatic causal inference."""
import argparse,collections,csv,hashlib,io,json,math,re
from pathlib import Path
from urllib.parse import quote
import skin_chinese_qtl_sources as prior
base=prior.base;require=base.require;ROOT=prior.ROOT
ARCH=ROOT/"docs/publication/evidence/skin_chinese_first_readback_20261125.json"
EXPECTED={"eQTL":("skin eQTL summary statistics.txt",70906052,"caca80b5a56f3fd9ceeb3b03b9b4cc07"),"pQTL":("skin pQTL summary statistics.txt",8085897,"5ed66ce34705c8aeff36b68162c73b64")}
GENE=prior.GENE;LEAD=prior.LEAD;LOW=prior.LOW;HIGH=prior.HIGH
def archive():
    a=json.loads(ARCH.read_text());require(a["run_id"]==45614087522 and len(a["records"])==3,"prior run identity")
    d={q["aggregate"]["mode"]:q["aggregate"] for q in a["records"]}
    require(d["eQTL"]["rows"]==24788 and d["eQTL"]["molecular_traits"]==282 and d["pQTL"]["rows"]==1765,"prior selected table denominator")
    f={v["name"]:(v["bytes"],v["checksum"]) for v in d["inventory"]["files"]}
    require(f=={name:(n,"md5:"+md5) for name,n,md5 in EXPECTED.values()},"public metadata mismatch")
    return d
def field(h,options):
    hits=[v for v in options if v in h];require(len(hits)==1,"missing/ambiguous column: "+"/".join(options));return hits[0]
def header(raw,mode):
    # Recognize explicit source header fields, not positional assumptions.
    choices=[]
    for delim in ("\t",","):
        rr=csv.reader(io.StringIO(raw.decode("utf-8-sig")),delimiter=delim)
        for i,line in enumerate(rr):
            if i>=8:break
            h=[v.strip() for v in line]
            if len(h)!=len(set(h)):continue
            try:
                cols={"variant":field(h,("variant_id","variant")),"trait":field(h,("gene_id",) if mode=="eQTL" else ("protein_id","pro_id","gene_id")),
                      "beta":field(h,("slope","beta")),"se":field(h,("slope_se","se","standard_error")),"p":field(h,("pval_nominal","pvalue","p_value")),"maf":field(h,("maf","MAF"))}
            except ValueError:continue
            choices.append((delim,i,h,cols))
    require(len(choices)==1,"unique explicit tab/comma summary header not established")
    return choices[0]
def scan(mode):
    archive();name,n,md5=EXPECTED[mode]
    url="https://zenodo.org/records/19638349/files/"+quote(name)+"?download=1"
    raw,src=base.get(url,n)
    require(len(raw)==n and hashlib.md5(raw).hexdigest()==md5,"public summary size/MD5 mismatch")
    out={"status":"CHINESE_SKIN_SUMMARY_SCHEMA","source":src,"md5":md5,"assay":mode}
    try:delim,offset,h,cols=header(raw,mode)
    except ValueError as e:
        # Header-like lines only; never export arbitrary data rows.
        first=raw.decode("utf-8-sig").splitlines()[:8]
        terms={"gene_id","protein_id","pro_id","gene_name","variant_id","variant","slope","slope_se","pval_nominal","maf","beta","se","pvalue","p_value","standard_error","MAF","ref","alt","chr_var","pos_var"}
        out.update(status="HOLD_SUMMARY_HEADER",reason=str(e),known_header_terms=sorted({v.strip() for l in first for v in re.split(r"[\t, ]+",l) if v.strip() in terms}))
        return out
    text=raw.decode("utf-8-sig");stream=io.StringIO(text)
    for i in range(offset):next(stream)
    rr=csv.DictReader(stream,delimiter=delim)
    require(rr.fieldnames==h,"trimmed header not exact")
    traits=set();nrows=0;valid=0;bad=0;target_count=0;reg_count=0;target=[];regional=[];exact=[];duplicates=0;seen=set()
    # The official pQTL file labels its trait column gene_id; establish its semantics
    # from the exact pinned significant protein table before interpreting the field.
    known_proteins={}
    if mode=="pQTL":
        fn,sha,_,_=prior.FILES["pQTL"]
        aux,auxsrc=base.get(f"https://raw.githubusercontent.com/zzwang1030/skin_pQTL/{prior.REF}/data/{fn}",1_000_000)
        require(prior.git_blob(aux)==sha,"author protein mapping drift")
        m=collections.defaultdict(set)
        for row in csv.DictReader(io.StringIO(aux.decode("utf-8-sig"))):m[row["protein_id"]].add(row["gene_name"])
        require(len(m)==36,"author protein denominator")
        known_proteins={k:next(iter(v)) for k,v in m.items() if len(v)==1}
        out["protein_mapping_source"]=auxsrc
    prior2=json.loads((ROOT/"docs/publication/evidence/skin_chinese_second_readback_20261125.json").read_text())
    require(prior2["run_id"]==45614624939 and len(prior2["records"])==2,"second readback identity")
    old=next(v["aggregate"] for v in prior2["records"] if v["aggregate"]["mode"]==mode)
    require(src["sha256"]==old["source"]["sha256"],"prior exact source SHA drift")
    for row in rr:
        if not any(row.values()):continue
        nrows+=1
        require(None not in row and all(v is not None for v in row.values()),"summary row width")
        trait=row[cols["trait"]].split(".")[0]
        variant=row[cols["variant"]]
        require(trait and variant,"blank keys")
        traits.add(trait)
        m=re.fullmatch(r"(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y):(\d+):([ACGT]+):([ACGT]+)",variant)
        if not m:bad+=1;continue
        ch,p,a,b=m.groups();pos=int(p);vid=f"{ch}:{p}:{a}:{b}"
        g=row.get("gene_name") or known_proteins.get(trait)
        is_target=trait==GENE or g=="GCC2"
        is_region=ch=="2" and LOW<=pos<=HIGH
        if is_target:target_count+=1
        if is_region:reg_count+=1
        try:
            beta,se,pv,maf=[float(row[cols[k]]) for k in ("beta","se","p","maf")]
            require(all(math.isfinite(v) for v in (beta,se,pv,maf)) and se>0 and 0<=pv<=1 and 0<maf<=.5,"invalid molecular statistic")
            if {"ref","alt"}<=set(h):require(row["ref"]==a and row["alt"]==b,"source variant allele disagreement")
            if {"chr_var","pos_var"}<=set(h):require(row["chr_var"].removeprefix("chr")==ch and int(row["pos_var"])==pos,"source coordinate disagreement")
        except (ValueError,KeyError,TypeError):bad+=1;continue
        valid+=1
        if not (is_target or is_region):continue
        k=(trait,vid)
        if k in seen:duplicates+=1
        seen.add(k)
        z={"trait":trait,"variant":vid,"reported_slope":beta,"se":se,"p":pv,"maf":maf}
        if g and re.fullmatch(r"[A-Za-z0-9_.-]{1,40}",g):z["gene_name"]=g
        if is_target:target.append(z)
        if is_region:regional.append(z)
        if vid==LEAD:exact.append(z)
    require(nrows>=1000,"implausibly small source")
    if mode=="pQTL":
        out["protein_trait_identity"]={"header":cols["trait"],"exact_author_protein_keys":len(traits&set(known_proteins)),"total_source_traits":len(traits),"all_traits_exact_author_protein_keys":traits<=set(known_proteins)}
        if not traits<=set(known_proteins):
            out.update(status="HOLD_PROTEIN_TRAIT_IDENTITY",rows=nrows,molecular_traits=len(traits),valid_statistic_rows=valid,invalid_rows=bad,reason="Not all source traits have exact author protein keys; no inferred protein names")
            return out
    out.update(status="CHINESE_SKIN_PUBLISHED_SUMMARY_READBACK",header=h,delimiter=repr(delim),rows=nrows,molecular_traits=len(traits),valid_statistic_rows=valid,invalid_rows=bad,
               GCC2_rows=target_count,GCC2_valid_rows=len(target),GCC2_region_rows=sum(v["variant"].startswith("2:") and LOW<=int(v["variant"].split(":")[1])<=HIGH for v in target),
               region_rows=reg_count,regional_valid_rows=len(regional),regional_traits=len({v["trait"] for v in regional}),duplicate_target_or_region_keys=duplicates,
               exact_rs4584990_rows=exact,GCC2_best_nominal=sorted(target,key=lambda z:z["p"])[:10],GCC2_summary_rows=target)
    out["limits"]=["Exact published Zenodo file MD5 verified, but may be selected/filtered summary statistics rather than all cis pairs; full assay universe not established by file size. No negative association claim from absent rows. Report slopes without assigning effect allele until methods confirmed. No causal inference or LD-aware regional colocalization. Chinese donors do not constitute independent Korean wrinkle outcome replication."]
    return out
def selftest():
    archive()
    raw=b"gene_id\tvariant_id\tslope\tslope_se\tpval_nominal\tmaf\nENSG1\t2:108872562:C:T\t.1\t.02\t.1\t.2\n"
    d,i,h,c=header(raw,"eQTL");assert d=="\t" and c["se"]=="slope_se"
    assert header(raw,"pQTL")[3]["trait"]=="gene_id"  # source label alone does not certify protein identity
    print("Published file identity, prior denominators, and explicit summary schema checks passed.")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=["eQTL","pQTL","selftest"]);ap.add_argument("--out",type=Path);a=ap.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=scan(a.mode)
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    z.update(mode=a.mode,scientific_status="HOLD_FULL_REGIONAL_AND_ANCESTRY_VALIDATION")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps({k:v for k,v in z.items() if k!="GCC2_summary_rows"},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
