#!/usr/bin/env python3
"""Independent East Asian skin QTL source qualification; significant-only tables are not allpairs."""
import argparse,collections,csv,hashlib,io,json,math,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"prior/scripts/publication"))
import early_track_acceleration as base
require=base.require
base.HOSTS.update({"raw.githubusercontent.com","zenodo.org"})
REF="74b17724f23f6de983a0b76320f95f85439a289a"
FILES={"eQTL":("TableS5_skin_eQTL_allSig_pairs.csv","f4e83aebb5b39d46948387294443d8eb2ad35b05",20999,282),
       "pQTL":("TableS4_skin_pQTL_allSig_pairs.csv","1b8080dad7ca60e6ad66b93d4b7a131b2c560e7c",1765,36)}
LEAD="2:108872562:C:T";GENE="ENSG00000135968";LOW=108372562;HIGH=109372562
def git_blob(raw):return hashlib.sha1(("blob %d\0"%len(raw)).encode()+raw).hexdigest()
def table(mode):
    name,sha,n0,g0=FILES[mode];url=f"https://raw.githubusercontent.com/zzwang1030/skin_pQTL/{REF}/data/{name}"
    raw,src=base.get(url,8_000_000);require(git_blob(raw)==sha,"pinned author repository blob mismatch")
    lines=raw.decode("utf-8-sig").splitlines()
    required={"variant_id","pval_nominal","slope","slope_se","maf"}
    ids=[i for i,line in enumerate(lines[:10]) if required<=set(next(csv.reader([line])))]
    require(len(ids)==1,"unique significant-table header")
    rr=csv.DictReader(io.StringIO("\n".join(lines[ids[0]:])))
    h=rr.fieldnames
    require(len(set(h))==len(h),"duplicate columns")
    trait="gene_id" if mode=="eQTL" else "protein_id"
    require(trait in h,"molecular trait key missing")
    count=0;traits=set();seen=set();target=[];region=[];exact=[];target_rows=0;duplicate=0;bad=0
    for row in rr:
        if not any(row.values()):continue
        count+=1;g=row[trait].split(".")[0];v=row["variant_id"]
        require(g and v,"blank molecular trait or variant")
        traits.add(g)
        if (g,v) in seen:duplicate+=1
        seen.add((g,v))
        is_target=(g==GENE or row.get("gene_name")=="GCC2")
        if is_target:target_rows+=1
        m=re.fullmatch(r"(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y):(\d+):([ACGT]+):([ACGT]+)",v)
        if not m:
            bad+=1;continue
        ch,p,a,b=m.groups();p=int(p);v2=f"{ch}:{p}:{a}:{b}"
        in_region=ch=="2" and LOW<=p<=HIGH
        if not (is_target or in_region or v2==LEAD):continue
        try:
            b0,se,pv,maf=[float(row[k]) for k in ("slope","slope_se","pval_nominal","maf")]
            require(all(math.isfinite(x) for x in (b0,se,pv,maf)) and se>0 and 0<=pv<=1 and 0<maf<=.5,"invalid summary statistic")
            if "ref" in row and "alt" in row:require(row["ref"]==a and row["alt"]==b,"variant allele columns differ")
        except (ValueError,KeyError,TypeError):bad+=1;continue
        z={"trait":g,"gene_name":row.get("gene_name"),"variant":v2,"reported_slope":b0,"se":se,"p":pv,"maf":maf}
        if is_target:target.append(z)
        if in_region:region.append(z)
        if v2==LEAD:exact.append(z)
    out={"status":"AUTHOR_SIGNIFICANT_TABLE_READBACK","assay":mode,"source":src,"git_blob":sha,"source_header":h,"rows":count,"molecular_traits":len(traits),
         "expected_paper_rows":n0,"expected_paper_traits":g0,"paper_denominator_match":count==n0 and len(traits)==g0,"duplicate_keys":duplicate,"invalid_or_unparsed_rows":bad,
         "GCC2_rows":target_rows,"GCC2_valid_summary_rows":len(target),"GCC2_region_rows":sum(LOW<=int(v["variant"].split(":")[1])<=HIGH and v["variant"].startswith("2:") for v in target),
         "GCC2_best_nominal":sorted(target,key=lambda v:v["p"])[:5],"exact_rs4584990_rows":exact,
         "regional_significant_rows":len(region),"regional_traits":len({v["trait"] for v in region})}
    out["limits"]=["Public significant-only supplementary table, pinned author Git blob; no individual genotype or expression data. Absence is not no association or no expression. Reported slopes as published, not assumed effect allele until source documentation; no colocalization or MR. Chinese molecular cohort not Korean-specific replication. Check paper denominator before scientific use."]
    return out
def inventory():
    raw,src=base.get("https://zenodo.org/api/records/19638349",4_000_000)
    d=json.loads(raw)
    require(d.get("id")==19638349 and d.get("doi")=="10.5281/zenodo.19638349","exact Zenodo record")
    fs=d.get("files");require(isinstance(fs,list) and fs,"record file inventory")
    result=[]
    for f in fs:
        name=f["key"];require(isinstance(name,str) and 0<len(name)<=240 and isinstance(f["size"],int) and f["size"]>=0,"file manifest entry")
        result.append({"name":name,"bytes":f["size"],"checksum":f.get("checksum")})
    return {"status":"CHINESE_SKIN_QTL_PUBLIC_RECORD_INVENTORY","source":src,"doi":d["doi"],"file_count":len(result),"files":result,"total_bytes":sum(f["bytes"] for f in result),
            "limits":["Only official record metadata; no genotype or participant-level downloads. Allpairs availability, exact compression/schema, effect orientation and full file integrity require separate qualification. No absence claim from significant tables."]}
def selftest():
    assert git_blob(b"")==hashlib.sha1(b"blob 0\0").hexdigest()
    assert len(FILES)==2 and FILES["eQTL"][2]==20999 and FILES["pQTL"][2]==1765
    assert LOW<int(LEAD.split(":")[1])<HIGH
    print("Pinned author reference and fixed GCC2 target checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["eQTL","pQTL","inventory","selftest"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=inventory() if a.mode=="inventory" else table(a.mode)
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    z.update(mode=a.mode,scientific_status="HOLD_FULL_REGIONAL_AND_ANCESTRY_VALIDATION")
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
