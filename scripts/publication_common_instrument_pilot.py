#!/usr/bin/env python3
"""Outcome-blind, one-per-autosome LDL-C pilot, using pinned original GIGASTROKE outcomes.

This is a deliberately sparse sensitivity analysis, not genome-wide LD-clumped
MR or independent validation. No private files, access tokens, or server.
"""
import argparse
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import itertools
import json
import math
import re
import statistics
import subprocess
import tempfile
import urllib.request
from pathlib import Path

LDL_URL="https://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz"
SBP_URL="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-024-01714-w/MediaObjects/41588_2024_1714_MOESM3_ESM.xlsx"
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"
OUTCOMES={
 "LAA":("GCST90104542","71ea9fd7c06791a831861053b6a1cfbe",230156887,5784788),
 "CE":("GCST90104541","c252262586a6ea84b985f30dd7ab32be",263057516,6659793),
 "SVD":("GCST90104543","f41bc5f55b3cab238e305e3e86363f69",237118943,5974028),
}
DNA=set("ACGT")
COMPLEMENT=dict(zip("ACGT","TGCA"))
FAMILY=6
MIN_N=10
FIELDS=("SNP_hg19","rsid","A1","A2","beta","se","N","P-value","Freq.A1.1000G.EUR")

def digest(path, algorithm):
    h=hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def download(url,path,limit,timeout,retry=0):
    subprocess.run(["curl","--retry",str(retry),"--retry-delay","1","--retry-max-time","300","--fail","--location","--silent","--show-error","--connect-timeout","25","--max-time",str(timeout),"--max-filesize",str(limit),"-o",str(path),url],check=True)
    size=path.stat().st_size
    if not 0<size<=limit:raise ValueError("Invalid transfer size")
    return {"bytes":size,"sha256":digest(path,"sha256"),"md5":digest(path,"md5"),"url":url}

def load_ldl(path):
    """Full scan, source P and F only, no outcome knowledge."""
    candidates=collections.defaultdict(list)
    invalid=collections.Counter()
    n=0
    with gzip.open(path,"rt",encoding="utf8") as f:
        header=f.readline().strip().split()
        if len(set(header))!=len(header) or not set(FIELDS)<=set(header):
            raise ValueError("Unexpected official GLGC 2013 header: "+str(header))
        for line in f:
            n+=1
            a=line.split()
            if len(a)!=len(header):
                invalid["width"]+=1
                continue
            r=dict(zip(header,a))
            try:
                match=re.fullmatch(r"(?:chr)?(\d+):(\d+)",r["SNP_hg19"],re.I)
                if not match:raise ValueError("coordinate")
                chrom,pos=map(int,match.groups())
                if not 1<=chrom<=22 or pos<1:raise ValueError("non_autosome")
                if chrom==6 and 25_000_000<=pos<=34_000_000:raise ValueError("extended_MHC")
                ea,oa=r["A1"].upper(),r["A2"].upper()
                if ea not in DNA or oa not in DNA or ea==oa:raise ValueError("non_snv")
                if {ea,oa} in ({"A","T"},{"C","G"}):raise ValueError("palindrome")
                b,se,p,af,N=map(float,[r["beta"],r["se"],r["P-value"],r["Freq.A1.1000G.EUR"],r["N"]])
                if not all(math.isfinite(x) for x in (b,se,p,af,N)) or se<=0 or not 0<=p<=1 or N<50000 or not .01<=af<=.99:
                    raise ValueError("numeric_frequency_or_N")
                z=b/se
                if p>=5e-8 or z*z<29.7168:raise ValueError("not_genomewide_significant")
                if p==0 and abs(z)<37:raise ValueError("unresolved_zero_p")
                rsid=r["rsid"]
                if not re.fullmatch(r"rs\d+",rsid):raise ValueError("missing_rsid")
                candidates[chrom].append({"chrom":chrom,"pos":pos,"rsid":rsid,"ea":ea,"oa":oa,"bx":b,"sx":se,
                   "p":p,"f":z*z,"source_N":N,"ref_eaf":af})
            except (ValueError,OverflowError) as e:
                label=str(e)
                invalid[label if label in {"coordinate","non_autosome","extended_MHC","non_snv","palindrome","numeric_frequency_or_N","not_genomewide_significant","unresolved_zero_p","missing_rsid"} else "parse"]+=1
    panel=[]
    for chrom,rows in sorted(candidates.items()):
        rows.sort(key=lambda r:(r["p"],-r["f"],r["pos"],r["rsid"]))
        panel.append(rows[0])
    if len({r["rsid"] for r in panel})!=len(panel):raise ValueError("Duplicate selected rsID")
    # Verify selected coordinates uniquely identify one original source row.
    counts=collections.Counter()
    selected={(r["chrom"],r["pos"]):r for r in panel}
    with gzip.open(path,"rt",encoding="utf8") as f:
        h=f.readline().split()
        col=h.index("SNP_hg19")
        for line in f:
            a=line.split()
            if len(a)!=len(h):continue
            m=re.fullmatch(r"(?:chr)?(\d+):(\d+)",a[col],re.I)
            if m:
                key=tuple(map(int,m.groups()))
                if key in selected:counts[key]+=1
    panel=[r for r in panel if counts[(r["chrom"],r["pos"])]==1]
    if len(panel)<MIN_N:raise ValueError("Insufficient sparse exposure panel")
    return panel,{"source_rows":n,"eligible_genomewide_significant_rows":sum(map(len,candidates.values())),
                  "selected_before_coordinate_uniqueness":len(selected),"selected":len(panel),
                  "excluded_rows":dict(invalid),"header":header,
                  "selection":"One strongest reported P per autosome after exposure-only filters; ties by decreasing F then position/rsID. No outcome-dependent replacements. Exclude A/T, C/G and extended MHC. P=0 allowed only if |z|>=37; treated as underflow, not exact zero. No within-chromosome pair retained."}

def load_sbp(path):
    """Select only published primary SBP sentinel associations from ST3."""
    from openpyxl import load_workbook
    wb=load_workbook(path,read_only=True,data_only=True)
    ws=wb["ST3"]
    rows=ws.iter_rows(values_only=True)
    title=next(rows)[0]
    if "Known sentinel SNPs" not in title:raise ValueError("Unexpected ST3 title")
    next(rows)
    h=next(rows)
    required=("rsID","CHR","POS","Trait","Effect allele","Other Allele","EAF","Effect","SE","P-value","Neff")
    if not set(required)<=set(h):raise ValueError("Unexpected ST3 header")
    ix={k:h.index(k) for k in required}
    pool=collections.defaultdict(list)
    n=0
    invalid=collections.Counter()
    for a in rows:
        if not a or not a[0]:continue
        n+=1
        if a[ix["Trait"]]!="SBP":continue
        try:
            chrom,pos=int(a[ix["CHR"]]),int(a[ix["POS"]])
            rsid=str(a[ix["rsID"]])
            ea,oa=str(a[ix["Effect allele"]]).upper(),str(a[ix["Other Allele"]]).upper()
            b,se,p,af,N=(float(a[ix[k]]) for k in ("Effect","SE","P-value","EAF","Neff"))
            if not 1<=chrom<=22 or pos<1 or (chrom==6 and 25000000<=pos<=34000000):
                raise ValueError("chrom_pos")
            if not re.fullmatch(r"rs\d+",rsid) or ea not in DNA or oa not in DNA or ea==oa or {ea,oa} in ({"A","T"},{"C","G"}):
                raise ValueError("variant_identity")
            if not all(math.isfinite(x) for x in (b,se,p,af,N)) or se<=0 or not 0<=p<5e-8 or (b/se)**2<29.7168 or N<50000 or not .01<=af<=.99:
                raise ValueError("numeric")
            if p==0 and abs(b/se)<37:raise ValueError("zero_p")
            pool[chrom].append({"chrom":chrom,"pos":pos,"rsid":rsid,"ea":ea,"oa":oa,
                 "bx":b,"sx":se,"p":p,"f":(b/se)**2,"source_Neff":N,"source_eaf":af})
        except (ValueError,TypeError,OverflowError):
            invalid["excluded_primary_sbp"]+=1
    wb.close()
    panel=[min(a,key=lambda r:(r["p"],-r["f"],r["pos"],r["rsid"])) for _,a in sorted(pool.items())]
    if len({r["rsid"] for r in panel})!=len(panel) or len(panel)<MIN_N:
        raise ValueError("Insufficient unique SBP source panel")
    meta={"sheet":"ST3","source_rows":n,"eligible_primary_sbp_rows":sum(map(len,pool.values())),
          "selected_before_grch37_check":len(panel),"excluded":dict(invalid),
          "selection":"One strongest primary SBP sentinel per autosome; no outcome-based replacement. Require a single consistent observed source assembly and a unique GRCh37 rsID/allele mapping. No inference from PRS weights.",
          "original_selected_panel":panel}
    refs={}
    meta["reference_sources"]={}
    for assembly,host in (("GRCh37","https://grch37.rest.ensembl.org"),("GRCh38","https://rest.ensembl.org")):
        url=host+"/variation/homo_sapiens"
        try:
            req=urllib.request.Request(url,data=json.dumps({"ids":[r["rsid"] for r in panel]}).encode(),
                 headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"publication-pilot/1.0"})
            with urllib.request.urlopen(req,timeout=35) as response:
                raw=response.read(10*1024*1024)
            refs[assembly]=json.loads(raw)
            meta["reference_sources"][assembly]={"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"returned_ids":len(refs[assembly])}
        except (OSError,ValueError) as e:
            meta["reference_sources"][assembly]={"url":url,"error":type(e).__name__+": "+str(e)[:120]}
    if set(refs)!={"GRCh37","GRCh38"}:
        return [],{**meta,"selected":0,"status":"HOLD_REFERENCE_SOURCE"}
    verified,details=map_sbp_panel(panel,refs)
    return verified,{**meta,**details}

def map_sbp_panel(panel,refs):
    maps={}
    for build in ("GRCh37","GRCh38"):
        maps[build]={}
        for r in panel:
            hits=[]
            for m in refs[build].get(r["rsid"],{}).get("mappings",[]):
                allele=str(m.get("allele_string","")).upper().split("/")
                if m.get("assembly_name")!=build or str(m.get("seq_region_name"))!=str(r["chrom"]) or m.get("start")!=m.get("end") or len(allele)!=2:
                    continue
                if orientation(r["ea"],r["oa"],*allele) is None:
                    continue
                hits.append(int(m["start"]))
            if len(hits)==1:maps[build][r["rsid"]]=hits[0]
    matches={build:{r["rsid"] for r in panel if maps[build].get(r["rsid"])==r["pos"]}
             for build in maps}
    only37=matches["GRCh37"]-matches["GRCh38"]
    only38=matches["GRCh38"]-matches["GRCh37"]
    build="GRCh37" if only37 and not only38 else "GRCh38" if only38 and not only37 else None
    diagnostics={"source_coordinate_match_counts":{k:len(v) for k,v in matches.items()},
                 "unique_allele_mapping_counts":{k:len(v) for k,v in maps.items()},
                 "source_build_only37":len(only37),"source_build_only38":len(only38),
                 "source_build":build,
                 "mapping_rows":[{"rsid":r["rsid"],"source_pos":r["pos"],
                     "GRCh37":maps["GRCh37"].get(r["rsid"]),"GRCh38":maps["GRCh38"].get(r["rsid"])}
                     for r in panel]}
    if build is None or len(matches[build])<max(MIN_N,math.ceil(.8*len(panel))):
        return [],{**diagnostics,"selected":0,"status":"HOLD_UNRESOLVED_SOURCE_BUILD"}
    verified=[{**r,"source_pos":r["pos"],"source_build":build,"pos":maps["GRCh37"][r["rsid"]]}
              for r in panel if r["rsid"] in matches[build] and r["rsid"] in maps["GRCh37"]]
    if len(verified)<MIN_N:
        return [],{**diagnostics,"selected":len(verified),"status":"HOLD_TOO_FEW_GRCH37_MAPS"}
    return verified,{**diagnostics,"selected":len(verified),"status":"SOURCE_BUILD_AND_GRCH37_MAPS_VERIFIED"}

def orientation(ea,oa,a,b):
    if a not in DNA or b not in DNA or a==b or {ea,oa} in ({"A","T"},{"C","G"}):return None
    if (a,b)==(ea,oa):return (1,"same")
    if (a,b)==(oa,ea):return (-1,"swap")
    if (COMPLEMENT[a],COMPLEMENT[b])==(ea,oa):return (1,"complement")
    if (COMPLEMENT[a],COMPLEMENT[b])==(oa,ea):return (-1,"complement_swap")
    return None

def match_outcome(path,panel,expected_rows):
    targets={(str(r["chrom"]),str(r["pos"])):r for r in panel}
    seen=collections.Counter()
    hits={}
    excluded=collections.Counter()
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        header=f.readline().rstrip("\r\n").split("\t")
        required=("chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error")
        if not set(required)<=set(header) or len(set(header))!=len(header):
            raise ValueError("Unexpected GIGASTROKE original header: "+str(header))
        ix={name:header.index(name) for name in required}
        n=0
        for line in f:
            n+=1
            if (ix["chromosome"],ix["base_pair_location"])==(0,1):
                a=line.split("\t",2)
                if len(a)<3:raise ValueError("Malformed coordinate line")
                key=(a[0],a[1])
                if key not in targets:continue
            row=line.rstrip("\r\n").split("\t")
            if len(row)!=len(header):raise ValueError("Malformed target row")
            key=(row[ix["chromosome"]],row[ix["base_pair_location"]])
            if key not in targets:continue
            seen[key]+=1
            r=targets[key]
            match=orientation(r["ea"],r["oa"],row[ix["effect_allele"]].upper(),row[ix["other_allele"]].upper())
            if match is None:
                excluded["allele_mismatch"]+=1
                continue
            b,se=float(row[ix["beta"]]),float(row[ix["standard_error"]])
            if not math.isfinite(b) or not math.isfinite(se) or se<=0:
                excluded["invalid_effect"]+=1
                continue
            hits[key]={"by":match[0]*b,"sy":se,"orientation":match[1]}
    if n!=expected_rows:raise ValueError("Original outcome row count mismatch")
    hits={targets[k]["rsid"]:v for k,v in hits.items() if seen[k]==1}
    return hits,{"rows":n,"matched":len(hits),"coordinate_absent":len(set(targets)-set(seen)),
                 "coordinate_multiple":sum(v>1 for v in seen.values()),"excluded_target_rows":dict(excluded),
                 "header":header}

def read_outcome(name,panels,work):
    acc,md5,size,n=OUTCOMES[name]
    fn=acc+"_buildGRCh37.tsv.gz"
    p=work/fn
    source=download(BASE+acc+"/"+fn,p,size+1,240,retry=1)
    if source["bytes"]!=size or source["md5"]!=md5:
        raise ValueError(name+" original published MD5/size mismatch")
    hits={}
    meta={}
    for exposure,panel in panels.items():
        hits[exposure],meta[exposure]=match_outcome(p,panel,n)
    p.unlink()
    return name,hits,{"accession":acc,"source":source,"matched_by_exposure":meta}

def fit(x,y,se):
    if len(x)!=len(y) or len(y)!=len(se) or len(x)<2 or any(s<=0 for s in se):
        raise ValueError("Invalid regression inputs")
    w=[1/s**2 for s in se]
    den=math.fsum(a*a*b for a,b in zip(x,w))
    b=math.fsum(a*c*d for a,c,d in zip(x,y,w))/den
    q=math.fsum(d*(c-a*b)**2 for a,c,d in zip(x,y,w))
    phi=max(1,q/(len(x)-1))
    sigma=math.sqrt(phi/den)
    p=math.erfc(abs(b/sigma)/math.sqrt(2))
    return {"beta":b,"se":sigma,"or":math.exp(b),"ci95":[math.exp(b-1.95996398454*sigma),math.exp(b+1.95996398454*sigma)],
            "p":p,"q":q,"df":len(x)-1,"multiplicative_dispersion":phi}

def analyze(panel,outcomes):
    common=set.intersection(*(set(v) for v in outcomes.values()))
    rows=[r for r in panel if r["rsid"] in common]
    n=len(rows)
    if n<MIN_N:return {"status":"HOLD_TOO_FEW_COMMON_SNP","common":n}
    x=[r["bx"] for r in rows]
    ys={k:[v[r["rsid"]]["by"] for r in rows] for k,v in outcomes.items()}
    ss={k:[v[r["rsid"]]["sy"] for r in rows] for k,v in outcomes.items()}
    marginal={k:fit(x,ys[k],ss[k]) for k in OUTCOMES}
    contrasts=[]
    for a,b in itertools.combinations(OUTCOMES,2):
        diff=[u-v for u,v in zip(ys[a],ys[b])]
        # Unknown cross-outcome covariance per variant. Cauchy-Schwarz gives
        # Var(Y1-Y2) <= (SE1+SE2)^2; use same fixed weights in the difference regression.
        worst_se=[u+v for u,v in zip(ss[a],ss[b])]
        r=fit(x,diff,worst_se)
        loo=[fit(x[:i]+x[i+1:],diff[:i]+diff[i+1:],worst_se[:i]+worst_se[i+1:]) for i in range(n)]
        contrasts.append({"comparison":a+"-"+b,**r,"bonferroni_p":min(1,FAMILY*r["p"]),
                          "loo_beta_range":[min(z["beta"] for z in loo),max(z["beta"] for z in loo)],
                          "loo_max_bonferroni_p":min(1,FAMILY*max(z["p"] for z in loo))})
    return {"status":"SPARSE_COMMON_AUTOSOME_PILOT","common":n,"common_ids":[r["rsid"] for r in rows],
            "F_min":min(r["f"] for r in rows),"F_median":statistics.median(r["f"] for r in rows),
            "marginal_ivw":marginal,"paired_contrasts":contrasts,"family_n":FAMILY,
            "interpretation":"First-order zero-intercept inverse-variance regressions, multiplicative random effects with dispersion floor 1. Direct SNP-difference weights use worst-case per-SNP cross-outcome covariance, conditional on independent autosomal SNP estimates and exposure coefficients treated as fixed. Different estimand from subtracting separately weighted IVWs. One per chromosome prevents local LD pairs but does not prove independence, absence of pleiotropy, no cohort overlap, or validated causality. Three SBP contrasts untested. Exploratory alternative historical exposure; not original GLGC 2021/Evangelou replication."}

def inspect_sbp(out,work):
    from openpyxl import load_workbook
    p=work/"keaton2024.xlsx"
    result={"status":"HOLD_SBP_SOURCE_OR_SCHEMA"}
    try:
        source=download(SBP_URL,p,100*1024*1024,120)
        wb=load_workbook(p,read_only=True,data_only=True)
        sheets=list(wb.sheetnames)
        wanted=[]
        for name in sheets:
            # Inspect the published sentinel and score tables without assuming
            # that score weights are univariate GWAS coefficients.
            if re.search(r"(?:^|[^0-9])(?:3|4|26a)(?:$|[^0-9])",name,re.I):
                wanted.append(name)
        result={"status":"PUBLIC_WORKBOOK_READ_NEEDS_EFFECT_SCHEMA","source":source,"sheets":sheets,"candidate_tables":{}}
        for name in wanted:
            first=[]
            for i,row in enumerate(wb[name].iter_rows(values_only=True),1):
                if i>8:break
                first.append([str(v)[:160] if v is not None else "" for v in row[:70]])
            result["candidate_tables"][name]=first
        wb.close()
    except (OSError,ValueError,subprocess.CalledProcessError) as e:
        result["error"]=type(e).__name__+": "+str(e)[:150]
    write_json(out,result)

def write_json(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n",encoding="utf8")
    print(json.dumps(obj,allow_nan=False),flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("operation",choices=["ldl","sbp","both"])
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    with tempfile.TemporaryDirectory() as d:
        work=Path(d)
        if args.operation=="sbp":
            inspect_sbp(args.out,work)
            return
        result={"status":"HOLD_EXPOSURE_SOURCE","exposures":{}}
        panels={}
        for name in (("LDL-C","SBP") if args.operation=="both" else ("LDL-C",)):
            info={}
            result["exposures"][name]=info
            try:
                if name=="LDL-C":
                    p=work/"ldl.gz"
                    info["source"]=download(LDL_URL,p,200*1024*1024,150)
                    if info["source"]["sha256"]!="ddf75ed591ffcce8c3e81ae2f4c77f4428cd7dbf4628b24467aa8cf6e85e3dce":
                        raise ValueError("Original LDL source hash changed")
                    panel,meta=load_ldl(p)
                else:
                    p=work/"keaton.xlsx"
                    info["source"]=download(SBP_URL,p,100*1024*1024,120)
                    if info["source"]["sha256"]!="a21b3841e549e02f7321b500fabfc632537219795c455cab77fc1d032010ad18":
                        raise ValueError("Published SBP workbook hash changed")
                    panel,meta=load_sbp(p)
                info.update({"selection":meta,"panel":panel})
                if panel:panels[name]=panel
                else:info["error"]="Exposure selection or coordinate reference unresolved"
                p.unlink()
            except (ValueError,OSError,subprocess.CalledProcessError,KeyError,StopIteration) as e:
                info["error"]=type(e).__name__+": "+str(e)[:180]
        write_json(args.out,result)
        if panels:
            result["status"]="HOLD_OUTCOME_SOURCE"
            effects={}
            result["outcomes"]={}
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
                futures={pool.submit(read_outcome,k,panels,work):k for k in OUTCOMES}
                for future in concurrent.futures.as_completed(futures):
                    k=futures[future]
                    try:
                        _,hits,qc=future.result()
                        effects[k]=hits
                        result["matched_effects_by_outcome"]=effects
                        result["outcomes"][k]=qc
                    except (ValueError,OSError,subprocess.CalledProcessError) as e:
                        result["outcomes"][k]={"error":type(e).__name__+": "+str(e)[:180]}
                    write_json(args.out,result)
            if set(effects)==set(OUTCOMES):
                result["status"]="PUBLIC_SPARSE_PILOT_RESULTS_AVAILABLE"
                result["analyses"]={}
                for name,panel in panels.items():
                    matched={k:effects[k][name] for k in OUTCOMES}
                    info=analyze(panel,matched)
                    info["matched_effects"]=matched
                    info["interpretation"]=info.get("interpretation","").replace("Three SBP contrasts untested. Exploratory alternative historical exposure; not original GLGC 2021/Evangelou replication.","Six-comparison family includes both exposures even if one could not be executed. Exploratory alternate sources, not GLGC 2021/Evangelou replication.")
                    result["analyses"][name]=info
                if len(result["analyses"])!=2 or any(x["status"]!="SPARSE_COMMON_AUTOSOME_PILOT" for x in result["analyses"].values()):
                    result["status"]="PARTIAL_PUBLIC_SPARSE_PILOT"
        write_json(args.out,result)

if __name__=="__main__":
    main()
