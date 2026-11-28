#!/usr/bin/env python3
"""Fixed nine-gene CKB cis-sentinel pilot using published independent-cohort EAS renal GWAS."""
import argparse
import collections
import gzip
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path

GENES = ("ACP1","CPVL","F12","GSTA1","GSTA3","HLA-E","INHBC","SDCCAG8","UMOD")
PIN = "f04eabc9bbb35753b50c4bb80bc574700784b81fcf3ba4ea4f9a08569218cb89"
EXPOSURE_SOURCE = "c6a6f4cfa58a2ee886a86482d7e08913c5c2aed4e42c3aa07f17c625d4fd6764"
SOURCES = {"egfr":(42774049,111041990,"TWB2_BBJ_eGFR_hg19_METAL_FUMA_noNA.gz"),
           "bun":(43218600,132381963,"TWB2_BBJ_BUN_hg19_METAL_FUMA_noNA.gz")}
COMPLEMENT = str.maketrans("ACGT","TGCA")

def digest(path, algorithm="sha256"):
    with path.open("rb") as fh:
        return hashlib.file_digest(fh,algorithm).hexdigest()

def orient(ea,oa,out_ea,out_oa):
    if {ea,oa} in ({"A","T"},{"C","G"}):
        return None
    if (ea,oa)==(out_ea,out_oa):
        return 1
    if (ea,oa)==(out_oa,out_ea):
        return -1
    return None

def wald(bx,sx,by,sy):
    if not all(math.isfinite(x) for x in (bx,sx,by,sy)) or bx==0 or min(sx,sy)<=0:
        raise ValueError("Invalid effect estimates")
    b=by/bx
    se=math.sqrt((sy/bx)**2+(by*sx/bx**2)**2)
    p=math.erfc(abs(b)/se/math.sqrt(2))
    return dict(mr_beta=b,mr_se_delta_independence=se,mr_ci95=[b-1.959964*se,b+1.959964*se],mr_p=p,
                mr_p_first_order=math.erfc(abs(by/sy)/math.sqrt(2)),
                mr_se_worst_covariance=abs(sy/bx)+abs(by*sx/bx**2))

def load_panel(path):
    if digest(path)!=PIN:
        raise ValueError("Frozen public derived artifact mismatch")
    with zipfile.ZipFile(path) as z:
        names=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="ckd.json"]
        if len(names)!=1 or z.testzip() is not None:
            raise ValueError("Invalid archive")
        p=json.loads(z.read(names[0]))
    if p["sha256"]!=EXPOSURE_SOURCE or set(p["candidate_genes"])!=set(GENES):
        raise ValueError("Unexpected exposure source")
    rows=p["cis_sentinels"]
    if len(rows)!=10 or collections.Counter(r["platform"] for r in rows)!={"Olink":5,"SomaScan":5}:
        raise ValueError("Frozen sentinel counts differ")
    for r in rows:
        if not r["basic_numeric_allele_qc"] or not r["F_ge_10"] or r["candidate_label"] not in GENES:
            raise ValueError("Invalid fixed sentinel")
    return rows

def checked_chain(root,name,manifest):
    path=root/name
    lines=(root/manifest).read_text().splitlines()
    md5s=[s.split()[0] for s in lines if len(s.split())==2 and s.split()[1].lstrip("*")==name]
    if len(md5s)!=1 or digest(path,"md5")!=md5s[0].lower():
        raise ValueError("UCSC chain MD5 mismatch "+name)
    return path

def lift(rows,forward,reverse):
    results=[]
    for r in rows:
        m=re.fullmatch(r"chr(\d+):(\d+):([ACGT]):([ACGT])",r["sentinel"])
        out=dict(r)
        if not m:
            results.append({**out,"status":"HOLD_NON_SNV"})
            continue
        chrom,pos,ref,alt=m.groups()
        pos=int(pos)
        ea,oa=r["effect_allele"],r["other_allele"]
        if {ea,oa}!={ref,alt}:
            raise ValueError("Exposure allele pair mismatch")
        if {ea,oa} in ({"A","T"},{"C","G"}):
            results.append({**out,"status":"HOLD_PALINDROMIC_NO_OUTCOME_EAF"})
            continue
        hits=forward.convert_coordinate("chr"+chrom,pos-1)
        if not hits or len(hits)!=1:
            results.append({**out,"status":"HOLD_NONUNIQUE_LIFT"})
            continue
        c,p,s,_=hits[0]
        back=reverse.convert_coordinate(c,p)
        if c!="chr"+chrom or not back or len(back)!=1 or back[0][:3]!=("chr"+chrom,pos-1,s):
            results.append({**out,"status":"HOLD_RECIPROCAL_LIFT"})
            continue
        if s not in ("+","-"):
            raise ValueError("Invalid chain strand")
        if s=="-":
            ea,oa,ref,alt=(a.translate(COMPLEMENT) for a in (ea,oa,ref,alt))
        results.append({**out,"chr37":int(c[3:]),"pos37":p+1,"strand":s,
                        "ea37":ea,"oa37":oa,"ref37":ref,"alt37":alt,"status":"MAPPED"})
    return results

def source_meta(root,trait):
    fid,size,name=SOURCES[trait]
    meta=json.loads((root/"figshare.json").read_text())
    if meta.get("id")!=24356587:
        raise ValueError("Unexpected Figshare article identity")
    match=[r for r in meta["files"] if r["id"]==fid and r["name"]==name and r["size"]==size]
    if len(match)!=1:
        raise ValueError("Exact source file metadata not verified")
    md5=match[0].get("computed_md5") or match[0].get("supplied_md5")
    path=root/"outcome.gz"
    if not re.fullmatch(r"[a-fA-F0-9]{32}",str(md5)) or path.stat().st_size!=size or digest(path,"md5")!=md5.lower():
        raise ValueError("Outcome file size or publisher MD5 mismatch")
    return dict(file_id=fid,name=name,bytes=size,md5=md5,sha256=digest(path),article=24356587)

def scan(path,rows):
    targets={(r["chr37"],r["pos37"]) for r in rows if r["status"]=="MAPPED"}
    found=collections.defaultdict(list)
    malformed=0
    total=0
    with gzip.open(path,"rt",encoding="utf-8") as f:
        h=f.readline().strip().split()
        required=("CHR","BP","SNP","Allele1","Allele2","Effect","StdErr","P")
        if not set(required)<=set(h) or len(h)!=len(set(h)):
            raise ValueError("Outcome header mismatch: "+repr(h))
        ix={k:h.index(k) for k in required}
        for line in f:
            if not line.strip():
                continue
            total+=1
            parts=line.split()
            if len(parts)!=len(h):
                malformed+=1
                continue
            try:
                coord=(int(parts[ix["CHR"]]),int(parts[ix["BP"]]))
            except ValueError:
                malformed+=1
                continue
            if coord in targets:
                found[coord].append({k:parts[i] for k,i in ix.items()})
    if total<1000000 or malformed:
        raise ValueError("Incomplete or malformed outcome: "+str((total,malformed)))
    for r in rows:
        if r["status"]!="MAPPED":
            continue
        rr=found[r["chr37"],r["pos37"]]
        matches=[]
        for q in rr:
            sign=orient(r["ea37"],r["oa37"],q["Allele1"].upper(),q["Allele2"].upper())
            if sign is not None:
                matches.append((q,sign))
        if len(matches)!=1:
            r["status"]="HOLD_MISSING_OR_AMBIGUOUS_OUTCOME"
            r["outcome_coordinate_rows"]=len(rr)
            continue
        q,sign=matches[0]
        try:
            by,sy,p=map(float,(q["Effect"],q["StdErr"],q["P"]))
            if not all(math.isfinite(x) for x in (by,sy,p)) or sy<=0 or not 0<=p<=1:
                raise ValueError("Invalid outcome")
            result = wald(float(r["beta"]),float(r["se"]),sign*by,sy)
            if set(result) & set(r):
                raise ValueError("Derived field collision with source evidence")
            r.update(result)
        except (ValueError,OverflowError):
            r["status"]="HOLD_OUTCOME_NUMERIC"
            continue
        # Nine pre-existing genes in primary Olink family; five frozen SomaScan rows
        # as separately labelled sensitivity family, including untestable rows.
        family=9 if r["platform"]=="Olink" else 5
        r.update(status="EXPLORATORY_WALD",rsid=q["SNP"],outcome_beta_aligned=sign*by,
                 outcome_se=sy,outcome_source_p=p,bonferroni_family=family,
                 bonferroni_p=min(1,r["mr_p"]*family))
    return total,malformed,h,rows

def selftest():
    assert orient("A","G","A","G")==1
    assert orient("A","G","G","A")==-1
    assert orient("A","T","A","T") is None
    assert orient("A","G","T","C") is None
    x=wald(0.5,0.02,-0.1,0.01)
    assert abs(x["mr_beta"]+0.2)<1e-12
    assert x["mr_se_delta_independence"]>0.02
    class Fake:
        def __init__(self,p): self.p=p
        def convert_coordinate(self,c,p): return [(c,p+self.p,"+",10)]
    r=lift([dict(sentinel="chr1:10:G:A",effect_allele="A",other_allele="G")],Fake(10),Fake(-10))
    assert r[0]["pos37"]==20
    assert r[0]["status"]=="MAPPED"
    assert not ({"beta", "se", "source_p"} & set(x))
    print("Nine synthetic contract assertions passed (including source field protection).")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("trait",choices=["egfr","bun","selftest"])
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--out",type=Path,default=Path("out.json"))
    a=p.parse_args()
    if a.trait=="selftest":
        selftest()
        return
    from pyliftover import LiftOver
    rows=load_panel(a.root/"panel.zip")
    meta=source_meta(a.root,a.trait)
    f=checked_chain(a.root,"hg38ToHg19.over.chain.gz","hg38.md5")
    b=checked_chain(a.root,"hg19ToHg38.over.chain.gz","hg19.md5")
    mapped=lift(rows,LiftOver(str(f)),LiftOver(str(b)))
    n,bad,h,results=scan(a.root/"outcome.gz",mapped)
    out=dict(status="RETROSPECTIVE_EXPLORATORY_WALD_NOT_COLOC",trait=a.trait,source=meta,
             source_rows=n,malformed_rows=bad,header=h,chain_sha256=[digest(f),digest(b)],
             candidate_family=list(GENES),exposure_source_sha256=EXPOSURE_SOURCE,
             exposure_artifact_sha256=PIN,results=results,
             limits=["Same EAS outcomes were used in earlier candidate prioritization.",
                     "CKB and TWB/BBJ are separate named cohorts, individual overlap not audited.",
                     "No outcome EAF, palindromes excluded; source-to-target reciprocal lift and exact alleles, no study-specific LD.",
                     "Platform estimates not independent; delta SE assumes zero exposure-outcome covariance.",
                     "No regional pQTL, colocalization, pleiotropy, protein-altering binding or clinical causality established."])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":
    main()
