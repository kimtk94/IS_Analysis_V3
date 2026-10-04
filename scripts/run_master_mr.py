#!/usr/bin/env python3
"""Generic two-sample cis-MR engine for standardized summary-statistic TSV files.

Exposure schema (required):
protein_id, gene_symbol, rsid, effect_allele, other_allele, beta, se
Optional: eaf, n

Outcome schema (required):
rsid, effect_allele, other_allele, beta, se
Optional: eaf, n, effect_scale

Input exposure instruments must already be LD-pruned/independent.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from collections import Counter,defaultdict
from pathlib import Path

COMP={"A":"T","T":"A","C":"G","G":"C"}
PAL={frozenset(("A","T")),frozenset(("C","G"))}

def fnum(x):
    try:
        v=float(str(x).strip())
        return v if math.isfinite(v) else None
    except Exception:
        return None

def allele(x): return str(x or "").strip().upper()
def comp(a): return COMP.get(a) if len(a)==1 else None
def pal(a,b): return len(a)==1 and len(b)==1 and frozenset((a,b)) in PAL
def pnorm(z): return math.erfc(abs(z)/math.sqrt(2.0))

def read_tsv(path):
    op=gzip.open if Path(path).suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def write_tsv(path,rows,fields):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"wt",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def harmonize(x,y,maf_threshold=0.42,freq_tolerance=0.10):
    ea,oa=allele(x["effect_allele"]),allele(x["other_allele"])
    ya,yo=allele(y["effect_allele"]),allele(y["other_allele"])
    if not all((ea,oa,ya,yo)): return "missing_allele",None
    if pal(ea,oa):
        if frozenset((ya,yo))!=frozenset((ea,oa)): return "palindrome_allele_mismatch",None
        fx,fy=fnum(x.get("eaf")),fnum(y.get("eaf"))
        if fx is None or fy is None: return "palindrome_no_frequency",None
        if min(fx,1-fx)>maf_threshold: return "palindrome_ambiguous_maf",None
        keep=abs(fy-fx); flip=abs(fy-(1-fx))
        if min(keep,flip)>freq_tolerance: return "palindrome_frequency_mismatch",None
        if abs(keep-flip)<0.05: return "palindrome_frequency_ambiguous",None
        return ("palindrome_keep",1) if keep<flip else ("palindrome_flip",-1)
    if (ya,yo)==(ea,oa): return "direct",1
    if (ya,yo)==(oa,ea): return "swapped",-1
    cya,cyo=comp(ya),comp(yo)
    if cya is not None and cyo is not None:
        if (cya,cyo)==(ea,oa): return "strand",1
        if (cya,cyo)==(oa,ea): return "strand_swapped",-1
    return "allele_mismatch",None

def ivw(insts):
    den=sum(r["bx"]**2/r["sy"]**2 for r in insts)
    if den<=0:return None
    b=sum(r["bx"]*r["by"]/r["sy"]**2 for r in insts)/den
    se=math.sqrt(1/den)
    return b,se,pnorm(b/se)

def bh(ps):
    out=[None]*len(ps); valid=sorted((p,i) for i,p in enumerate(ps) if p is not None)
    m=len(valid); prev=1.0
    for rank in range(m,0,-1):
        p,i=valid[rank-1]
        q=min(prev,p*m/rank,1.0);out[i]=q;prev=q
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--exposure",type=Path,required=True)
    ap.add_argument("--outcome",type=Path,required=True)
    ap.add_argument("--harmonized-output",type=Path,required=True)
    ap.add_argument("--mr-output",type=Path,required=True)
    ap.add_argument("--ancestry",required=True)
    ap.add_argument("--phenotype",required=True)
    ap.add_argument("--min-fstat",type=float,default=10.0)
    ap.add_argument("--palindrome-maf-threshold",type=float,default=0.42)
    ap.add_argument("--freq-tolerance",type=float,default=0.10)
    a=ap.parse_args()

    required_x={"protein_id","gene_symbol","rsid","effect_allele","other_allele","beta","se"}
    required_y={"rsid","effect_allele","other_allele","beta","se"}

    exposures=list(read_tsv(a.exposure))
    outcomes=list(read_tsv(a.outcome))
    if exposures and not required_x.issubset(exposures[0]):
        raise SystemExit("exposure schema missing: "+",".join(sorted(required_x-set(exposures[0]))))
    if outcomes and not required_y.issubset(outcomes[0]):
        raise SystemExit("outcome schema missing: "+",".join(sorted(required_y-set(outcomes[0]))))

    out_by_rsid={r["rsid"]:r for r in outcomes if r.get("rsid")}
    groups=defaultdict(list); reasons=Counter(); harmonized=[]
    for x in exposures:
        bx,sx=fnum(x.get("beta")),fnum(x.get("se"))
        if bx is None or sx in (None,0) or bx==0:
            reasons["invalid_exposure"]+=1;continue
        fstat=(bx/sx)**2
        if fstat<a.min_fstat:
            reasons["weak_f"]+=1;continue
        y=out_by_rsid.get(x.get("rsid",""))
        if y is None:
            reasons["missing_outcome"]+=1;continue
        by,sy=fnum(y.get("beta")),fnum(y.get("se"))
        if by is None or sy in (None,0):
            reasons["invalid_outcome"]+=1;continue
        status,sign=harmonize(x,y,a.palindrome_maf_threshold,a.freq_tolerance)
        reasons[status]+=1
        if sign is None: continue
        row={
          "protein_id":x["protein_id"],"gene_symbol":x["gene_symbol"],"rsid":x["rsid"],
          "ancestry":a.ancestry,"phenotype":a.phenotype,
          "effect_allele":allele(x["effect_allele"]),"other_allele":allele(x["other_allele"]),
          "beta_exposure":bx,"se_exposure":sx,"exposure_eaf":x.get("eaf",""),"exposure_n":x.get("n",""),
          "beta_outcome":by*sign,"se_outcome":sy,"outcome_eaf":y.get("eaf",""),"outcome_n":y.get("n",""),
          "effect_scale":y.get("effect_scale",""),"f_stat":fstat,"harmonization":status
        }
        harmonized.append(row);groups[(x["protein_id"],x["gene_symbol"])].append(row)

    hfields=["protein_id","gene_symbol","rsid","ancestry","phenotype","effect_allele","other_allele",
             "beta_exposure","se_exposure","exposure_eaf","exposure_n","beta_outcome","se_outcome",
             "outcome_eaf","outcome_n","effect_scale","f_stat","harmonization"]
    write_tsv(a.harmonized_output,harmonized,hfields)

    results=[]
    for (protein,gene),insts in sorted(groups.items()):
        if len(insts)==1:
            r=insts[0]; b=r["beta_outcome"]/r["beta_exposure"]; se=r["se_outcome"]/abs(r["beta_exposure"])
            method="Wald ratio"; p=pnorm(b/se)
        else:
            x=ivw(insts)
            if x is None:continue
            b,se,p=x;method="IVW"
        scale=insts[0].get("effect_scale","")
        results.append({
          "protein_id":protein,"gene_symbol":gene,"ancestry":a.ancestry,"phenotype":a.phenotype,
          "n_instruments":len(insts),"method":method,"beta":b,"se":se,"p":p,
          "or":math.exp(b) if scale=="log_odds" else "",
          "min_f_stat":min(r["f_stat"] for r in insts)
        })
    qs=bh([r["p"] for r in results])
    for r,q in zip(results,qs):r["fdr_bh"]=q
    fields=["protein_id","gene_symbol","ancestry","phenotype","n_instruments","method","beta","se","p","fdr_bh","or","min_f_stat"]
    write_tsv(a.mr_output,results,fields)
    print(f"PASS proteins={len(results)} harmonized={len(harmonized)}")
    print("QC "+ " ".join(f"{k}={v}" for k,v in sorted(reasons.items())))

if __name__=="__main__":
    main()
