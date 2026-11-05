#!/usr/bin/env python3
"""Dependency-free MR robustness analysis for harmonized cis-instrument tables.

Expected columns:
protein_id, gene_symbol, ancestry, phenotype, beta_exposure, beta_outcome,
se_outcome. Optional: rsid, se_exposure, outcome_n, exposure_n.

Methods:
- fixed-effect IVW through origin
- Cochran Q and chi-square survival P
- weighted median of ratio estimates
- MR-Egger weighted regression with multiplicative residual variance
- leave-one-out IVW influence summary
- optional Steiger direction when both exposure_n and outcome_n are available

This module is intended for harmonized, approximately independent instruments.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


def fnum(x):
    try:
        v=float(str(x).strip())
        return v if math.isfinite(v) else None
    except Exception:
        return None


def read_tsv(path: Path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as fh:
        yield from csv.DictReader(fh,delimiter="\t")


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    path.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"wt",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def p_norm(z):
    return math.erfc(abs(z)/math.sqrt(2.0))


def gammaincc(a, x):
    """Regularized upper incomplete gamma Q(a,x), Numerical Recipes style."""
    if x < 0 or a <= 0:
        return float("nan")
    if x == 0:
        return 1.0
    eps=3e-14; tiny=1e-300; itmax=1000
    gln=math.lgamma(a)
    if x < a+1.0:
        ap=a; summ=1.0/a; delta=summ
        for _ in range(itmax):
            ap += 1.0
            delta *= x/ap
            summ += delta
            if abs(delta) < abs(summ)*eps:
                break
        p=summ*math.exp(-x+a*math.log(x)-gln)
        return max(0.0,min(1.0,1.0-p))
    b=x+1.0-a; c=1.0/tiny; d=1.0/max(b,tiny); h=d
    for i in range(1,itmax+1):
        an=-i*(i-a)
        b += 2.0
        d=an*d+b
        if abs(d)<tiny: d=tiny
        c=b+an/c
        if abs(c)<tiny: c=tiny
        d=1.0/d
        delta=d*c
        h*=delta
        if abs(delta-1.0)<eps:
            break
    q=math.exp(-x+a*math.log(x)-gln)*h
    return max(0.0,min(1.0,q))


def chi2_sf(q, df):
    if df <= 0 or q < 0:
        return None
    return gammaincc(df/2.0,q/2.0)


def ivw(insts):
    good=[r for r in insts if r["bx"] not in (None,0) and r["by"] is not None and r["sy"] not in (None,0)]
    if not good: return None
    den=sum((r["bx"]**2)/(r["sy"]**2) for r in good)
    if den<=0: return None
    beta=sum(r["bx"]*r["by"]/(r["sy"]**2) for r in good)/den
    se=math.sqrt(1.0/den)
    q=sum(((r["by"]-beta*r["bx"])**2)/(r["sy"]**2) for r in good)
    df=len(good)-1
    return {"beta":beta,"se":se,"p":p_norm(beta/se),"q":q,"q_df":df,"q_p":chi2_sf(q,df)}


def weighted_median(insts):
    vals=[]
    for r in insts:
        if r["bx"] in (None,0) or r["by"] is None or r["sy"] in (None,0): continue
        ratio=r["by"]/r["bx"]
        se=r["sy"]/abs(r["bx"])
        if se<=0: continue
        vals.append((ratio,1.0/(se*se)))
    if not vals: return None
    vals.sort()
    total=sum(w for _,w in vals)
    acc=0.0
    for ratio,w in vals:
        acc += w
        if acc >= total/2.0:
            return ratio
    return vals[-1][0]


def egger(insts):
    good=[r for r in insts if r["bx"] is not None and r["by"] is not None and r["sy"] not in (None,0)]
    if len(good)<3: return None
    sw=sum(1/r["sy"]**2 for r in good)
    sx=sum(r["bx"]/r["sy"]**2 for r in good)
    sy=sum(r["by"]/r["sy"]**2 for r in good)
    sxx=sum(r["bx"]**2/r["sy"]**2 for r in good)
    sxy=sum(r["bx"]*r["by"]/r["sy"]**2 for r in good)
    det=sw*sxx-sx*sx
    if det<=0: return None
    intercept=(sxx*sy-sx*sxy)/det
    slope=(sw*sxy-sx*sy)/det
    rss=sum(((r["by"]-intercept-slope*r["bx"])**2)/(r["sy"]**2) for r in good)
    df=len(good)-2
    phi=max(rss/df,1.0) if df>0 else 1.0
    var_i=phi*sxx/det
    var_s=phi*sw/det
    se_i=math.sqrt(max(var_i,0))
    se_s=math.sqrt(max(var_s,0))
    return {
        "intercept":intercept,"intercept_se":se_i,
        "intercept_p":p_norm(intercept/se_i) if se_i>0 else None,
        "slope":slope,"slope_se":se_s,
        "slope_p":p_norm(slope/se_s) if se_s>0 else None,
        "rss":rss,"df":df,
    }


def loo(insts):
    full=ivw(insts)
    if full is None or len(insts)<3: return None
    estimates=[]
    for i,r in enumerate(insts):
        x=ivw(insts[:i]+insts[i+1:])
        if x is not None:
            estimates.append((r.get("rsid",""),x["beta"]))
    if not estimates: return None
    betas=[x[1] for x in estimates]
    max_item=max(estimates,key=lambda x:abs(x[1]-full["beta"]))
    sign_ok=sum((b==0 and full["beta"]==0) or b*full["beta"]>0 for b in betas)
    return {
      "loo_min_beta":min(betas),"loo_max_beta":max(betas),
      "loo_max_delta":abs(max_item[1]-full["beta"]),
      "loo_max_delta_rsid":max_item[0],
      "loo_sign_concordance":sign_ok/len(betas),
    }


def approx_r2(beta,se,n):
    if beta is None or se in (None,0) or n is None or n<=2: return None
    return beta*beta/(beta*beta + se*se*n)


def steiger(insts):
    votes=[]
    for r in insts:
        rx=approx_r2(r["bx"],r.get("sx"),r.get("nx"))
        ry=approx_r2(r["by"],r.get("sy"),r.get("ny"))
        if rx is not None and ry is not None:
            votes.append(rx>ry)
    if not votes: return None
    return {"n":len(votes),"forward_fraction":sum(votes)/len(votes),"forward":sum(votes)>len(votes)/2}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--min-instruments",type=int,default=2)
    a=ap.parse_args()

    groups=defaultdict(list)
    for row in read_tsv(a.input):
        key=(row.get("protein_id",""),row.get("gene_symbol",""),row.get("ancestry",""),row.get("phenotype",""))
        groups[key].append({
          "rsid":row.get("rsid",""),
          "bx":fnum(row.get("beta_exposure")),
          "sx":fnum(row.get("se_exposure")),
          "by":fnum(row.get("beta_outcome")),
          "sy":fnum(row.get("se_outcome")),
          "nx":fnum(row.get("exposure_n")),
          "ny":fnum(row.get("outcome_n")),
        })

    out=[]
    for (protein,gene,anc,pheno),insts in sorted(groups.items()):
        usable=[r for r in insts if r["bx"] not in (None,0) and r["by"] is not None and r["sy"] not in (None,0)]
        if len(usable)<a.min_instruments: continue
        x=ivw(usable); e=egger(usable); l=loo(usable); s=steiger(usable)
        row={
          "protein_id":protein,"gene_symbol":gene,"ancestry":anc,"phenotype":pheno,
          "n_instruments":len(usable),
          "ivw_beta":x["beta"] if x else "",
          "ivw_se":x["se"] if x else "",
          "ivw_p":x["p"] if x else "",
          "cochran_q":x["q"] if x else "",
          "cochran_q_df":x["q_df"] if x else "",
          "cochran_q_p":x["q_p"] if x else "",
          "weighted_median_beta":weighted_median(usable),
          "egger_slope":e["slope"] if e else "",
          "egger_slope_se":e["slope_se"] if e else "",
          "egger_slope_p":e["slope_p"] if e else "",
          "egger_intercept":e["intercept"] if e else "",
          "egger_intercept_se":e["intercept_se"] if e else "",
          "egger_intercept_p":e["intercept_p"] if e else "",
          "loo_min_beta":l["loo_min_beta"] if l else "",
          "loo_max_beta":l["loo_max_beta"] if l else "",
          "loo_max_delta":l["loo_max_delta"] if l else "",
          "loo_max_delta_rsid":l["loo_max_delta_rsid"] if l else "",
          "loo_sign_concordance":l["loo_sign_concordance"] if l else "",
          "steiger_n":s["n"] if s else 0,
          "steiger_forward_fraction":s["forward_fraction"] if s else "",
          "steiger_forward":int(s["forward"]) if s else "",
        }
        out.append(row)

    fields=[
      "protein_id","gene_symbol","ancestry","phenotype","n_instruments",
      "ivw_beta","ivw_se","ivw_p","cochran_q","cochran_q_df","cochran_q_p",
      "weighted_median_beta","egger_slope","egger_slope_se","egger_slope_p",
      "egger_intercept","egger_intercept_se","egger_intercept_p",
      "loo_min_beta","loo_max_beta","loo_max_delta","loo_max_delta_rsid","loo_sign_concordance",
      "steiger_n","steiger_forward_fraction","steiger_forward"
    ]
    write_tsv(a.output,out,fields)
    print(f"PASS groups={len(out)} output={a.output}")

if __name__=="__main__":
    main()
