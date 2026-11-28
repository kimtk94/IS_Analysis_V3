#!/usr/bin/env python3
"""Frozen 241-gene panel: external 7-day CPM robustness, not gene-wise replication p-values."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import math
import re
import statistics
import subprocess
import tempfile
import zipfile
from pathlib import Path

PIN="5ad39d6fb293130492e66c713ff7d155d03d752f31cac84c98f78c42671bd606"
CPM_PIN="717f4c9861de4d6b169439ef56a781f50d93ec918c456a9d0e61873e4a6e5f09"
CPM="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301517/suppl/GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"
STRAINS={"TM62":4,"TM70":2,"TM72":1}
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_panel(path):
    if sha(path)!=PIN:raise ValueError("Frozen candidate archive hash mismatch")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError("Archive CRC")
        names=z.namelist()
        def table(name):
            match=[n for n in names if n.rsplit("/",1)[-1]==name]
            if len(match)!=1:raise ValueError("Missing or ambiguous archive table "+name)
            return list(csv.DictReader(io.StringIO(z.read(match[0]).decode()),delimiter="\t"))
        panel=table("strong_241.tsv")
        paths=table("reactome_all_terms.tsv")
    result={}
    for r in panel:
        gid=r["gene_id"]
        if not re.fullmatch(r"ENSG\d+",gid) or gid in result:raise ValueError("Duplicate or invalid candidate")
        b=[float(r[k]) for k in ("donor_log2fc","pool_log2fc","external_log2fc")]
        f=[float(r[k]) for k in ("donor_fdr","pool_fdr","external_fdr")]
        if any(not math.isfinite(x) for x in b+f) or not all(0<=x<.05 for x in f) or min(abs(b[0]),abs(b[1]))<1 or b[0]*b[1]<=0 or b[0]*b[2]<=0:
            raise ValueError("Candidate threshold changed")
        result[gid]={"gene_id":gid,"symbol":r["donor_symbol"],"acute_log2fc":b[0],"pool_log2fc":b[1],"published_48h_log2fc":b[2]}
    if len(result)!=241:raise ValueError("Expected frozen 241 genes")
    ecm=[r for r in paths if r["pathway_id"]=="R-HSA-1474244"]
    if len(ecm)!=1:raise ValueError("Expected exact ECM organization pathway")
    symbols=set(ecm[0]["overlap_symbols"].split(";"))
    if len(symbols)!=16:raise ValueError("Expected frozen 16 ECM genes")
    for r in result.values():
        r["ecm_member"]=r["symbol"].strip().upper() in symbols
    if sum(r["ecm_member"] for r in result.values())!=16:raise ValueError("ECM symbol mapping not unique")
    return result

def parse_cpm(path,panel):
    if sha(path)!=CPM_PIN:raise ValueError("Original 7-day CPM hash mismatch")
    results={}
    seen=set()
    n=0
    with gzip.open(path,"rt",encoding="utf8",newline="") as h:
        rd=csv.DictReader(h,delimiter="\t")
        required={""}
        for s,k in STRAINS.items():
            required.update(f"{s}_{arm}_{i}" for arm in ("Control","Dex") for i in range(1,k+1))
        if not required<=set(rd.fieldnames or []) or len(rd.fieldnames)!=len(set(rd.fieldnames)):
            raise ValueError("Original CPM column contract")
        for row in rd:
            n+=1
            gid=row[""].split(".")[0]
            if gid not in panel:continue
            if gid in seen:raise ValueError("Duplicated candidate ID")
            seen.add(gid)
            means={}
            for s,k in STRAINS.items():
                for arm in ("Control","Dex"):
                    v=[float(row[f"{s}_{arm}_{i}"]) for i in range(1,k+1)]
                    if any(not math.isfinite(x) or x<0 for x in v):raise ValueError("Invalid CPM")
                    means[s,arm]=statistics.mean(v)
            effects={s:math.log2((means[s,"Dex"]+1)/(means[s,"Control"]+1)) for s in STRAINS}
            b=panel[gid]["acute_log2fc"]
            aligned=[v*(1 if b>0 else -1) for v in effects.values()]
            # Both-arm equal-strain mean expression; do not select strains by outcome.
            min_cpm=min(statistics.mean(means[s,arm] for s in STRAINS) for arm in ("Control","Dex"))
            results[gid]={**panel[gid],"status":"observed","min_arm_mean_cpm":min_cpm,
                "mean_7d_log2_cpm_ratio":statistics.mean(effects.values()),
                "n_strains_same_sign":sum(x>0 for x in aligned),
                "n_strains_same_sign_abs_ge_025":sum(x>=.25 for x in aligned),
                "all_leave_one_strain_out_same_sign":all(statistics.mean(aligned[:i]+aligned[i+1:])>0 for i in range(len(aligned))),
                "min_signed_strain_effect":min(aligned),
                **{s+"_log2_cpm_ratio":v for s,v in effects.items()}}
    if n!=85298:raise ValueError("Unexpected original CPM row count")
    return results,n

def summarize(panel,observed):
    groups={"all_241":list(panel),"ECM_16":[k for k,v in panel.items() if v["ecm_member"]]}
    summary={}
    for name,ids in groups.items():
        mapped=[observed[k] for k in ids if k in observed]
        stats=[]
        for c in (0,1,5):
            rows=[r for r in mapped if r["min_arm_mean_cpm"]>=c]
            stats.append({"min_both_arm_mean_cpm":c,"evaluable":len(rows),
                "aggregate_same_sign":sum(r["mean_7d_log2_cpm_ratio"]*r["acute_log2fc"]>0 for r in rows),
                "all_three_strains_same_sign":sum(r["n_strains_same_sign"]==3 for r in rows),
                "all_three_strains_same_sign_abs_ge_025":sum(r["n_strains_same_sign_abs_ge_025"]==3 for r in rows),
                "all_leave_one_strain_out_same_sign":sum(r["all_leave_one_strain_out_same_sign"] for r in rows)})
        summary[name]={"frozen_total":len(ids),"mapped":len(mapped),"unavailable":len(ids)-len(mapped),"strata":stats}
    return summary

def get_source(p):
    attempts=[]
    for url in (CPM,"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE301517&format=file&file=GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"):
        try:
            subprocess.run(["curl","-fLsS","--connect-timeout","20","--max-time","80","--max-filesize","26000000","-o",str(p),url],check=True)
            if sha(p)!=CPM_PIN:raise ValueError("Source hash mismatch")
            return url,attempts
        except (OSError,ValueError,subprocess.CalledProcessError) as e:
            attempts.append({"url":url,"error":type(e).__name__})
    return None,attempts

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--panel",type=Path,required=True)
    a.add_argument("--cpm", type=Path, help="Replay the original SHA-256-pinned GEO gzip from disk.")
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    panel=read_panel(args.panel)
    summary={"status":"HOLD_7DAY_SOURCE","frozen_panel_sha256":PIN,"frozen_cpm_sha256":CPM_PIN,
       "interpretation":"Frozen earlier 241-gene/16-ECM panel tested against published GSE301517 normalized CPM. Technical repeats averaged within each of three strains and strains weighted equally. This extends a previously inspected cohort, not a wholly unexamined holdout; strain donor independence not verified. Gene-level FDR and causal time/dose inference unavailable. Expression cutoffs are descriptive sensitivity strata, not selected thresholds. No gene-independence binomial p-value."}
    with tempfile.TemporaryDirectory() as d:
        if args.cpm is None:
            p = Path(d) / "cpm.gz"
            url, attempts = get_source(p)
            summary["cpm_input"] = "download"
        else:
            p = args.cpm
            url, attempts = CPM, []
            summary["cpm_input"] = "local_verified_copy"
        summary["attempts"]=attempts
        if url is not None:
            obs,n=parse_cpm(p,panel)
            summary.update({"status":"FROZEN_PANEL_DIRECTIONAL_AUDIT_COMPLETE","source_url":url,"source_rows":n,"groups":summarize(panel,obs)})
            summary["ecm_gene_results"]=[{k:v for k,v in r.items() if k in ("symbol","acute_log2fc","mean_7d_log2_cpm_ratio","n_strains_same_sign","n_strains_same_sign_abs_ge_025","min_arm_mean_cpm")} for r in obs.values() if r["ecm_member"]]
            # Keep all frozen candidates, including missing rows, in the output.
            fields=list(next(iter(obs.values()))) if obs else list(next(iter(panel.values())))+["status"]
            with (args.out/"frozen_241_directional_results.tsv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
                w.writeheader()
                for k in sorted(panel):
                    w.writerow(obs.get(k,{**panel[k],"status":"not_mapped"}))
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
