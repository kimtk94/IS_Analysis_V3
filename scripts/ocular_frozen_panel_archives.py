#!/usr/bin/env python3
"""Analyze the frozen panel against real archived alternate-annotation and 7-day aggregates."""
import argparse
import collections
import csv
import gzip
import io
import json
import math
import zipfile
from pathlib import Path
from ocular_frozen_241_external import read_panel,sha

SOURCES={
 "score050":("83dab7522d27db251b305024f02c47c4fc06620f57f85c58848f7e8a510af009","13041944372"),
 "dual55":("59790368f2b81c094c6ffd30e3e93b22a1bcf150adb59f8e1d47e05d2a8fda29","13044762180"),
 "seven_day":("77396ceb334557806177f8a4d512c9aec967ef2ccc6d17b38bbdcf76b4d5cbd2","13037335888"),
}
def archive(path,role):
    if sha(path)!=SOURCES[role][0]:raise ValueError(role+" SHA256 mismatch")
    z=zipfile.ZipFile(path)
    if z.testzip() is not None:
        z.close()
        raise ValueError(role+" CRC mismatch")
    return z

def table(z,suffix):
    names=[n for n in z.namelist() if n==suffix or n.endswith("/"+suffix)]
    if len(names)!=1:raise ValueError("Missing/ambiguous "+suffix+": "+str(z.namelist()[:25]))
    data=z.read(names[0])
    if suffix.endswith(".gz"):data=gzip.decompress(data)
    return list(csv.DictReader(io.StringIO(data.decode("utf8")),delimiter="\t"))

def gene_stats(rows,n):
    d={}
    for r in rows:
        g=r["gene_id"].split(".")[0]
        b,f=float(r["logFC"]),float(r["FDR"])
        if not g.startswith("ENSG") or g in d or not math.isfinite(b) or not math.isfinite(f) or not 0<=f<=1:
            raise ValueError("Invalid alternate-annotation gene stats")
        d[g]=(b,f)
    if len(d)!=n:raise ValueError("Alternate model gene count mismatch")
    return d

def chronic_stats(z,panel):
    """Preserve missing panel genes, not an extrapolation from previously selected EMT/P9."""
    d={}
    a=table(z,"GSE301517_EMT_shared_3strain_direction.tsv")
    b=table(z,"GSE301517_16_p9_panel_3strain_vs_48h.tsv")
    if len(a)!=194 or len(b)!=16:raise ValueError("Historical external set mismatch")
    for rows,kind in ((a,"EMT"),(b,"P9")):
        seen=set()
        for r in rows:
            gid=r["gene_id"].split(".")[0]
            if gid in seen:raise ValueError("Duplicate external ID")
            seen.add(gid)
            if gid not in panel:continue
            if kind=="EMT":
                primary=float(r["primary_48h_log2FC"])
                effect=float(r["third_7d_CPM_log2ratio"])
                n=int(r["strains_3_matching_primary"])
            else:
                primary=float(r["primary_500nM_48h_log2FC"])
                effect=float(r["third_500nM_7d_CPM_log2ratio"])
                n=int(r["three_strain_sign_match_n"])
            if not math.isclose(primary,panel[gid]["acute_log2fc"],rel_tol=0,abs_tol=1e-5) or not math.isfinite(effect) or n not in range(4):
                raise ValueError("Frozen external reference mismatch "+gid)
            if gid in d:
                old=d[gid]
                if old["n_strains_same_sign"]!=n or not math.isclose(old["chronic_log2ratio"],effect,abs_tol=1e-5,rel_tol=0):
                    raise ValueError("Overlapping historical external source discrepancy")
                old["source_subset"]="EMT+P9"
            else:
                d[gid]={"chronic_log2ratio":effect,"n_strains_same_sign":n,"source_subset":kind}
    return d

def compare(panel,models,chronic):
    rows=[]
    for g,p in sorted(panel.items()):
        b=p["acute_log2fc"]
        record=dict(p)
        for name,model in models.items():
            v=model.get(g)
            record[name+"_log2fc"]=v[0] if v else None
            record[name+"_fdr"]=v[1] if v else None
        values=[m[g] for m in models.values() if g in m]
        record["all_three_alternate_models_tested"]=len(values)==3
        record["all_alternate_signed_fdr05"]=len(values)==3 and all(b*x>0 and f<.05 for x,f in values)
        record["all_alternate_signed_fdr05_abs1"]=record["all_alternate_signed_fdr05"] and all(abs(x)>=1 for x,f in values)
        record["chronic_covered"]=g in chronic
        record.update(chronic.get(g,{"chronic_log2ratio":None,"n_strains_same_sign":None,"source_subset":None}))
        record["chronic_mean_same_sign"]=g in chronic and chronic[g]["chronic_log2ratio"]*b>0
        record["chronic_3of3_same_sign"]=g in chronic and chronic[g]["n_strains_same_sign"]==3
        rows.append(record)
    out={}
    for name,group in (("all_241",rows),("ECM_16",[r for r in rows if r["ecm_member"]])):
        m={}
        for k in models:
            tested=[r for r in group if r[k+"_fdr"] is not None]
            sign=[r for r in tested if r[k+"_log2fc"]*r["acute_log2fc"]>0]
            sig=[r for r in sign if r[k+"_fdr"]<.05]
            m[k]={"tested":len(tested),"sign_matches":len(sign),"sign_and_fdr05":len(sig),
                  "sign_fdr05_abs1":sum(abs(r[k+"_log2fc"])>=1 for r in sig)}
        out[name]={"n":len(group),"per_model":m,**{k:sum(r[k] for r in group) for k in ("all_three_alternate_models_tested","all_alternate_signed_fdr05","all_alternate_signed_fdr05_abs1","chronic_covered","chronic_mean_same_sign","chronic_3of3_same_sign")},
          "both_annotation_and_chronic_3of3":sum(r["all_alternate_signed_fdr05"] and r["chronic_3of3_same_sign"] for r in group)}
    return rows,out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("panel","score050","dual55","seven-day","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    panel=read_panel(a.panel)
    with archive(a.score050,"score050") as z:
        s=gene_stats(table(z,"p5_all_gene_results.tsv.gz"),16728)
    with archive(a.dual55,"dual55") as z:
        b=gene_stats(table(z,"BeamA_075/p5_all_gene_results.tsv.gz"),16051)
        f=gene_stats(table(z,"TM_Fibro/p5_all_gene_results.tsv.gz"),14003)
    with archive(a.seven_day,"seven_day") as z:
        chronic=chronic_stats(z,panel)
    models={"score050_57":s,"beam55":b,"fibro55":f}
    rows,summary=compare(panel,models,chronic)
    result={"status":"ARCHIVED_GENE_LEVEL_SENSITIVITY_COMPLETE","sources":SOURCES,"groups":summary,
       "interpretation":"All alternate annotation models use the same source cohort, not independent replication. 7-day coverage is the exact intersection with previously exported EMT194 and P9-16, not the full 241-panel test; no extrapolation or gene-wise P/FDR. Source library/treatment confounding and strain donor independence remain unresolved. No new wet-lab, clinical or causal inference.",
       "ecm_genes":[r for r in rows if r["ecm_member"]]}
    with (a.out/"candidate_evidence.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (a.out/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
