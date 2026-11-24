#!/usr/bin/env python3
"""P5m: exploratory opposite-sign TM evidence tiering + official 1h/6h/48h trajectories.

Reads only P5l gene-level aggregate audit and 2026 IOVS public workbook.
No donor counts, no raw FASTQ, no interpretation as proven dose reversal.
"""
import argparse
import csv
import io
import json
import math
import urllib.request
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from openpyxl import load_workbook

ZIP_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/supplementaryFiles"
XLSX="iovs-67-2-13_s003.xlsx"
SHEETS={
    "1h":"SuppTable64_TM1Hr_All_DEGs",
    "6h":"SuppTable41_TM6Hrs_All_DEGs",
    "48h":"SuppTable4_TM2Days_All_DEGs",
}
EXPECTED={"1h":17077,"6h":17132,"48h":13313}

def number(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):
        return None

def boolean(x):
    if x not in ("True","False"):raise ValueError("EXPECTED_TRUE_OR_FALSE "+str(x))
    return x=="True"

def sign(x):
    return 1 if x>0 else -1 if x<0 else 0

def load_discordant(filename):
    with open(filename,encoding="utf8",newline="") as f:
        dr=csv.DictReader(f,delimiter="\t")
        req={"gene_id","gene_symbol","primary_log2fc","primary_FDR",
             "external_2026_TM_log2fc","external_2026_TM_FDR",
             "poolunit_log2fc","poolunit_FDR","poolunit_same_as_primary",
             "triangulated_conflict"}
        if not req<=set(dr.fieldnames or []):
            raise ValueError("P5L_SOURCE_SCHEMA_MISMATCH")
        rows=list(dr)
    if len(rows)!=96 or len({r["gene_id"] for r in rows})!=96:
        raise ValueError("EXPECTED_96_UNIQUE_OPPOSITE_SIGN_GENES")
    for r in rows:
        for c in ("primary_log2fc","primary_FDR","external_2026_TM_log2fc",
                  "external_2026_TM_FDR","poolunit_log2fc","poolunit_FDR"):
            n=number(r[c])
            if n is None:raise ValueError("NONNUMERIC_"+c)
            r[c]=n
        for c in ("poolunit_same_as_primary","triangulated_conflict"):
            r[c]=boolean(r[c])
        if not (r["primary_FDR"]<.05 and r["external_2026_TM_FDR"]<.05
                and r["primary_log2fc"]*r["external_2026_TM_log2fc"]<0):
            raise ValueError("DISCORDANT_CONTRACT_CHANGED_"+r["gene_id"])
    return rows

def fetch_original_xlsx():
    req=urllib.request.Request(ZIP_URL,
          headers={"User-Agent":"ocular-public-tm-dose-time-audit/1.0"})
    with urllib.request.urlopen(req,timeout=120) as resp:
        data=resp.read(65*1024*1024+1)
    if len(data)>65*1024*1024:raise ValueError("SUPPLEMENTARY_ARCHIVE_TOO_LARGE")
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        paths=[p for p in z.namelist() if p.rsplit("/",1)[-1]==XLSX]
        if len(paths)!=1:raise ValueError("OFFICIAL_SUPPLEMENT_WORKBOOK_MISSING")
        wb=z.read(paths[0])
    if len(wb)>35*1024*1024:raise ValueError("WORKBOOK_TOO_LARGE")
    return wb,len(data)

def read_timepoint(sheet, subset):
    it=sheet.iter_rows(values_only=True)
    colnames=[str(x).strip() if x is not None else "" for x in next(it)]
    needed={"ID","log2FoldChange","padj"}
    if not needed<=set(colnames):raise ValueError(sheet.title+" REQUIRED_FIELDS_MISSING")
    idx={v:i for i,v in enumerate(colnames)}
    stored={}
    count=0
    for row in it:
        if not row or row[idx["ID"]] is None:continue
        count+=1
        gene=str(row[idx["ID"]]).strip().split(".")[0].upper()
        if gene not in subset:continue
        if gene in stored:raise ValueError("DUPLICATE_IN_TIME_SHEET "+gene)
        fc=number(row[idx["log2FoldChange"]])
        padj=number(row[idx["padj"]])
        stored[gene]={"lfc":fc,"fdr":padj}
    return stored,count

def classify(d):
    p=abs(d["primary_log2fc"])
    e=abs(d["external_2026_TM_log2fc"])
    q=abs(d["poolunit_log2fc"])
    if not (d["triangulated_conflict"] and d["poolunit_same_as_primary"]):
        return "C_POOL_NOT_SIGNIFICANT_OR_SIGN_CHANGED"
    if p>=1 and e>=.5 and q>=.5:
        return "A_STRONG_CONFLICT_HYPOTHESIS"
    if p>=.5 and e>=.5 and q>=.5:
        return "B_MODERATE_CONFLICT_HYPOTHESIS"
    if p>=.25 and e>=.5 and q>=.25:
        return "C_SMALLER_BUT_STABLE_CONFLICT"
    return "D_SMALL_EFFECT_OR_WEAK_EXTERNAL_EFFECT"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--p5l",required=True)
    ap.add_argument("--p5b-loo",required=True)
    ap.add_argument("--output-dir",required=True)
    args=ap.parse_args()
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    result={"status":"HOLD","created_utc":datetime.now(timezone.utc).isoformat()}
    try:
        disc=load_discordant(args.p5l)
        # Most P5l conflicts are outside P5b top-100; do NOT fabricate LOO verification.
        with open(args.p5b_loo,encoding="utf8",newline="") as f:
            it=csv.DictReader(f,delimiter="\t")
            if not {"gene_id","excluded_pool"}<=set(it.fieldnames or []):
                raise ValueError("LOO_FIELDS_MISSING")
            loo_rows=list(it)
        loo_coverage=Counter(r["gene_id"] for r in loo_rows)
        candidates={x["gene_id"] for x in disc}
        overlap=candidates&set(loo_coverage)
        for gene in overlap:
            if loo_coverage[gene]!=9:raise ValueError("PARTIAL_LOO_GENE_"+gene)
        official_bytes,source_bytes=fetch_original_xlsx()
        wb=load_workbook(io.BytesIO(official_bytes),read_only=True,data_only=True)
        if not set(SHEETS.values())<=set(wb.sheetnames):
            raise ValueError("OFFICIAL_TM_TIMES_NOT_ALL_PRESENT")
        by_time={}
        for time,sheet in SHEETS.items():
            by_time[time],observed=read_timepoint(wb[sheet],candidates)
            if observed!=EXPECTED[time]:
                raise ValueError("OFFICIAL_SHEET_ROW_DRIFT_"+time+"_"+str(observed))
        wb.close()
        output=[]
        for g in disc:
            gid=g["gene_id"];t={time:by_time[time].get(gid) for time in SHEETS}
            tier=classify(g)
            enriched={
                "gene_id":gid,"gene_symbol":g["gene_symbol"],
                "priority_tier":tier,
                "primary_BeamA_500nM_48h_log2fc":g["primary_log2fc"],
                "primary_BeamA_FDR":g["primary_FDR"],
                "conservative_poolunit_log2fc":g["poolunit_log2fc"],
                "conservative_poolunit_FDR":g["poolunit_FDR"],
                "external_TM_100nM_48h_log2fc":g["external_2026_TM_log2fc"],
                "external_TM_100nM_48h_FDR":g["external_2026_TM_FDR"],
                "external_TM_100nM_1h_log2fc":t["1h"]["lfc"] if t["1h"] else "",
                "external_TM_100nM_1h_FDR":t["1h"]["fdr"] if t["1h"] and t["1h"]["fdr"] is not None else "",
                "external_TM_100nM_6h_log2fc":t["6h"]["lfc"] if t["6h"] else "",
                "external_TM_100nM_6h_FDR":t["6h"]["fdr"] if t["6h"] and t["6h"]["fdr"] is not None else "",
                "1h_sign_matches_external_48h":(
                    sign(t["1h"]["lfc"])==sign(t["48h"]["lfc"])
                    if t["1h"] and t["48h"] and t["1h"]["lfc"] is not None and t["48h"]["lfc"] is not None else ""),
                "6h_sign_matches_external_48h":(
                    sign(t["6h"]["lfc"])==sign(t["48h"]["lfc"])
                    if t["6h"] and t["48h"] and t["6h"]["lfc"] is not None and t["48h"]["lfc"] is not None else ""),
                "has_formal_9x_donor_LOO":gid in overlap,
                "mechanistic_novelty":"UNVERIFIED",
                "causal_dose_reversal":"UNVERIFIED",
            }
            if t["48h"] and t["48h"]["lfc"] is not None:
                if abs(t["48h"]["lfc"]-g["external_2026_TM_log2fc"])>1e-5:
                    raise ValueError("REPEATED_48H_VALUE_CHANGED_"+gid)
            output.append(enriched)
        grades=["A_STRONG_CONFLICT_HYPOTHESIS","B_MODERATE_CONFLICT_HYPOTHESIS",
                "C_SMALLER_BUT_STABLE_CONFLICT","D_SMALL_EFFECT_OR_WEAK_EXTERNAL_EFFECT",
                "C_POOL_NOT_SIGNIFICANT_OR_SIGN_CHANGED"]
        output.sort(key=lambda v:(grades.index(v["priority_tier"]),
                                  -abs(v["primary_BeamA_500nM_48h_log2fc"]),
                                  v["gene_id"]))
        counts=Counter(x["priority_tier"] for x in output)
        sign_paths={}
        for t in ("1h","6h"):
            seen=[x for x in output if isinstance(x[t+"_sign_matches_external_48h"],bool)]
            sign_paths[t]={
                "all_96_present_with_finite_fc":len(seen),
                "same_direction_as_external_48h":sum(x[t+"_sign_matches_external_48h"] for x in seen),
                "opposite_direction_to_external_48h":sum(not x[t+"_sign_matches_external_48h"] for x in seen),
                "FDR_significant_at_timepoint":sum(
                    isinstance(x["external_TM_100nM_"+t+"_FDR"],float)
                    and x["external_TM_100nM_"+t+"_FDR"]<.05 for x in seen)
            }
        strong=[{k:v for k,v in item.items() if k not in ("mechanistic_novelty","causal_dose_reversal")}
                for item in output if item["priority_tier"]==grades[0]]
        result={
           "status":"PASS_EXPLORATORY_96_TM_DOSE_TIME_TIERING",
           "source_P5l_run":41158660082,
           "source_P5b_run":41085982185,
           "source_original_workbook":XLSX,
           "external_100nM_TM_times":["1h","6h","48h"],
           "n_discordant":len(output),
           "priority_counts":dict(counts),
           "strong_A_genes":[x["gene_symbol"] for x in strong],
           "strong_A_detailed":strong,
           "9x_formal_donor_LOO_coverage_in_96":len(overlap),
           "note_LOO":"P5b leave-one-pool-out was limited to original top100 genes: absence is NOT evidence of instability",
           "external_temporal_patterns":sign_paths,
           "source_zip_download_bytes":source_bytes,
           "interpretation":[
              "Tier thresholds were specified post hoc for hypothesis generation, not preregistered inference",
              "The public 100nM TM time-course is within one independent bulk study; not 500nM dose response",
              "The 500nM versus 100nM difference is confounded with cell state, donor source and protocol",
              "1h/6h/48h reported logFC not matched donor-by-donor from original raw reads",
              "Technical library-by-treatment confounding persists in GSE313736",
              "Clinical intraocular pressure effect, new regulatory mechanism and causal dose reversal are NOT established"]
        }
        with (out/"P5m_96_conflicts_timecourse_priorities.tsv").open("w",newline="",encoding="utf8") as f:
            writer=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
            writer.writeheader();writer.writerows(output)
        print("P5M_PASS",counts,"external_temporal",sign_paths,
              "full_9x_LOO_covered",len(overlap),flush=True)
    except Exception as e:
        result={"status":"HOLD_P5M_ERROR","error":type(e).__name__+": "+str(e)[:350]}
        print("P5M_ERROR",result["error"],flush=True)
    (out/"P5m_timecourse_tier_summary.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    if result["status"]!="PASS_EXPLORATORY_96_TM_DOSE_TIME_TIERING":
        raise SystemExit(2)
if __name__=="__main__":main()
