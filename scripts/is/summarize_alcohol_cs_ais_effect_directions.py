#!/usr/bin/env python3
"""Descriptive ALT-aligned Japanese alcohol versus AIS effect comparison.

No causal Wald-ratio/MR estimator is computed, and no promotion of loci.
ALDH2 alcohol p=0 is treated as numerical underflow, not literal zero probability.
"""
import argparse,csv,json,math
from collections import Counter
from pathlib import Path

def load(path):
    with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))

def run():
    p=argparse.ArgumentParser()
    p.add_argument("--ais-overlap",type=Path,required=True)
    p.add_argument("--alcohol-source",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    rows=load(a.ais_overlap)
    if len(rows)!=22: raise ValueError("Expect exactly 22 AIS results (11 variants x 2 datasets)")
    exposure={}
    for locus in ("ADH1B","ALDH2"):
        path=a.alcohol_source/locus/"variants.tsv"
        with path.open(newline="") as f:
            for r in csv.DictReader(f,delimiter="\t"):
                if r["ID"] in exposure: raise ValueError("Duplicate alcohol variant")
                exposure[r["ID"]]=r

    final=[]
    for r in rows:
        e=exposure.get(r["variant_id"])
        if not e: raise ValueError("Alcohol input not found: "+r["variant_id"])
        ref,alt=r["variant_id"].split(":")[2:]
        if e["reference_ALT"].upper()!=alt.upper(): raise ValueError("Exposure ALT conflict")
        beta_e=float(e["source_beta_ALT"])
        try: beta_o=float(r["beta_alt"])
        except (ValueError,TypeError): beta_o=None
        p_o=float(r["p"]) if r["p"] else None
        status=r["status"]
        if status not in ("MATCH_EXACT","REF_ALT_SWAP_REVIEW","MISSING_FROM_DATASET"):
            raise ValueError("Unreviewed overlap status: "+status)
        direction=("SAME_SIGN" if beta_o is not None and beta_e*beta_o>0 else
                   "OPPOSITE_SIGN" if beta_o is not None and beta_e*beta_o<0 else "NOT_ASSESSABLE")
        sig=("P_LT_5E_8" if p_o is not None and p_o<5e-8 else
             "P_GE_5E_8" if p_o is not None else "MISSING")
        final.append({
            "locus":r["locus"],"variant_id":r["variant_id"],"rsid":r["rsid_verified"],
            "dataset":r["dataset"],"AIS_overlap_status":status,
            "alcohol_trait":"Koyanagi2024_log2_grams_day_plus1",
            "alcohol_beta_ALT":beta_e,"alcohol_se":e["source_se"],
            "alcohol_source_p":e["source_original_p"],
            "alcohol_source_p_underflow":str(e["source_original_p"]) in ("0","0.0"),
            "ais_beta_ALT":r["beta_alt"],"ais_se":r["se"],"ais_p":r["p"],
            "ais_p_category":sig,"direction_descriptive_only":direction,
            "AIS_sample_n":r["sample_n"],"ALCOHOL_sample_n":e["source_n"],
            "quality_note":("GWAS_CANONICAL_REF_ALT_SWAP_REVIEW" if status=="REF_ALT_SWAP_REVIEW"
                else "NO_AIS_SUMMARY_STAT" if status=="MISSING_FROM_DATASET" else ""),
        })
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {a.out_dir}")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    dest=a.out_dir/"ALCOHOL_CS_AIS_ALIGNED_EFFECTS_EXPLORATORY.tsv"
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(final[0]),delimiter="\t");w.writeheader();w.writerows(final)
    manifest={"result_status":"DESCRIPTIVE_ASSOCIATION_ONLY_NOT_CAUSAL_INFERENCE",
        "comparison":"ALCOHOL_LOG2_GRAMS_DAY_PLUS1 vs BBJ_IS_AND_GIGASTROKE_EAS_AIS",
        "genome_build":"GRCh37_FASTA_CHECKED_COORDINATES_BUT_ORIGINAL_ALCOHOL_GWAS_BUILD_OFFICIALLY_UNRESOLVED",
        "effect_alignment":"INPUT_ALT","gwas_sample_overlap":"UNASSESSED","cohort_matched_LD":"UNAVAILABLE",
        "total_rows":len(final),"quality_status":dict(Counter(x["AIS_overlap_status"] for x in final))}
    (a.out_dir/"ALCOHOL_CS_AIS_EFFECT_DIRECTION_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("ALCOHOL_AIS_EFFECTS_COMPLETE",json.dumps(manifest["quality_status"]))
    for x in final:
        if x["rsid"] in ("rs671","rs1229984"):
            print("TOP",x["dataset"],x["rsid"],"alcohol_beta",x["alcohol_beta_ALT"],
                  "AIS_beta",x["ais_beta_ALT"],"AIS_p",x["ais_p"],"QC",x["AIS_overlap_status"])

if __name__=="__main__":run()
