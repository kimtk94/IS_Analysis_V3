#!/usr/bin/env python3
"""Merge two AIS GWAS lookups after independent GIGASTROKE raw+FASTA audit.

Correct four legacy REF/ALT-name reversals in *derived reporting only*.
Output is descriptive association evidence; no MR, independent replication,
causal inference, colocalization, or SuSiE reliability promotion.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def load_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def number(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):
        return None


def approximately_equal(a,b):
    x,y=number(a),number(b)
    return x is not None and y is not None and math.isclose(x,y,rel_tol=1e-9,abs_tol=1e-12)


def make_matrix(overlap_rows, provenance_rows, exposure_rows):
    prov={}
    for r in provenance_rows:
        if r["variant_id_fasta"] in prov: raise ValueError("Duplicate provenance variant")
        prov[r["variant_id_fasta"]]=r
    expo={}
    for r in exposure_rows:
        if r["ID"] in expo: raise ValueError("Duplicate exposure source variant")
        expo[r["ID"]]=r
    if len(prov)!=11: raise ValueError("Expected 11 provenance records")
    expected={(v,ds) for v in prov for ds in ("BBJ_JAPAN_IS","GIGASTROKE_EAS_AIS")}
    actual={(r["variant_id"],r["dataset"]) for r in overlap_rows}
    if len(overlap_rows)!=22 or actual!=expected: raise ValueError("Unexpected AIS cross-source matrix")
    result=[]
    for r in overlap_rows:
        vid=r["variant_id"]
        p=prov[vid]
        e=expo.get(vid)
        if e is None: raise ValueError("Exposure SNP not in source: "+vid)
        if e["reference_ALT"].upper()!=vid.split(":")[-1]:
            raise ValueError("Exposure ALT not aligned")
        ds=r["dataset"]
        if ds=="GIGASTROKE_EAS_AIS":
            if p["status"]=="MISSING_BOTH":
                if r["status"]!="MISSING_FROM_DATASET":raise ValueError("Missing GIGA source not missing in overlap")
                status="MISSING_FROM_DATASET"; beta=se=pvalue=eaf=""
            elif p["status"] in ("PASS_FASTA_REF_ALT","CANONICAL_ID_IS_OTHER_EFFECT_NOT_REF_ALT"):
                if p["source_canonical_stats_equal"]!="True" or p["source_md5_pass"]!="True":
                    raise ValueError("Source identity or raw-to-canonical statistic mismatch")
                if p["validated_fasta_variant_id"]!=vid:
                    raise ValueError("Reference corrected GIGA variant id conflict")
                if r["status"] not in ("MATCH_EXACT","REF_ALT_SWAP_REVIEW"):
                    raise ValueError("Unexpected GIGA status")
                if not (approximately_equal(r["beta_alt"],p["ais_beta_alt"])
                        and approximately_equal(r["se"],p["ais_se"])
                        and approximately_equal(r["p"],p["ais_p"])):
                    raise ValueError("Previous overlap stats conflict with raw GIGA")
                status=("VERIFIED_SOURCE_FASTA_ID_ALREADY_CORRECT" if p["status"]=="PASS_FASTA_REF_ALT"
                        else "VERIFIED_SOURCE_FASTA_ID_REPAIRED_IN_DERIVATIVE")
                beta,se,pvalue,eaf=p["ais_beta_alt"],p["ais_se"],p["ais_p"],p["ais_eaf_alt"]
            else:raise ValueError("Unreviewed GIGASTROKE status")
        elif ds=="BBJ_JAPAN_IS":
            if r["status"]=="MATCH_EXACT":
                status="VERIFIED_POSITION_REF_ALT_GWAS"
                beta,se,pvalue,eaf=r["beta_alt"],r["se"],r["p"],r["eaf_alt"]
            elif r["status"]=="MISSING_FROM_DATASET":
                status="MISSING_FROM_DATASET";beta=se=pvalue=eaf=""
            else:
                raise ValueError("Unreviewed BBJ status "+r["status"])
        else:raise ValueError("Unknown AIS source "+ds)

        pb=number(beta);pe=number(e["source_beta_ALT"])
        direction=("SAME_SIGN" if pb is not None and pe is not None and pb*pe>0 else
                   "OPPOSITE_SIGN" if pb is not None and pe is not None and pb*pe<0 else
                   "NO_DIRECTIONAL_EVIDENCE")
        result.append({
            "gene":r["locus"],"rsid":r["rsid_verified"],
            "variant_id_GRCh37_REF_ALT":vid,"AIS_source":ds,
            "source_quality_status":status,
            "alcohol_beta_ALT":e["source_beta_ALT"],
            "alcohol_se":e["source_se"],
            "alcohol_p_as_recorded":e["source_original_p"],
            "alcohol_p_is_numerical_zero":str(e["source_original_p"]) in ("0","0.0"),
            "ais_beta_ALT":beta,"ais_se":se,"ais_p":pvalue,"ais_eaf_ALT":eaf,
            "ais_source_original_id":p.get("original_canonical_id","") if ds=="GIGASTROKE_EAS_AIS" else r["source_variant_id"],
            "effect_direction_descriptive_only":direction,
            "independent_replication":"NOT_INDEPENDENT_BBJ_CONTRIBUTION_OR_RISK",
            "interpretation":"OBSERVATIONAL_ASSOCIATION_ONLY_NOT_MR"
        })
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--overlap",required=True,type=Path)
    p.add_argument("--giga-provenance",required=True,type=Path)
    p.add_argument("--alcohol-source",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    exposure=[]
    for gene in ("ADH1B","ALDH2"):
        exposure+=load_tsv(a.alcohol_source/gene/"variants.tsv")
    result=make_matrix(load_tsv(a.overlap),load_tsv(a.giga_provenance),exposure)
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Do not overwrite existing analysis outputs")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    dest=a.out_dir/"ALCOHOL_CS_AIS_VERIFIED_ASSOCIATION_MATRIX.tsv"
    with dest.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t")
        writer.writeheader();writer.writerows(result)
    info={"analysis":"ALCOHOL_CS_AIS_VERIFIED_ASSOCIATION_MATRIX",
          "status":"SOURCE_AUDITED_ASSOCIATION_ONLY",
          "utc":datetime.now(timezone.utc).isoformat(),
          "n_rows":len(result),"counts":dict(Counter(x["source_quality_status"] for x in result)),
          "gwas_coordinates":"Koyanagi_source_hg19_GRCh37_JCGE_official",
          "outcome_metadata":"GIGASTROKE_EAS_AIS_GRCh37",
          "BBJ_overlap":"POSSIBLE_DEFINITE_SHARED_COHORT_NOT_EXACT_PARTICIPANT_COUNT",
          "GIGASTROKE_includes_BBJ":"YES_NOT_INDEPENDENT_OF_BBJ_ISR_GWAS",
          "GIGASTROKE_source_metadata_md5":"904e12e33ec3c2eeefe4dc93bec4995c",
          "gwas_effect_orientation":"FASTA_VALIDATED_ALT",
          "LD_model_reliability":"NOT_RESOLVED",
          "ALDH2_fine_mapping":"BLOCKED",
          "ADH1B_fine_mapping":"EXPLORATORY",
          "MR":"NOT_PERFORMED",
          "coloc":"NOT_PERFORMED"}
    (a.out_dir/"ALCOHOL_CS_AIS_VERIFIED_MATRIX_MANIFEST.json").write_text(json.dumps(info,indent=2)+"\n")
    print("AIS_VERIFIED_MATRIX_COMPLETE",json.dumps(info["counts"]),flush=True)
    for x in result:
        if x["rsid"] in ("rs671","rs1229984"):
            print("TOP",x["AIS_source"],x["rsid"],x["source_quality_status"],
                 "alcohol_beta",x["alcohol_beta_ALT"],
                 "ais_beta",x["ais_beta_ALT"],"ais_p",x["ais_p"],flush=True)


if __name__=="__main__":main()
