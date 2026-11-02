#!/usr/bin/env python3
"""Third-cohort cross-study gene-sign audit: GSE301517 500 nM 7d TM (3 strains)
versus GSE313736 500 nM 48h paired 57-donor BeamA edgeR.

GEO normalized CPM is not raw counts: never fit edgeR or assign false FDR.
Per-strain technical samples are averaged WITHIN strain then 3 strains equally
weighted. No individual donor profiles, reads, private metadata or raw CPM exported.
"""
import argparse,csv,gzip,hashlib,io,json,math
from pathlib import Path
from statistics import mean
from ocular_gse301517_third_cohort_preflight import CPM,fetch

EXPECTED_SHA="717f4c9861de4d6b169439ef56a781f50d93ec918c456a9d0e61873e4a6e5f09"
SAMPLE_COLUMNS={
 "TM62":{"Control":[f"TM62_Control_{i}" for i in range(1,5)],
         "Dex":[f"TM62_Dex_{i}" for i in range(1,5)]},
 "TM70":{"Control":[f"TM70_Control_{i}" for i in range(1,3)],
         "Dex":[f"TM70_Dex_{i}" for i in range(1,3)]},
 "TM72":{"Control":["TM72_Control_1"],"Dex":["TM72_Dex_1"]},
}
PANEL=["MYOC","TGFBI","MMP2","ACTA2","SDC1","EGR3","LPAR1","TFPI",
       "LUM","IGFBP3","VCAM1","NTM","FKBP5","CHAC1","ERFE","PARM1"]
CORE=set(PANEL[:4]);CONTEXT=set(PANEL[4:8]);ECM=set(PANEL[8:12])
def readtsv(path,compressed=False):
    op=gzip.open if compressed else open
    with op(path,"rt",encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def value(s):
    try:
        x=float(s)
        return x if math.isfinite(x) and x>=0 else None
    except (TypeError,ValueError):return None
def out_tsv(p,rows):
    if not rows:raise ValueError("EMPTY_ROWS")
    with p.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--primary",required=True,type=Path)
    parser.add_argument("--emt",required=True,type=Path)
    parser.add_argument("--markers",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    primary=readtsv(args.primary,True)
    if len(primary)!=16050:raise ValueError("P5B_EXPECTED_16050_GENES")
    byid={}
    for row in primary:
        gid=row["gene_id"].split(".")[0]
        if gid in byid or not gid.startswith("ENSG"):raise ValueError("PRIMARY_ID_DUPLICATED_OR_INVALID")
        byid[gid]={"gene_symbol":row["gene_symbol"],"primary_logFC":float(row["logFC"]),
            "primary_FDR":float(row["FDR"]),"primary_logCPM":float(row["logCPM"])}
    emt=readtsv(args.emt)
    if len(emt)!=194:raise ValueError("P7A_EMT_EXPECTED_194")
    emt_ids={x["gene_id"].split(".")[0] for x in emt}
    emt_source_byid={x["gene_id"].split(".")[0]:x for x in emt}
    if len(emt_ids)!=194:raise ValueError("P7A_EMT_ENSG_DUPLICATED")
    cpm=fetch(CPM,26_000_000)
    actual_sha=hashlib.sha256(cpm).hexdigest()
    if actual_sha!=EXPECTED_SHA:raise ValueError("GSE301517_OFFICIAL_CPM_SHA256_CHANGED")
    rows=0;matched=0;non_numeric=0;allresults={};zero_genes=0
    with gzip.open(io.BytesIO(cpm),"rt",encoding="utf8",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        hdr=set(rd.fieldnames or [])
        need={"","external_gene_name","gene_biotype"}
        for d in SAMPLE_COLUMNS.values():
            for arm in d.values():need.update(arm)
        if not need.issubset(hdr):
            raise ValueError("GSE301517_EXACT_COLUMNS_NOT_AVAILABLE "+str(sorted(need-hdr)))
        for row in rd:
            rows+=1
            gid=(row.get("") or "").split(".")[0]
            if gid not in byid:continue
            matched+=1
            if gid in allresults:raise ValueError("DUP_MATCHED_GENE_ID")
            strain={}
            invalid=False
            for name,arms in SAMPLE_COLUMNS.items():
                measures={}
                for arm,colnames in arms.items():
                    nums=[value(row[k]) for k in colnames]
                    if any(x is None for x in nums):invalid=True;break
                    measures[arm]=mean(nums)
                if invalid:break
                ctl=measures["Control"];dex=measures["Dex"]
                lfc=math.log2((dex+1)/(ctl+1))
                strain[name]={"log2ratio":lfc,"cpm_control":ctl,"cpm_dex":dex}
            if invalid:
                non_numeric+=1;continue
            pooled_effect=mean(strain[s]["log2ratio"] for s in sorted(SAMPLE_COLUMNS))
            if all(strain[s]["cpm_control"]==0 and strain[s]["cpm_dex"]==0 for s in strain):
                zero_genes+=1
            info=byid[gid]
            effect=info["primary_logFC"]
            agreement=sum((strain[s]["log2ratio"]*effect>0) for s in strain)
            allresults[gid]={
              "gene_id":gid,"gene_symbol":info["gene_symbol"],
              "primary_500nM_48h_log2fc":effect,
              "primary_FDR":info["primary_FDR"],
              "primary_logCPM":info["primary_logCPM"],
              "third_500nM_7d_equal_strain_log2_CPM_ratio":pooled_effect,
              "third_vs_primary_sign_match":pooled_effect*effect>0,
              "three_strain_sign_match_n":agreement,
              "minimum_7d_3_strain_mean_CPM":min(mean(strain[s][arm] for s in strain) for arm in ("cpm_control","cpm_dex")),
              **{"third_"+s+"_log2_CPM_ratio":strain[s]["log2ratio"] for s in sorted(SAMPLE_COLUMNS)}
            }
    if rows!=85298 or matched<12000 or non_numeric>10:
        raise ValueError("SOURCE_ROWS_MATCHED_OR_NONNUMERIC_DRIFT "+str((rows,matched,non_numeric)))
    if not (CORE.issubset({x["gene_symbol"] for x in allresults.values()})):
        raise ValueError("CORE_GENE_NOT_ACCESSIBLE")
    selected=[]
    targetset=set(PANEL)
    for x in allresults.values():
        if x["gene_symbol"] in targetset:
            selected.append({"gene_symbol":x["gene_symbol"],"gene_id":x["gene_id"],
             "tier":"CORE" if x["gene_symbol"] in CORE else "CONTEXT" if x["gene_symbol"] in CONTEXT else "ECM_SECONDARY" if x["gene_symbol"] in ECM else "OTHER",
             "primary_500nM_48h_log2FC":round(x["primary_500nM_48h_log2fc"],5),
             "primary_48h_FDR":x["primary_FDR"],
             "third_500nM_7d_log2_CPM_ratio":round(x["third_500nM_7d_equal_strain_log2_CPM_ratio"],5),
             "third_TM62_log2_CPM_ratio":round(x["third_TM62_log2_CPM_ratio"],5),
             "third_TM70_log2_CPM_ratio":round(x["third_TM70_log2_CPM_ratio"],5),
             "third_TM72_log2_CPM_ratio":round(x["third_TM72_log2_CPM_ratio"],5),
             "third_vs_primary_sign_match":x["third_vs_primary_sign_match"],
             "three_strain_sign_match_n":x["three_strain_sign_match_n"],
             "nominal_third_7d_p_or_FDR":"NOT_ESTIMATED_NORMALIZED_CPM",
             "is_independent_clinical_phenotype":False})
    if len(selected)!=16 or len({x["gene_symbol"] for x in selected})!=16:
        raise ValueError("PANEL_16_UNIQUE_SYMBOLS_EXPECTED "+str(len(selected)))
    out_tsv(args.out/"GSE301517_16_p9_panel_3strain_vs_48h.tsv",
       sorted(selected,key=lambda r:(["CORE","ECM_SECONDARY","CONTEXT","OTHER"].index(r["tier"]),r["gene_symbol"])))
    common_emt=[x for k,x in allresults.items() if k in emt_ids]
    emt_rows=[{"gene_id":x["gene_id"],"gene_symbol":x["gene_symbol"],
        "primary_48h_log2FC":round(x["primary_500nM_48h_log2fc"],6),
        "third_7d_CPM_log2ratio":round(x["third_500nM_7d_equal_strain_log2_CPM_ratio"],6),
        "sign_match":x["third_vs_primary_sign_match"],
        "strains_3_matching_primary":x["three_strain_sign_match_n"]}
        for x in sorted(common_emt,key=lambda x:x["gene_symbol"])]
    out_tsv(args.out/"GSE301517_EMT_shared_3strain_direction.tsv",emt_rows)
    # The 100nM 48h external published effects are derived from the original
    # P7a paper supplement; they are NOT a third randomized dosing arm.
    threeway=[]
    for x in sorted(common_emt,key=lambda g:g["gene_symbol"]):
        ref=emt_source_byid[x["gene_id"]]
        external=float(ref["external_100nM_log2fc"])
        acute=x["primary_500nM_48h_log2fc"]
        chronic=x["third_500nM_7d_equal_strain_log2_CPM_ratio"]
        all_same=(acute*external>0 and external*chronic>0)
        threeway.append({"gene_symbol":x["gene_symbol"],"gene_id":x["gene_id"],
            "acute_500nM_48h_log2FC":round(acute,6),
            "external_100nM_48h_log2FC":round(external,6),
            "third_500nM_7d_CPM_log2ratio":round(chronic,6),
            "three_studies_all_same_direction":all_same,
            "external_48h_vs_7d_same_direction":external*chronic>0,
            "independent_three_way_dose_time_test":False})
    out_tsv(args.out/"GSE301517_EMT_194_3study_signed_direction.tsv",threeway)
    markers=readtsv(args.markers)
    external_markers={x["symbol"]:x for x in markers}
    if not CORE.issubset(external_markers):raise ValueError("P7_EXTERNAL_CORE4_GENE_MISSING")
    trio=[]
    for x in sorted(selected,key=lambda r:r["gene_symbol"]):
        sym=x["gene_symbol"]
        if sym not in CORE:continue
        original=external_markers[sym]
        ref=float(original["external_lfc"])
        acute=x["primary_500nM_48h_log2FC"]
        chronic=x["third_500nM_7d_log2_CPM_ratio"]
        trio.append({"gene_symbol":sym,"acute_500nM_48h_log2FC":acute,
            "published_100nM_48h_log2FC":ref,
            "chronic_500nM_7d_CPM_log2ratio":chronic,
            "all_three_same_direction":acute*ref>0 and ref*chronic>0,
            "three_strains_same_sign_as_acute":x["three_strain_sign_match_n"]==3,
            "not_a_causal_dose_or_time_test":True})
    out_tsv(args.out/"GSE301517_core4_3study_direction.tsv",trio)
    expressed=[x for x in allresults.values() if x["minimum_7d_3_strain_mean_CPM"]>=1]
    strong=[x for x in expressed if abs(x["third_500nM_7d_equal_strain_log2_CPM_ratio"])>=.25
              and abs(x["primary_500nM_48h_log2fc"])>=.25]
    summ={
      "status":"PASS_SOURCE_INTEGRITY_AND_DESCRIPTIVE_CROSS_STUDY_CPM",
      "source_study":"GSE301517","source_official_CPM_sha256":actual_sha,
      "source_official_rows":rows,"source_primary_edgeR_tested":len(primary),
      "matched_primary_Ensembl_ids":matched,"matched_rows_with_all_14_numeric_CPM":len(allresults),
      "n_missing_numeric":non_numeric,"all_arms_zero_genes":zero_genes,
      "GSE301517_n_treatment_control_runs":14,"GSE301517_strains":list(SAMPLE_COLUMNS),
      "GSE301517_technical_replicates_per_strain":{"TM62":4,"TM70":2,"TM72":1},
      "GSE301517_dose":"500 nM","GSE301517_exposure":"7 days",
      "GSE313736_dose":"500 nM","GSE313736_exposure":"48 hours",
      "not_a_matched_time_effect_experiment":True,
      "not_a_new_RNAseq_edgeR_model":True,
      "normalized_CPM_pseudocount":1,
      "independent_unit":"3 TM cultured strains; donor identities/independence not independently verified",
      "panel16_sign_matched_n":sum(x["third_vs_primary_sign_match"] for x in selected),
      "core4_sign_matched_n":sum(x["third_vs_primary_sign_match"] for x in selected if x["tier"]=="CORE"),
      "core4_sign_matched_3of3_n":sum(x["three_strain_sign_match_n"]==3 for x in selected if x["tier"]=="CORE"),
      "emt_matched_n":len(common_emt),
      "emt_sign_match_n":sum(x["third_vs_primary_sign_match"] for x in common_emt),
      "emt_all_three_studies_same_signed_n":sum(x["three_studies_all_same_direction"] for x in threeway),
      "emt_external_48h_vs_7d_same_signed_n":sum(x["external_48h_vs_7d_same_direction"] for x in threeway),
      "core4_all_three_studies_same_signed_n":sum(x["all_three_same_direction"] for x in trio),
      "gene_universe_expressed_both_7d_arms_meanCPM_ge_1_n":len(expressed),
      "gene_universe_both_studies_abs_effect_ge_point25_n":len(strong),
      "gene_universe_strong_signed_agree_n":sum(x["third_vs_primary_sign_match"] for x in strong),
      "reproducibility_not_clinical_causality":True,
      "cannot_infer_7d_vs_48h_causal_time_effect":True,
      "no_gene_level_p_values_in_third_study":True}
    (args.out/"GSE301517_3strain_vs_48h_evidence_gate.json").write_text(json.dumps(summ,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    print("THIRD_COHORT",json.dumps(summ,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
