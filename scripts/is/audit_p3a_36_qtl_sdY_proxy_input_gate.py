#!/usr/bin/env python3
"""IS P3A: source-only empirical QTL SE/MAF/AN sample and sdY-route diagnostics.

Read the three original QTL Catalogue indexed whole reference union streams,
anchor their digests to the earlier complete original source CI, reproduce
exact SNP-by-gene AN/2 histograms independently, and describe an uncalibrated
per-SNP residual-variance proxy SE*sqrt(2*MAF*(1-MAF)*(AN/2)).

DO NOT confuse an empirical residual proxy with the direct phenotype sdY
required by the quantitative beta/varbeta coloc branch. Never output a
fabricated N_eff, sdY=1, prior, H4, or runnable coloc input. Publisher
GIGASTROKE beta/SE with type cc has a distinct effect-size route that need
not use dataset N/s as p-only inference would; scientific readiness still
depends on valid log-odds effect units, credible region and common SNP set.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import json
import math
import statistics
import time
from pathlib import Path

QTL=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
MATCH=importlib.import_module("scripts.is.audit_p3a_36_real_gwas_qtl_reciprocal_snp")
UNION=importlib.import_module("scripts.is.audit_p3a_official_cis_union_36")

STUDY_N={"QTD000131":387,"QTD000171":205,"QTD000434":560}
ANOMALY_DELTA=0.01
PROXY_COMMON_MAF_MIN=0.05
PROXY_COMMON_MAF_MAX=0.45
MAX_QTL_INDEXED_ROWS=1_200_000


def require(test,reason):
    if not test:raise ValueError(reason)


def pin_evidence(manifest,union,allele_n,orig_match,orig_cis,official):
    grp=manifest.get("groups",[])
    keys={(x["QTL_dataset"],x["gene_id"]) for x in grp}
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status")=="HOLD" and len(grp)==36
            and len(keys)==36 and len({x["gene_id"] for x in grp})==20,
            "genuine 36 original gene×tissue group source lost")
    require(union.get("schema")=="is-p3a-36-expanded-original-cis-reference-union-v1"
            and union.get("technical_status")==
                 "PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS"
            and union.get("scientific_status")=="HOLD"
            and union.get("n_approved_coloc_abf")==0
            and union.get("cis_union_reference_Ensembl105_GRCh38_1based")==
                UNION.EXPECTED_WINDOW and len(union.get("groups",[]))==36,
            "source validated 20-gene reference full cis union missing")
    require(allele_n.get("schema")=="is-p3a-primary-gwas-yaml-QTL-allele-N-36-v1"
            and allele_n.get("technical_status")==
                "PASS_PRIMARY_PUBLISHER_METADATA_AND_36_QTL_AN_SOURCE_READ"
            and len(allele_n.get("groups",[]))==36
            and allele_n.get("formal_coloc_abf_eligible")==0,
            "original exact SNP-specific AN/2 evidence not available")
    require(orig_match.get("schema")==
                "is-p3a-36-real-source-reciprocal-SNP-GWAS-harmonization-v1"
            and orig_match.get("technical_status")==
                 "PASS_QTL_GWAS_SOURCE_SNP_HARMONIZATION_AGGREGATES_ONLY"
            and len(orig_match.get("rows",[]))==36
            and orig_match.get("n_formal_coloc_abf_eligible")==0,
            "original two-ancestry signed-effect GWAS-QTL source baseline missing")
    require(orig_cis.get("schema")==
                "is-p3a-36-original-eQTL-permuted-n_variants-nominal-completeness-v1"
            and orig_cis.get("technical_status")==
                "PASS_3_ORIGINAL_PERMUTATION_AND_NOMINAL_REGION_SOURCES_READ"
            and orig_cis.get("n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT")==36
            and orig_cis.get("n_formal_coloc_abf_eligible")==0
            and len(orig_cis.get("groups",[]))==36,
            "independent publisher original full cis n_variants not proved")
    require(official.get("technical_status")==
                "PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE"
            and set(official["QTL_datasets"])==set(STUDY_N),
            "real original QTL public release source identity missing")
    for name,vals in (
        ("original-full-cis",union["groups"]),
        ("sample-AN",allele_n["groups"]),
        ("GWAS-allele",orig_match["rows"]),
        ("permuted-cis",orig_cis["groups"])):
        found={(x["QTL_dataset"],x["gene_id"]) for x in vals}
        require(found==keys and len(found)==36,name+" source group key divergence")
    for q,n in STUDY_N.items():
        require(official["QTL_datasets"][q]["nominal_catalogue_sample_size"]==n,
                "catalogue cohort N identity changed")
    return keys


def new_stats():
    return {"count":0,"n_even_positive":0,"n_invalid":0,
       "n_more_than_study":0,"n_hist":collections.Counter(),
       "qtl_maf_ac_an_delta_gt_0p01":0,
       "nominal_SE_proxy":[],
       "nominal_SE_proxy_common_MAF":[],
       "p_lt_5e_minus_8":0,"p_equal_zero":0}


def proxy_value(se,maf,n):
    """Simple-model descriptive residual SD proxy; NOT coloc sdY."""
    require(isinstance(n,int) and n>0 and 0<maf<1
            and math.isfinite(se) and se>0 and math.isfinite(maf),
            "invalid empirical proxy numeric input")
    val=se*math.sqrt(2*maf*(1-maf)*n)
    require(math.isfinite(val) and val>0,"invalid SE-variance proxy")
    return val


def quantile(vals,p):
    require(vals and 0<=p<=1,"quantile sample empty")
    ix=(len(vals)-1)*p
    low=int(math.floor(ix));high=int(math.ceil(ix))
    return vals[low]*(high-ix)+vals[high]*(ix-low) if low!=high else vals[low]


def summarize(vals):
    if not vals:return {"n":0,"median":None,"p10":None,"p90":None}
    x=sorted(vals)
    return {"n":len(x),"median":round(statistics.median(x),8),
            "p10":round(quantile(x,.10),8),"p90":round(quantile(x,.90),8)}


def scan_one(dataset,url,genes,open_tabix=None):
    if open_tabix is None:
        import pysam
        open_tabix=pysam.TabixFile
    region=UNION.EXPECTED_WINDOW
    stats={g:new_stats() for g in genes}
    seen={g:set() for g in genes}
    dig=hashlib.sha256(("\t".join(QTL.OFFICIAL_COLUMN_ORDER)+"\n").encode())
    total=0
    with open_tabix(url) as tabix:
        chrom="10" if "10" in tabix.contigs else "chr10" if "chr10" in tabix.contigs else None
        require(chrom is not None,"original QTL Tabix lacks chromosome 10")
        for line in tabix.fetch(chrom,region["start"]-1,region["end"]):
            if isinstance(line,bytes):line=line.decode("utf8")
            raw=line.rstrip("\r\n")
            dig.update((raw+"\n").encode("utf8"))
            total+=1
            require(total<=MAX_QTL_INDEXED_ROWS,"source indexed QTL scan exceeded cap")
            cols=raw.split("\t")
            require(len(cols)==len(QTL.OFFICIAL_COLUMN_ORDER),
                    "source QTL official 19 physical columns changed")
            fields=dict(zip(QTL.OFFICIAL_COLUMN_ORDER,cols))
            gene=fields["gene_id"]
            if gene not in stats:continue
            item=MATCH.qtl_record(fields)
            if item is None:continue
            key=(item["pos38"],item["ref"],item["alt"],item["molecular_trait_id"])
            if key in seen[gene]:continue
            seen[gene].add(key)
            x=stats[gene]
            x["count"]+=1
            require(item["molecular_trait_id"].split(".")[0]==gene,
                    "QTL per-SNP molecular phenotype ID not candidate gene")
            try:
                an=int(fields["an"]);ac=float(fields["ac"])
                good=(an>0 and an%2==0 and math.isfinite(ac)
                      and 0<=ac<=an)
            except (TypeError,ValueError,OverflowError):
                good=False
            if not good:
                x["n_invalid"]+=1
                continue
            n=an//2
            x["n_even_positive"]+=1
            x["n_hist"][n]+=1
            if n>STUDY_N[dataset]:x["n_more_than_study"]+=1
            observed_minor=min(ac/an,1-ac/an)
            if abs(observed_minor-item["maf"])>ANOMALY_DELTA:
                x["qtl_maf_ac_an_delta_gt_0p01"]+=1
            pv=proxy_value(item["se"],item["maf"],n)
            x["nominal_SE_proxy"].append(pv)
            if PROXY_COMMON_MAF_MIN<=item["maf"]<=PROXY_COMMON_MAF_MAX:
                x["nominal_SE_proxy_common_MAF"].append(pv)
            # This is diagnostic only, do not use p to filter full cis source.
            try:originalp=float(fields["pvalue"])
            except (TypeError,ValueError):continue
            if originalp<5e-8:x["p_lt_5e_minus_8"]+=1
            if originalp==0:x["p_equal_zero"]+=1
    result={}
    for gene,x in stats.items():
        hist=x.pop("n_hist")
        sample=x.pop("nominal_SE_proxy")
        common=x.pop("nominal_SE_proxy_common_MAF")
        size=sum(hist.values())
        result[gene]={
            "original_strict_nonpal_SNPs":x["count"],
            "positive_even_genotyped_AN_count":x["n_even_positive"],
            "invalid_AN_AC_count":x["n_invalid"],
            "AN_over_2_higher_than_catalogue_study_N":x["n_more_than_study"],
            "min_observed_SNP_N":min(hist) if hist else None,
            "max_observed_SNP_N":max(hist) if hist else None,
            "n_distinct_observed_SNP_N":len(hist),
            "n_records_at_dataset_nominal_N":hist[STUDY_N[dataset]],
            "fraction_SNP_records_at_dataset_nominal_N":
                round(hist[STUDY_N[dataset]]/size,7) if size else None,
            "most_common_N_observed":[{"variant_n":n,"records":v}
                for n,v in hist.most_common(8)],
            "maf_vs_source_ac_an_delta_gt_0p01":x["qtl_maf_ac_an_delta_gt_0p01"],
            "per_SNP_regression_SE_MAF_AN_residual_proxy_all":summarize(sample),
            "per_SNP_regression_SE_MAF_AN_residual_proxy_common_MAF":
                summarize(common),
            "p_lt_5e_minus_8_gene_assay_records":x["p_lt_5e_minus_8"],
            "p_equal_zero_gene_assay_records":x["p_equal_zero"],
            "quantitative_source_phenotype_sdY_independently_observed":False,
            "source_study_per_phenotype_model_regression_N_approved":False,
        }
    return result,{"original_indexed_source_url":url,
                   "original_full_region_rows_read":total,
                   "original_complete_region_stream_SHA256":dig.hexdigest(),
                   "source_cohort_study_n_not_regression_N":STUDY_N[dataset]}


def evaluate(manifest,union,allele_n,orig_match,orig_cis,official,bydataset,source):
    keys=pin_evidence(manifest,union,allele_n,orig_match,orig_cis,official)
    pna={(g["QTL_dataset"],g["gene_id"]):
         g["original_variant_allele_N_summary"] for g in allele_n["groups"]}
    pmatch={(g["QTL_dataset"],g["gene_id"]):g
            for g in orig_match["rows"]}
    pfull={(g["QTL_dataset"],g["gene_id"]):g
           for g in orig_cis["groups"]}
    rows=[]
    for g in manifest["groups"]:
        ds,gene=g["QTL_dataset"],g["gene_id"]
        key=(ds,gene)
        x=bydataset[ds][gene]
        old=pna[key]
        require(x["original_strict_nonpal_SNPs"]==
                    old["count_valid_nonpal_original_QTL"]
                and x["positive_even_genotyped_AN_count"]==
                    old["count_an_valid_even_autosomal"]
                and x["invalid_AN_AC_count"]==old["count_an_missing_or_invalid"]
                and x["min_observed_SNP_N"]==old["min_variant_an_over_2"]
                and x["max_observed_SNP_N"]==old["max_variant_an_over_2"]
                and x["n_distinct_observed_SNP_N"]==
                    old["n_distinct_variant_an_over_2"]
                and x["most_common_N_observed"]==
                    old["top_variant_an_over_2_count"]
                and x["maf_vs_source_ac_an_delta_gt_0p01"]==
                    old["count_QTL_maf_vs_minor_ac_an_delta_gt_0_01"],
                "original per-SNP sample N histogram/Maf independent check drift "+str(key))
        require(pfull[key]["published_permutation_vs_recovered_nominal_unique_agree"]
                and pfull[key]["original_permutation_rows_for_trait"]==1,
                "original full cis tested variants not exactly source-certified")
        matched=pmatch[key]
        require(matched["QTL_source_indexed_nonpal_valid"]==
                    x["original_strict_nonpal_SNPs"]
                and matched["QTL_uniquely_reciprocally_lifted_REF_valid"]==
                    x["original_strict_nonpal_SNPs"],
                "original QTL-GWAS SNP direct source baseline drift")
        ancestry={}
        for a in ("EAS","EUR"):
            orig=matched["per_ancestry"][a]
            nn=orig["n_nonpal_QTL_GWAS_recipient_allele_matches"]
            require(nn>0 and orig["dataset"]==
                  ("GCST90104545" if a=="EAS" else "GCST90104540"),
                  "not original actual EAS/EUR stroke signed GWAS source overlap")
            ancestry[a]={
                "original_signed_GWAS_QTL_nonpal_overlap_assay_rows":nn,
                "signed_beta_original_SE_field_technically_available":True,
                "coloc_case_control_beta_varbeta_route_does_not_require_invented_N_eff":
                    True,
                "case_control_log_odds_effect_units_verified_per_meta_analysis":False,
                "allele_freq_ancestry_and_full_cis_overlap_preapproved":False,
                "scientific_GWAS_abf_component_approved":False,
            }
        rows.append({"QTL_dataset":ds,"region_id":g["region_id"],
            "gene_id":gene,"gene_symbol":g["gene_symbol"],
            "QTL_source_variant_N_and_residual_proxy":x,
            "GWAS_ancestry_effect_size_route":ancestry,
            "QTL_inverse_normal_transformation_pipeline_documented":True,
            "this_exact_gene_QTL_source_sdY_measured":False,
            "may_set_sdY_equal_1_without_per_gene_verification":False,
            "QTL_exact_phenotype_regression_effective_N_source_verified":False,
            "QTL_quantitative_ABF_input_sdY_approved":False,
            "GWAS_original_AIS_N_eff_automatically_assigned":False,
            "case_fraction_used_in_effect_size_cc_route":False,
            "formal_coloc_abf_eligible":False,
            "formal_coloc_susie_eligible":False,
            "scientific_status":"HOLD_QTL_SD_Y_OR_N_MODEL_AND_GWAS_EFFECT_UNITS",
        })
    require(len(rows)==36 and len({x["gene_id"] for x in rows})==20,
            "36 real gene×tissue assays/20 genes missing")
    num=sum(x["QTL_source_variant_N_and_residual_proxy"]["original_strict_nonpal_SNPs"]
            for x in rows)
    require(num==148344,"original full nonpal source SNP assay denominator drift")
    varN=sum(x["QTL_source_variant_N_and_residual_proxy"]["n_distinct_observed_SNP_N"]>1
             for x in rows)
    return {
       "schema":"is-p3a-36-primary-qtlsdY-proxy-GWAS-case-control-effect-size-path-v1",
       "technical_status":"PASS_3_REAL_ORIGINAL_QTL_SD_Y_PROXY_INPUT_DIAGNOSTICS",
       "scientific_status":"HOLD",
       "n_real_gene_tissue_assays":36,"n_unique_genes":20,
       "n_independent_stroke_GWAS_loci":1,
       "n_strict_nonpal_source_gene_assay_SNP_records":num,
       "n_groups_with_variant_N_heterogeneity":varN,
       "n_groups_with_all_source_SNP_N_constant":36-varN,
       "source_original_QTL_pub_QTL_phenotype_inverse_normal_process_reference":
          "https://www.ebi.ac.uk/eqtl/Methods/",
       "coloc_case_control_beta_varbeta_documentation":
          "https://chr1swallace.github.io/coloc/articles/a02_data.html",
       "proxy_not_direct_sdY_and_not_approved_input":True,
       "QTL_source_direct_phenotype_sdY_approved":False,
       "GWAS_original_per_SNP_N_eff_fabricated":False,
       "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
       "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
       "source_provenance":source,"groups":rows,
    }


def main():
    ap=argparse.ArgumentParser()
    for f in ("manifest","union","allele-n","original-gwas-match",
              "original-cis-permuted","official","out"):
        ap.add_argument("--"+f,type=Path,required=True)
    a=ap.parse_args()
    result={"schema":"is-p3a-36-primary-qtlsdY-proxy-GWAS-case-control-effect-size-path-v1",
            "technical_status":"HOLD_SOURCE_NOT_READ",
            "scientific_status":"HOLD",
            "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
            "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0}
    try:
        def load(s):
            return json.loads(getattr(a,s.replace("-","_")).read_text(encoding="utf8"))
        manifest,union,allele_n,match,permuted,official=[
            load(s) for s in ("manifest","union","allele-n","original-gwas-match",
                              "original-cis-permuted","official")]
        keys=pin_evidence(manifest,union,allele_n,match,permuted,official)
        bydataset={};source={}
        for i,dataset in enumerate(STUDY_N):
            if i:time.sleep(8)
            genes={gene for q,gene in keys if q==dataset}
            url=official["QTL_datasets"][dataset]["source_https_probe"]
            groups,meta=scan_one(dataset,url,genes)
            prior=union["QTL_datasets"][dataset]
            require(meta["original_full_region_rows_read"]==
                        prior["interval_rows"]["region_rows_scanned"]
                    and meta["original_complete_region_stream_SHA256"]==
                        prior["region_stream_SHA256"],
                    "original QTL region stream SHA-256 differs from main verified source")
            bydataset[dataset]=groups;source[dataset]=meta
        result=evaluate(manifest,union,allele_n,match,permuted,official,
                        bydataset,source)
    except (ValueError,OSError,KeyError,TypeError,UnicodeError,TimeoutError,
            ImportError,OverflowError) as ex:
        result["technical_status"]="HOLD_PRIMARY_QTL_SOURCE_OR_PROXY_DIAGNOSTICS"
        result["error"]=str(ex)[:360]
        result["error_type"]=type(ex).__name__
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:result.get(k) for k in
           ("technical_status","scientific_status","n_real_gene_tissue_assays",
            "n_strict_nonpal_source_gene_assay_SNP_records",
            "n_groups_with_variant_N_heterogeneity",
            "n_groups_with_all_source_SNP_N_constant",
            "n_formal_coloc_abf_eligible","error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
