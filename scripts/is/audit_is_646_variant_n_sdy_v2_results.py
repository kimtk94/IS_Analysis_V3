#!/usr/bin/env python3
"""Independent V2 Colab output gate: 8956 native-GWAS SNPs + 646 sdY tests.

Refuses any missing manifest, divergent source hash, incomplete status,
contradictory posterior, or false native-liftover coverage claim. Never
promotes causality; evaluates computational + source provenance only.
"""
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
from statistics import median

KEY=("locus","dataset_key","gene_base")
ALLOWED_NATIVE="PASS_NATIVE_BBJ_ALLELE2_BETA_EAF_GRCH38_CHAIN"
NATIVE_CATALOG="IS_FULL_646_GWAS_VARIANT_NATIVE_STATUS.tsv"
STATUS="IS_VARIANT_N_SDY_CALIBRATION_STATUS_646.tsv"
INPUT_HASHES="IS_646_VARIANT_N_SDY_SOURCE_SHA256.tsv"
MANIFEST="IS_646_VARIANT_N_SDY_COLAB_MANIFEST.json"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def tsv(p):
    with p.open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def keys(rows,desc):
    mapped={}
    for row in rows:
        k=tuple(row[x] for x in KEY)
        if k in mapped:raise ValueError("DUPLICATE_"+desc+": "+str(k))
        mapped[k]=row
    return mapped

def finite(z):
    v=float(z)
    if not math.isfinite(v):raise ValueError("NONFINITE_NUMERIC_VALUE")
    return v

def check(inputs,g0_dir,reference,index,master,runner,out_dir):
    if out_dir.exists():raise FileExistsError("OUTPUT_AUDIT_EXISTS_NO_OVERWRITE")
    required=[MANIFEST,STATUS,INPUT_HASHES,NATIVE_CATALOG,
              "IS_VARIANT_N_SDY_CALIBRATION_SUMMARY.json",
              "G1_ORIGINAL_ALLELE_STRUCTURAL_V1/IS_G1_SOURCE_STRUCTURAL_STATUS_646.tsv",
              "G1_ORIGINAL_ALLELE_STRUCTURAL_V1/IS_G1_SOURCE_STRUCTURAL_SUMMARY.json"]
    for f in required:
        if not (inputs/f).is_file():raise FileNotFoundError("MISSING_V2_SOURCE:"+f)
    m=json.loads((inputs/MANIFEST).read_text())
    if m.get("schema")!="IS_PHASE10_11_646_VARIANT_N_SDY_COLAB_V2":
        raise ValueError("MANIFEST_SCHEMA_MUST_BE_V2")
    if m.get("status")!="FULL_646_VARIANT_N_SDY_CALIBRATION_COMPLETE":
        raise ValueError("COLAB_COMPLETE_STATUS_NOT_CONFIRMED")
    if m.get("counts")!={"PASS":646,"FAIL":0,"MISSING_INPUT":0}:
        raise ValueError("COLAB_RUN_INCOMPLETE")
    if m.get("inputs_available_before_run")!=646:
        raise ValueError("NOT_646_ORIGINAL_INPUTS")
    if m.get("R_package_coloc")!="5.2.3" or m.get("original_code_baseline_priors")!={"p1":1e-4,"p2":1e-4,"p12":1e-5}:
        raise ValueError("HISTORICAL_STATISTICAL_CONTRACT_CHANGED")
    if m.get("native_BBJ_8956_validation_status")!=ALLOWED_NATIVE or \
       m.get("native_BBJ_unique_variants_full_646")!=8956 or \
       m.get("native_BBJ_new_variants_audited")!=647 or \
       m.get("native_BBJ_30_previous_union")!=8309:
        raise ValueError("NATIVE_SNP_PROVENANCE_NOT_FULL")
    if m["source_sha256"]["original_index"]!=sha(index) or \
       m["source_sha256"]["original_master"]!=sha(master) or \
       m["source_sha256"]["embedded_R_script"]!=sha(runner) or \
       m["source_sha256"]["output_status"]!=sha(inputs/STATUS):
        raise ValueError("MANIFEST_CRYPTOSOURCE_FAILED")
    if m["native_BBJ_reference_fingerprint_SHA256"]!=sha(reference) or \
       m["native_BBJ_full_8956_status_SHA256"]!=sha(inputs/NATIVE_CATALOG):
        raise ValueError("NATIVE_CATALOG_SOURCE_HASH_CHANGED")
    all648=(tsv(inputs/STATUS),tsv(g0_dir/"IS_LEGACY_646_SNP_REPLAY_STATUS.tsv"),
       tsv(inputs/INPUT_HASHES),tsv(g0_dir/"IS_646_SOURCE_FILE_SHA256.tsv"),
       tsv(index),tsv(master),
       tsv(inputs/"G1_ORIGINAL_ALLELE_STRUCTURAL_V1/IS_G1_SOURCE_STRUCTURAL_STATUS_646.tsv"))
    if any(len(x)!=646 for x in all648):
        raise ValueError("DENOMINATOR_646_MISMATCH")
    rt,gt,nh,gh,idx,ms,g1=map(lambda pair:keys(*pair),
       zip(all648,["sdY","G0","new_hash","G0_hash","original_index","original_master","G1"]))
    if any(set(rt)!=set(other) for other in (gt,nh,gh,idx,ms,g1)):
        raise ValueError("ORIGINAL_646_KEY_UNIVERSE_DIVERGES")
    delta=[];cross05=[];cross075=[];an_fail=[]
    for k,row in rt.items():
        g=gt[k];n=nh[k];h=gh[k];nidx=idx[k];old=ms[k];z=g1[k]
        if row["status"]!="PASS" or g["status"]!="PASS" or \
           z["status"]!="PASS_STRUCTURAL_ONLY" or \
           n["status"]!="PASS" or h["status"]!="PASS":
            raise ValueError("NONPASS_ORIGINAL_INPUT:"+str(k))
        filename=row["source_filename"]
        if filename!=g["filename"] or filename!=n["filename"] or filename!=h["filename"] or \
           filename!=Path(nidx["file"]).name:
            raise ValueError("SOURCE_NAME_NOT_IDENTICAL:"+str(k))
        if n["source_sha256"]!=h["source_sha256"]:
            raise ValueError("G0_V2_SOURCE_SNP_BYTES_DIFFER:"+str(k))
        if int(row["n_snps"])!=int(g["actual_nsnps"]) or \
           int(row["n_snps"])!=int(nidx["nsnps"]) or \
           int(row["n_snps"])!=int(old["nsnps"]):
            raise ValueError("SNP_DENOMINATOR_NOT_HISTORICAL:"+str(k))
        if abs(finite(row["qtl_n_first"])-finite(nidx["qtl_n_first"]))>1e-9:
            raise ValueError("FIRST_QTL_N_NOT_HISTORICAL:"+str(k))
        for field,target in (("H3_original","PP_H3"),("H4_original","PP_H4")):
            if abs(finite(row[field])-finite(g[target]))>1e-8:
                raise ValueError("G0_REPRODUCIBILITY_CHANGED:"+field+":"+str(k))
        h4=finite(row["H4_original"]); v=finite(row["H4_variantN_sdY"]); d=finite(row["H4_delta"])
        if not 0<=v<=1 or not 0<=h4<=1 or abs((v-h4)-d)>1e-10:
            raise ValueError("INVALID_H4_DELTA:"+str(k))
        if finite(row["sdY_variant_n"])<=0 or finite(row["sdY_scalar_first"])<=0:
            raise ValueError("NONPOSITIVE_SDY:"+str(k))
        delta.append(d)
        if (h4>=.5)!=(v>=.5):cross05.append(k)
        if (h4>=.75)!=(v>=.75):cross075.append(k)
        if row["qtl_an_equals_2n_all"]!="TRUE":an_fail.append(k)
    if len(delta)!=646:raise ValueError("SOURCE_ROW_COUNT")
    reference_map={}
    with gzip.open(reference,"rt",newline="",encoding="utf8") as f:
        for z in csv.DictReader(f,delimiter="\t"):
            key=(z["locus"],z["variant_id"])
            if key in reference_map:raise ValueError("NATIVE_REFERENCE_DUPLICATE")
            reference_map[key]=z["native_variant_fingerprint_sha256_160bits"]
    if len(reference_map)!=56189:raise ValueError("NATIVE_REFERENCE_COUNT_INVALID")
    snps=tsv(inputs/NATIVE_CATALOG)
    if len(snps)!=8956:raise ValueError("NATIVE_VARIANTS_NOT_8956")
    unique=set(); old=0;new=0;per_locus={}
    for z in snps:
        key=(z["locus"],z["variant_id"])
        if key in unique:raise ValueError("NATIVE_SNP_DUPLICATE")
        unique.add(key)
        if z["native_BBJ_effect_chain_validation"]!=ALLOWED_NATIVE:
            raise ValueError("NATIVE_SNP_NOT_PASS:"+str(key))
        if reference_map.get(key)!=z["source_fingerprint_sha256_160bits"]:
            raise ValueError("NATIVE_SOURCE_HASH_MISMATCH:"+str(key))
        if z["original_subset"]=="CACHED30":old+=1
        elif z["original_subset"]=="NEW647":new+=1
        else:raise ValueError("NATIVE_SNP_SUBSET_LABEL_INVALID")
        per_locus[z["locus"]]=per_locus.get(z["locus"],0)+1
    if (len(unique),old,new)!=(8956,8309,647):
        raise ValueError("NATIVE_SOURCE_SET_PARTITION_INVALID")
    summary=json.loads((inputs/"IS_VARIANT_N_SDY_CALIBRATION_SUMMARY.json").read_text())
    if (summary.get("original_assays"),summary.get("pass"),summary.get("fail"),summary.get("missing"))!=(646,646,0,0):
        raise ValueError("R_SUMMARY_NOT_FULL_646_PASS")
    if summary.get("causal_evidence")!="NOT_ESTABLISHED":
        raise ValueError("CAUSE_STATUS_UNSUPPORTED")
    g1meta=json.loads((inputs/"G1_ORIGINAL_ALLELE_STRUCTURAL_V1/IS_G1_SOURCE_STRUCTURAL_SUMMARY.json").read_text())
    if g1meta.get("status_counts")!={"PASS_STRUCTURAL_ONLY":646}:
        raise ValueError("G1_SOURCE_646_NOT_FULL_PASS")
    verdict={
       "schema":"IS_646_VARIANT_N_SDY_V2_INDEPENDENT_SOURCE_PROOF_V1",
       "status":"COMPLETE_646_SDY_AND_8956_NATIVE_BBJ_PROVENANCE_COMPUTATIONAL",
       "n_tests":646,"n_G1_structural_pass":646,"n_native_GWAS_allele_lifted_pass":8956,
       "n_native_prior_cached":old,"n_native_newly_checked":new,
       "native_variants_by_locus":dict(sorted(per_locus.items())),
       "median_abs_delta_h4":median(abs(d) for d in delta),
       "max_abs_delta_h4":max(abs(d) for d in delta),
       "h4_ge0_5_old":sum(finite(x["H4_original"])>=.5 for x in rt.values()),
       "h4_ge0_5_new":sum(finite(x["H4_variantN_sdY"])>=.5 for x in rt.values()),
       "h4_threshold_0_5_crossing_keys":[list(k) for k in cross05],
       "h4_threshold_0_75_crossing_keys":[list(k) for k in cross075],
       "qtl_an_diff_2n_gene_tissue_count":len(an_fail),
       "source_status_sha256":sha(inputs/STATUS),
       "native_catalog_sha256":sha(inputs/NATIVE_CATALOG),
       "native_reference_SHA256":sha(reference),
       "causal_mechanism":"NOT_ESTABLISHED",
       "actual_GTEx_genotype_cohort_LD":"NOT_PROVIDED",
       "variant_level_SE_RECALCULATION":"NOT_PERFORMED",
       "original_four_loci_only":True,
    }
    out_dir.mkdir(parents=True)
    (out_dir/"IS_646_VARIANT_N_SDY_V2_INDEPENDENT_AUDIT.json").write_text(
       json.dumps(verdict,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(verdict,ensure_ascii=False))
    return verdict

if __name__=="__main__":
    a=argparse.ArgumentParser()
    for name in ("inputs","g0_dir","reference","index","master","runner","out_dir"):
        a.add_argument("--"+name,required=True,type=Path)
    x=a.parse_args()
    check(x.inputs,x.g0_dir,x.reference,x.index,x.master,x.runner,x.out_dir)
