#!/usr/bin/env python3
"""Read-only all-SNP GWAS/GTEx-v8 AF/allele provenance QC (646 coloc assays).

This explicitly DISAMBIGUATES:
  1. GWAS EAF vs GWAS mathematical MAF (min(EAF,1-EAF))
  2. GTEx input maf vs AC/AN and mathematical MAF
  3. ALLELE LABEL IDENTITY vs effect-direction provenance
  4. locus-level unique SNPs vs repeated SNP*gene*tissue records
  5. frequency discordance vs strand flip (NOT equivalent)
  6. palindromic site risk and uncertain effect-orientation semantics.

Official eQTL Catalogue columns: ac=alternative allele count; maf=minor allele
frequency, alt=effect allele. In this specific archived local subset, AC/AN
is <=0.5 for all source rows; that alone does NOT establish a minor allele
effect direction. Source catalog & original GWAS effect allele metadata must
be assessed separately.

NO filtering, no allele complementation, no change to existing coloc or
GWAS/QTL data. Outputs ONLY to new audit directory.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter,defaultdict
from pathlib import Path

REQUIRED = {
    "locus","dataset_key","gene_base","variant_id_hg38","match_key",
    "ref","alt","gwas_eaf","gwas_maf","eqtl_maf","ac","an",
    "gwas_beta","eqtl_beta","harmonization",
}
PAIRS={("A","T"),("T","A"),("C","G"),("G","C")}

def number(value):
    n=float(value)
    if not math.isfinite(n):raise ValueError("Nonfinite numeric field")
    return n

def classify_variant(r):
    parts=r["match_key"].split(":")
    if len(parts)!=4 or parts[1]=="":
        raise ValueError("Unexpected coordinate/REF/ALT variant key")
    ch,pos,ref,alt=parts
    if r["variant_id_hg38"]!=r["match_key"]:
        raise ValueError("Inconsistent GRCh38 SNP keys")
    if r["ref"]!=ref or r["alt"]!=alt:
        raise ValueError("Input REF/ALT inconsistent with match key")
    if r["harmonization"]!="EXACT_REF_ALT_GRCH38":
        raise ValueError("Original harmonization label not exact-match")
    if not ch.isdigit() or not pos.isdigit() or int(pos)<=0:
        raise ValueError("Invalid chromosome or position")
    if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT" or ref==alt:
        return {"type":"NON_SNV_OR_INVALID","palindromic":False}
    return {"type":"SNV","palindromic":(ref,alt) in PAIRS}

def audit_row(r):
    v=classify_variant(r)
    eaf=number(r["gwas_eaf"])
    gm=number(r["gwas_maf"])
    qtl=number(r["eqtl_maf"])
    ac=number(r["ac"])
    an=number(r["an"])
    if not (0<=eaf<=1 and 0<=gm<=0.5 and 0<=qtl<=0.5 and
            an>0 and 0<=ac<=an):
        raise ValueError("Bad AF/count range")
    ac_ratio=ac/an
    g_maf_from_eaf=min(eaf,1-eaf)
    q_maf_from_ac=min(ac_ratio,1-ac_ratio)
    if abs(gm-g_maf_from_eaf)>1e-9:
        raise ValueError("GWAS MAF derived from EAF disagrees")
    if abs(qtl-q_maf_from_ac)>1e-4:
        raise ValueError("QTL MAF differs from AC/AN minor AF >1e-4")
    return {
      "gwas_eaf":eaf,"gwas_maf":gm,"eqtl_maf":qtl,
      "qtl_ac_an_ratio":ac_ratio,
      "maf_abs_delta":abs(gm-qtl),
      "eaf_vs_AC_AN_delta":abs(eaf-ac_ratio),
      "eaf_complement_vs_AC_AN_delta":abs((1-eaf)-ac_ratio),
      "gwas_eaf_gt_half":eaf>0.5,
      "qtl_AC_AN_gt_half":ac_ratio>0.5,
      "palindromic":v["palindromic"],
      "palindromic_near_half":v["palindromic"] and
                             (abs(eaf-0.5)<=0.08 or abs(ac_ratio-0.5)<=0.08),
      "snv":v["type"]=="SNV",
      "site_id":r["locus"]+"|"+r["match_key"],
    }

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def write(path,rows,cols):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t")
        w.writeheader()
        for r in rows:w.writerow(r)

def summarize(input_dir,master):
    manifest={(r["locus"],r["dataset_key"],r["gene_base"]):r for r in master}
    if len(master)!=646 or len(manifest)!=646:
        raise ValueError("Must have exactly 646 original coloc records")
    per_assay=[]
    unique={}
    total=Counter()
    all_files=[]
    for key in sorted(manifest):
        locus,ds,gene=key
        name=f"{locus}__{ds}__{gene}.tsv"
        path=input_dir/name
        if not path.is_file():raise FileNotFoundError(path)
        all_files.append(path)
        expected=int(manifest[key]["nsnps"])
        n=0; counts=Counter();maf_deltas=[]
        with path.open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f,delimiter="\t")
            if not REQUIRED.issubset(reader.fieldnames or []):
                raise ValueError("Missing schema columns "+str(path))
            for r in reader:
                n+=1
                if (r["locus"],r["dataset_key"],r["gene_base"])!=key:
                    raise ValueError("Source locus/tissue/gene mismatch")
                x=audit_row(r)
                if not x["snv"]:raise ValueError("Non-SNV encountered; inspect before interpreting")
                site=x["site_id"]
                if site not in unique:
                    unique[site]={
                      "locus":locus,"match_key":r["match_key"],
                      "ref":r["ref"],"alt":r["alt"],
                      "gwas_eaf":x["gwas_eaf"],
                      "gwas_maf":x["gwas_maf"],
                      "palindromic":x["palindromic"],
                      "gwas_EAF_GT_0p5":x["gwas_eaf_gt_half"],
                      "n_assay_records":0,"n_tissues":set(),"qtl_AC_AN_min":1.,
                      "qtl_AC_AN_max":0.,
                      "qtl_maf_min":1.,"qtl_maf_max":0.,
                      "max_maf_abs_delta":0.,
                      "minor_maf_delta_gt_0p1_in_any_assay":False,
                      "QTL_AC_AN_GT_0p5_in_any_assay":False,
                    }
                v=unique[site]
                if abs(v["gwas_eaf"]-x["gwas_eaf"])>1e-8:
                    raise ValueError("Same locus SNP has different archived GWAS EAF")
                v["n_assay_records"]+=1
                v["n_tissues"].add(ds)
                v["qtl_AC_AN_min"]=min(v["qtl_AC_AN_min"],x["qtl_ac_an_ratio"])
                v["qtl_AC_AN_max"]=max(v["qtl_AC_AN_max"],x["qtl_ac_an_ratio"])
                v["qtl_maf_min"]=min(v["qtl_maf_min"],x["eqtl_maf"])
                v["qtl_maf_max"]=max(v["qtl_maf_max"],x["eqtl_maf"])
                v["max_maf_abs_delta"]=max(v["max_maf_abs_delta"],x["maf_abs_delta"])
                v["minor_maf_delta_gt_0p1_in_any_assay"]|=x["maf_abs_delta"]>0.1
                v["QTL_AC_AN_GT_0p5_in_any_assay"]|=x["qtl_AC_AN_gt_half"]
                counts["rows"]+=1
                counts["minor_maf_diff_gt_0p1"]+=x["maf_abs_delta"]>0.1
                counts["minor_maf_diff_gt_0p2"]+=x["maf_abs_delta"]>0.2
                counts["qtl_AC_AN_gt_half"]+=x["qtl_AC_AN_gt_half"]
                counts["gwas_EAF_gt_half"]+=x["gwas_eaf_gt_half"]
                counts["palindromic"]+=x["palindromic"]
                counts["palindromic_near_half"]+=x["palindromic_near_half"]
                counts["eaf_matches_qtl_AC_AN_better_than_complement"]+=(x["eaf_vs_AC_AN_delta"] < x["eaf_complement_vs_AC_AN_delta"])
                counts["complement_matches_AC_AN_better_than_eaf"]+=(x["eaf_complement_vs_AC_AN_delta"] < x["eaf_vs_AC_AN_delta"])
                maf_deltas.append(x["maf_abs_delta"])
        if n!=expected:raise ValueError("Expected historical SNP count changed")
        total.update(counts)
        per_assay.append({
           "locus":locus,"dataset_key":ds,"gene_base":gene,
           "nsnps":n,
           "fraction_MAF_difference_GT_0p1":counts["minor_maf_diff_gt_0p1"]/n,
           "fraction_MAF_difference_GT_0p2":counts["minor_maf_diff_gt_0p2"]/n,
           "median_MAF_absolute_difference":statistics.median(maf_deltas),
           "n_palindromic":counts["palindromic"],
           "n_palindromic_freq_near_half":counts["palindromic_near_half"],
           "n_qtl_AC_AN_GT_half":counts["qtl_AC_AN_gt_half"],
           "n_gwas_EAF_GT_half":counts["gwas_EAF_gt_half"],
           "n_EAF_closer_to_QTL_AC_AN":counts["eaf_matches_qtl_AC_AN_better_than_complement"],
           "n_complement_closer_to_QTL_AC_AN":counts["complement_matches_AC_AN_better_than_eaf"],
           "scientific_status":"FREQUENCY_PROVENANCE_REVIEW_NOT_STRAND_FLIP_INFERENCE",
        })
        if len(per_assay)%100==0:print("AF_AUDIT_PROGRESS",len(per_assay),"/646",flush=True)
    flat=[]
    for v in unique.values():
        z={**v,"n_distinct_GTEx_tissues":len(v["n_tissues"])}
        del z["n_tissues"]
        flat.append(z)
    flat.sort(key=lambda x:(x["locus"],x["match_key"]))
    report={
       "status":"AF_SOURCE_CONSISTENCY_DIAGNOSTIC_ONLY",
       "assays":len(per_assay),"total_SNP_by_assay_rows":total["rows"],
       "unique_locus_variant_keys":len(flat),
       "unique_variant_locus_counts":dict(Counter(x["locus"] for x in flat)),
       "assay_counts_MAF_diff_gt_0p1_over_half":sum(x["fraction_MAF_difference_GT_0p1"]>0.5 for x in per_assay),
       "assay_counts_MAF_diff_gt_0p1_over_quarter":sum(x["fraction_MAF_difference_GT_0p1"]>0.25 for x in per_assay),
       "total_rows_MAF_diff_gt_0p1":total["minor_maf_diff_gt_0p1"],
       "total_rows_MAF_diff_gt_0p2":total["minor_maf_diff_gt_0p2"],
       "total_rows_palindromic":total["palindromic"],
       "total_rows_palindromic_near_half":total["palindromic_near_half"],
       "unique_sites_palindromic":sum(x["palindromic"] for x in flat),
       "qtl_AC_AN_gt_half_rows":total["qtl_AC_AN_gt_half"],
       "sites_qtl_AC_AN_gt_half_in_any_assay":sum(x["QTL_AC_AN_GT_0p5_in_any_assay"] for x in flat),
       "rows_GWAS_EAF_gt_half":total["gwas_EAF_gt_half"],
       "rows_gwas_eaf_closer_to_qtl_ac_an":total["eaf_matches_qtl_AC_AN_better_than_complement"],
       "rows_gwas_eaf_complement_closer_to_qtl_ac_an":total["complement_matches_AC_AN_better_than_eaf"],
       "gwas_maf_computed_as_min_eaf_1_minus_eaf":"PASS_WITH_1e-9_TOLERANCE",
       "eqtl_maf_vs_min_ACAN_1minusACAN":"PASS_WITH_1e-4_TOLERANCE",
       "GTEx_ac_documented_as_ALT_count":"EQTL_CATALOGUE_DOCUMENTATION_SAYS_YES",
       "GTEx_alt_documented_as_effect_allele":"EQTL_CATALOGUE_AND_GTEX_FAQ",
       "GTEx_AC_AN_le_half_in_all_local_rows":"OBSERVED_NOT_PROOF_AC_MEANS_MINOR",
       "MAF_differences_are_NOT_automatic_allele_flips":True,
       "allele_effect_direction_cross_study_automatically_flipped":False,
       "original_QTL_and_GWAS_inputs_modified":False,
       "new_causal_gene_asserted":False,
       "interpretation":"ALT/REF reference ID matching is distinct from phenotype effect-allele convention, especially for palindromic SNVs; MAF cannot orient alleles.",
       "documentation":[
          "https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md",
          "https://www.ebi.ac.uk/eqtl/Data_access/",
          "https://www.gtexportal.org/home/faq",
       ],
    }
    return per_assay,flat,report,all_files

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input-dir",required=True,type=Path)
    p.add_argument("--master",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Refuse to overwrite prior audit")
    with a.master.open(newline="") as f:
        master=list(csv.DictReader(f,delimiter="\t"))
    assays,sites,report,paths=summarize(a.input_dir,master)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write(a.out_dir/"IS_646_COLOC_AF_PER_ASSAY.tsv",assays,list(assays[0]))
    write(a.out_dir/"IS_COLOC_AF_UNIQUE_LOCUS_SITES.tsv",sites,list(sites[0]))
    report["source_input_names_SHA256"]=hashlib.sha256("\n".join(sorted(p.name for p in paths)).encode()).hexdigest()
    report["archived_master_sha256"]=sha256(a.master)
    (a.out_dir/"IS_646_COLOC_ALLELE_AF_PROVENANCE.json").write_text(json.dumps(report,indent=2)+"\n")
    print("IS_646_ALLELE_AF_PROVENANCE_AUDIT_PASS",
          json.dumps({k:report[k] for k in (
           "assays","total_SNP_by_assay_rows","unique_locus_variant_keys",
           "qtl_AC_AN_gt_half_rows","rows_GWAS_EAF_gt_half",
           "assay_counts_MAF_diff_gt_0p1_over_half",
           "total_rows_palindromic")},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
