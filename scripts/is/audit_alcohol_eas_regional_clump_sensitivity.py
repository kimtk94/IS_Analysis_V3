#!/usr/bin/env python3
"""Evidence gates for Koyanagi Japanese alcohol GWAS 7 locus EAS clumping.

PLINK2 2.00a6: p1=5e-8, p2=1e-4, kb=1000, r2=.1;
EAS504 n, +/-500kb lead windows. Three MAF sensitivities 0/.01/.05.
Original p=0 underflow => 1e-300 ties; clump index selection unstable.
Clumps are LD-tagged *association sets*, not independently causal signals.
"""
import argparse,csv,json,math
from collections import Counter
from pathlib import Path
SRC=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_gws_regional_clumping")
REF=Path("/srv/is-analysis/data/is/ld_reference/alcohol_region_clump_eas_20261010")
REGIONS={"2":[("GCKR",27730940)],
         "4":[("KLB",39413780),("ADH1B",100239319)],
         "9":[("ALDH1B1",38395928),("ALDH1A1",75461066)],
         "12":[("chr12_106Mb",106750302),("ALDH2",112241766)]}
def read_tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def assign(c,pos):
    return min(((name,abs(int(pos)-target)) for name,target in REGIONS[c]),key=lambda x:x[1])
def file_for(c,maf,r2=.1):
    if maf==0:return f"chr{c}.EAS504.alcohol_r2_0p1.clumps"
    if r2==.1:return f"chr{c}.EAS504.alcohol_maf_{maf:.2f}_r2_0p1.clumps"
    return f"chr{c}.EAS504.alcohol_maf_0.01_r2_{r2}.clumps"
def run(root,ref):
    audit=json.loads((root/"IS_ALCOHOL_REGIONAL_GWAS_REFERENCE_COVERAGE.json").read_text())
    if audit["source_original_total_rows"]<7000000:
        raise ValueError("Incomplete original exposure scan")
    if not (audit["p1_index_threshold"]==5e-8 and audit["p2_secondary_threshold"]==1e-4):
        raise ValueError("Clumping source thresholds changed")
    matched=read_tsv(root/"IS_ALCOHOL_SOURCE_REF_MATCHED_REGIONAL_ASSOCIATIONS.tsv")
    source_index={x["ID"]:x for x in matched}
    if len(source_index)!=len(matched):
        raise ValueError("Duplicate association reference IDs")
    counts=[];leads=[];baseline=[]
    af={}
    for c in REGIONS:
        q=json.loads((ref/f"chr{c}.region_source_audit.json").read_text())
        if not q["validated"] or q["plink_sample_count"]!=504:
            raise ValueError("EAS reference not verified "+c)
        freqs=read_tsv(ref/f"chr{c}.EAS504.regions.afreq")
        af.update({x["ID"]:float(x["ALT_FREQS"]) for x in freqs if x["ALT_FREQS"]!="."})
    for maf in (0,.01,.05):
        for c in REGIONS:
            p=root/file_for(c,maf)
            if not p.exists():raise FileNotFoundError("Missing PLINK clump output "+str(p))
            entries=read_tsv(p)
            grouped=Counter();low=Counter()
            for e in entries:
                idx=e["ID"]
                if idx not in source_index:
                    raise ValueError("PLINK index not in original allele-matched GWAS")
                name,distance=assign(c,e["POS"])
                if distance>500000:raise ValueError("Clump index out of validated locus interval")
                grouped[name]+=1
                maf_ref=min(af[idx],1-af[idx])
                if maf and maf_ref+1e-9<maf:
                    raise ValueError("PLINK clump index violates reference MAF threshold")
                if maf_ref<.01:low["MAF_LT_1pct"]+=1
                elif maf_ref<.05:low["MAF_1_5pct"]+=1
                else:low["MAF_GE_5pct"]+=1
                a=source_index[idx]
                lead={"variant_id":idx,"chrom":c,"pos":int(e["POS"]),
                      "locus":name,"distance_from_initial_sentinel_bp":distance,
                      "MAF_threshold":maf,"r2_threshold":.1,"max_window_kb":1000,
                      "source_p_for_plink":e["P"],"source_p_original":a["P_original"],
                      "p_below_numeric_precision_clamped":a["P_clamped_from_zero"],
                      "reference_MAF_EAS504":maf_ref,
                      "clump_total_nearby_members_excl_index":int(e["TOTAL"]),
                      "reference_sample_size":504,
                      "index_is_unique_credible_set_or_independent_causal_variant":False,
                      "aldh2_protein_causal_effect_proven":False}
                leads.append(lead)
                if maf==0:baseline.append((e,a))
            for name,_ in REGIONS[c]:
                counts.append({"MAF_cutoff":maf,"r2_cutoff":.1,"chrom":c,
                    "locus":name,"clump_index_count":grouped[name],
                    "MAF_lt_0p01_index_count":low["MAF_LT_1pct"],
                    "MAF_0p01_to_0p05_index_count":low["MAF_1_5pct"],
                    "MAF_ge_0p05_index_count":low["MAF_GE_5pct"],
                    "source_GWS_in_chrom_after_harmonization":audit["chrom_source_stats"][c]["clump_input_GWS"],
                    "regional_not_genomewide":True})
    # r2 alternate threshold sensitivity is restricted to chr4/12, and should
    # not be added into whole-genome or whole-7-region totals.
    r2s={}
    for c in ("4","12"):
        for r2 in (.2,.5):
            path=root/file_for(c,.01,r2)
            if not path.exists():raise FileNotFoundError(path)
            rows=read_tsv(path)
            r2s[f"chr{c}_MAF1pct_r2_{r2}"]=len(rows)
    def num(maf):
        return sum(x["clump_index_count"] for x in counts if x["MAF_cutoff"]==maf)
    if num(0)<num(.01) or num(.01)<num(.05):
        raise ValueError("Unexpected increasing LD clumps when filtering rare SNPs")
    if sum(x["clump_index_count"] for x in counts if x["MAF_cutoff"]==0)!=len(baseline):
        raise ValueError("Baseline index count inconsistent")
    # Detailed sentinel: index or clump member; named by source ID, not
    # assuming lead SNP selected if study p values underflow to zero.
    assignments={}
    for c,pairs in REGIONS.items():
        all_entries=read_tsv(root/file_for(c,0))
        lookup={}
        for r in all_entries:
            lookup[r["ID"]]=r["ID"]
            members=[x for x in r["SP2"].split(",") if x not in ("","NONE",".")]
            for v in members:lookup.setdefault(v,r["ID"])
        for name,pos in pairs:
            src_rows=[x for x in matched if x["chr"]==c and int(x["pos"])==pos]
            if len(src_rows)!=1:raise ValueError("Original sentinel source unmatched")
            vid=src_rows[0]["ID"]
            assignments[name]={"original_GWAS_sentinel":vid,
                 "clump_index_assigned":lookup.get(vid,"NOT_IN_CLUMP_INDEX_OR_SECONDARY"),
                 "is_selected_clump_index":lookup.get(vid)==vid,
                 "source_original_p":src_rows[0]["P_original"],
                 "source_p_numeric_underflow":bool(int(src_rows[0]["P_clamped_from_zero"]))}
    with (root/"IS_ALCOHOL_7REGION_CLUMP_SENSITIVITY.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(counts[0]),delimiter="\t");w.writeheader();w.writerows(counts)
    with (root/"IS_ALCOHOL_7REGION_CLUMP_INDEX_EVIDENCE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(leads[0]),delimiter="\t");w.writeheader();w.writerows(leads)
    result={
      "reference":"1000G_EAS504_GRCh37_Phase3_v5b",
      "source_GWAS_Koyanagi_alcohol_v1_original_rows":audit["source_original_total_rows"],
      "source_p_lt_5e8_in_seven_windows":sum(x["source_p_lt_5e8"] for x in audit["chrom_source_stats"].values()),
      "allele_harmonized_GWS_in_seven_windows":sum(x["clump_input_GWS"] for x in audit["chrom_source_stats"].values()),
      "allele_harmonized_p_lt_1e4_regional":len(matched),
      "p_zero_was_clamped_in_matching":audit["p_zero_clamped_for_PLINK"],
      "p_zero_clump_index_tie_unstable":True,
      "clump_threshold":{"p1":5e-8,"p2":1e-4,"r2":.1,"radius_kb":1000},
      "regional_window_half_width_bp":500000,
      "initial_loci":7,
      "base_clump_count_MAF_unrestricted":num(0),
      "sensitivity_clump_count_MAF_ge_0p01":num(.01),
      "sensitivity_clump_count_MAF_ge_0p05":num(.05),
      "chr4_chr12_other_r2_thresholds":r2s,
      "sentinel_to_index_or_member":assignments,
      "one_clump_one_causal_signal_assumption_justified":False,
      "conditioned_causal_independence_verified":False,
      "EAS_region_reference_504_sample_ld_matched":True,
      "full_genome_wide_exposure_clumping_completed":False,
      "measured_individual_study_cohort_overlap":False,
      "horizontal_pleiotropy_excluded":False,
      "causal_eligible_MR_IVs":0,
      "MR_effects_calculated":False,
      "status":"REGIONAL_CLUMP_QC_PASS_INDEPENDENT_CAUSAL_SIGNALS_NOT_COUNTABLE"}
    (root/"IS_ALCOHOL_7REGION_CLUMP_SENSITIVITY_AUDIT.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=SRC)
    p.add_argument("--reference",type=Path,default=REF)
    a=p.parse_args();run(a.root,a.reference)
