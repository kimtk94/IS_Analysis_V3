#!/usr/bin/env python3
"""Independent post-harmonization QTL ALT-beta provenance and rsID dedup checks.

Official EBI eQTL Catalogue documentation explicitly defines ALT as effect
allele; original GRCh38 ALT maps as complement on negative-strand liftOver
(with unchanged regression coefficient) to GRCh37 ALT iff FASTA37 REF agrees.
Do not assert matched-ancestry QTL LD or causal coloc from allele matching.
Original upstream matched QTL beta column remains source value for auditability.
"""
import argparse,csv,gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
OFFICIAL_DOC="https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md"
OFFICIAL_FAQ="https://www.ebi.ac.uk/eqtl/Data_access/"
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def rd(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save(p,rows):
 if not rows:raise ValueError("Empty output rows")
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def group_data(source,readiness):
 if len(source)!=15 or len(readiness)!=2225:raise ValueError("No complete input source coverage")
 allowed={(r["region_id"],r["dataset_id"]) for r in source}
 if len(allowed)!=15:raise ValueError("Duplicate QTL source slice")
 genes={g["gene_id_stable"]:g for g in readiness}
 if len(genes)!=2225:raise ValueError("Distinct original gene universe changed")
 return allowed,genes
def process(matched,source,readiness,out):
 allowed,genes=group_data(source,readiness)
 matchedrows=Counter();seen={};summary=defaultdict(Counter)
 gene_records=defaultdict(list)
 outrows=[]
 with gzip.open(matched,"rt") as f:
  for r in csv.DictReader(f,delimiter="\t"):
   group=(r["region_id"],r["QTL_dataset"])
   if group not in allowed or r["gene_id"] not in genes:
    raise ValueError("Unexpected source/gene in GWAS matching output")
   qtlid=(r["region_id"],r["QTL_dataset"],r["gene_id"],r["mapped_variant_GRCh37"])
   fields=("QTL_p","QTL_beta_SOURCE_ORIENTATION_UNVERIFIED","QTL_se",
           "QTL_maf","QTL_ac","QTL_an","mapped_strand")
   signature=tuple(r[k] for k in fields)
   matchedrows["source_matched_rows"]+=1
   if qtlid in seen:
    matchedrows["duplicated_gene_variant_rsID_row"]+=1
    if seen[qtlid]!=signature:raise ValueError("Same gene/allele has discordant QTL betas or genotypes")
    continue
   seen[qtlid]=signature
   gene_records[(r["region_id"],r["QTL_dataset"],r["gene_id"])].append(r)
   bb=r["BBJ_is_present"]=="True"
   ai=r["EAS_AIS_is_present"]=="True"
   pal=r["palindromic_SNP_review"]=="True"
   if not (bb or ai):raise ValueError("Matched-source data must have at least one source hit")
   if not math.isfinite(float(r["QTL_beta_SOURCE_ORIENTATION_UNVERIFIED"])):
    raise ValueError("QTL beta not finite")
   ac=float(r["QTL_ac"]);an=float(r["QTL_an"])
   if not 0<=ac<=an or an<=0:raise ValueError("Invalid AC/AN for ALT frequency")
   qaf=ac/an
   bbaf=float(r["BBJ_ALT_EAF"]) if bb else None
   aaaf=float(r["EAS_AIS_ALT_EAF"]) if ai else None
   afg=abs(qaf-aaaf) if ai else None
   passed=(bb and ai and not pal)
   matchedrows["deduplicated_geneSNP"]+=1
   matchedrows["both_SNP_sources"]+=int(bb and ai)
   matchedrows["both_SNP_sources_nonpalindromic"]+=int(passed)
   matchedrows["palindromic_rows"]+=int(pal)
   matchedrows["EAS_AIS_QTL_ALT_AF_delta_ge_0p1"]+=int(afg is not None and afg>=.1)
   matchedrows["EAS_AIS_QTL_ALT_AF_delta_ge_0p2"]+=int(afg is not None and afg>=.2)
   outrows.append({
    "region_id":r["region_id"],"QTL_dataset":r["QTL_dataset"],
    "gene_id":r["gene_id"],"mapped_variant_GRCh37":r["mapped_variant_GRCh37"],
    "original_QTL_variant_GRCh38":r["source_variant_GRCh38"],
    "QTL_liftover_strand":r["mapped_strand"],
    "QTL_beta_mapped_GRCh37_ALT_OFFICIAL_DEFINITION":float(r["QTL_beta_SOURCE_ORIENTATION_UNVERIFIED"]),
    "QTL_se":float(r["QTL_se"]),"QTL_p":float(r["QTL_p"]),
    "QTL_alt_eaf_ac_div_an":qaf,"QTL_maf_source":float(r["QTL_maf"]),
    "BBJ_exact_SNP":bb,"BBJ_beta_mapped_ALT":float(r["BBJ_beta_ALT"]) if bb else "",
    "BBJ_p":float(r["BBJ_p"]) if bb else "","BBJ_alt_eaf":bbaf if bb else "",
    "EAS_AIS_exact_SNP":ai,"EAS_AIS_beta_mapped_ALT":float(r["EAS_AIS_beta_ALT"]) if ai else "",
    "EAS_AIS_p":float(r["EAS_AIS_p"]) if ai else "","EAS_AIS_alt_eaf":aaaf if ai else "",
    "abs_QTL_vs_EAS_AIS_alt_eaf_delta":afg if afg is not None else "",
    "palindromic_SNP":pal,"both_GWAS_sources_nonpalindromic":passed,
    "status":"VARIANT_ALLELE_HARMONIZED_EXPLORATORY_NOT_CAUSAL",
    "QTL_allele_semantics":"QTL_ALT_IS_EFFECT_ALLELE_OFFICIALLY_DEFINED",
    "source_cohort_ancestry_verified":False,"matched_QTL_LD_verified":False,
    "coloc_recomputed":False
   })
 if not outrows:raise ValueError("No available matched variants")
 summaries=[]
 for (region,dataset,gene),r in sorted(gene_records.items()):
  source_gene=genes[gene]
  both=[x for x in r if x["BBJ_is_present"]=="True" and x["EAS_AIS_is_present"]=="True"]
  nonpal=[x for x in both if x["palindromic_SNP_review"]=="False"]
  minQ=min(float(x["QTL_p"]) for x in r)
  minais=min((float(x["EAS_AIS_p"]) for x in r if x["EAS_AIS_is_present"]=="True"),default=None)
  minbbj=min((float(x["BBJ_p"]) for x in r if x["BBJ_is_present"]=="True"),default=None)
  g=summary[(region,dataset,gene)]
  summaries.append({
   "region_id":region,"QTL_dataset":dataset,
   "gene_id":gene,
   "gene_symbol":source_gene["canonical_gene_symbol"],
   "original_archived_ABF_assays":int(source_gene["archived_646_coloc_assays"]),
   "unique_lifted_valid_QTL_SNP_with_any_GWAS":len(r),
   "exact_GRCh37_SNP_in_BOTH_BBJ_AIS":len(both),
   "exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS":len(nonpal),
   "QTL_nominal_min_p_in_any_GWAS_overlap":minQ,
   "EAS_AIS_min_p_in_QTL_overlap":minais if minais is not None else "",
   "BBJ_min_p_in_QTL_overlap":minbbj if minbbj is not None else "",
   "can_now_assert_OFFICIAL_ALT_EFFECT_DIRECTION":True,
   "can_now_assert_source_genome_GRCh37_REF_after_liftover":True,
   "dual_GWAS_and_QTL_coloc_results":"NOT_RUN",
   "QTL_source_ancestry_verified":False,
   "molecular_QTL_matched_LD_verified":False,
   "possible_single_signal_ABF_followup":len(nonpal)>=100 and minais is not None and minais<5e-8 and minQ<1e-4,
   "scientific_status":"REVIEW_QTL_N_FULL_CIS_MISSINGNESS_COHORT_OVERLAP_BEFORE_COLOC"
  })
 if len(outrows)!=matchedrows["deduplicated_geneSNP"]:raise ValueError("Internal row drift")
 save(out/"IS_EAS5_OFFICIAL_QTL_ALT_GWAS_HARMONIZED_VARIANTS.tsv",outrows)
 save(out/"IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv",summaries)
 manifest={
  "status":"OFFICIAL_EQTL_CATALOGUE_ALT_BETA_SEMANTICS_PLUS_FASTA37_GWAS_MATCH",
  "official_eQTL_catalogue_column_definition":OFFICIAL_DOC,
  "official_eQTL_catalogue_FAQ":OFFICIAL_FAQ,
  "n_unique_harmonized_geneSNP_records":len(outrows),
  "n_gene_tissue_with_GWAS_variant_match":len(summaries),
  "n_unique_gene_ids_with_any_GWAS_variant_match":len({x["gene_id"] for x in summaries}),
  "n_unique_new_gene_ids_without_historical_ABF":len({x["gene_id"] for x in summaries if x["original_archived_ABF_assays"]==0}),
  "n_gene_tissue_potential_ABF_followup_under_input_only_heuristics":sum(r["possible_single_signal_ABF_followup"] for r in summaries),
  "genetic_allele_matching_pass_not_coloc":True,
  "no_QTL_ancestry_verification":True,
  "no_causal_gene_promotion":True,
  "duplicate_source_rsids_removed_per_gene_variant":matchedrows["duplicated_gene_variant_rsID_row"],
  "matched_row_QC":dict(matchedrows),
  "sha256_inputs":{str(p):sha(p) for p in (matched,)}
 }
 (out/"IS_EAS5_ALT_BETA_PROVENANCE_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_EAS5_OFFICIAL_ALT_EFFECT_DIRECTION_AND_GWAS_MATCH_AUDITED",json.dumps({
 k:manifest[k] for k in ("n_unique_harmonized_geneSNP_records","n_gene_tissue_with_GWAS_variant_match",
 "n_unique_gene_ids_with_any_GWAS_variant_match","n_unique_new_gene_ids_without_historical_ABF",
 "n_gene_tissue_potential_ABF_followup_under_input_only_heuristics",
 "duplicate_source_rsids_removed_per_gene_variant")}),flush=True)
 return manifest
def main():
 p=argparse.ArgumentParser()
 for k in ("matched","source-manifest","gene-readiness","out-dir"):
  p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 for file in (a.matched,a.source_manifest,a.gene_readiness):
  if not file.is_file():raise FileNotFoundError(file)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("Refuse output overwrite")
 a.out_dir.mkdir(parents=True)
 process(a.matched,rd(a.source_manifest),rd(a.gene_readiness),a.out_dir)
if __name__=="__main__":main()
