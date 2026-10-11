#!/usr/bin/env python3
"""Exact GRCh38->GRCh37 QTL↔GWAS SNP overlap for all 98 EAS5 BH05 eGene assays.

This extends, rather than overwrites, previously audited five-region *narrow*
nominal QTL genomic slices. Each public molecular-QTL gene has a reconstructed
full cis-tested count consistent with EBI permutations, and the published
permutation lead genomic variant is present.

Crucial limitations:
- Official eQTL Catalogue ALT is the molecular-QTL beta effect allele.
- GIGASTROKE canonical variant_id is accepted as REF/ALT ONLY if independently
  checked against GRCh37 FASTA and its effect/other allele pair.
- Palindromic SNPs remain review only, not clean replicated matching.
- Count+lead inclusion do not prove all original tested SNP identities.
- GWAS groups are physically overlapping cohorts, not independent replication.
- Source EBI eQTL donor ancestry/LD are unverified; NO colocalization.
"""
import argparse,csv,gzip,hashlib,importlib.util,json,math,sys
from collections import Counter,defaultdict
from pathlib import Path

GENES=240
ASSAYS=98
EGENES=66
REGIONS={
 "IS_XDATA_G0015":("10",101403919,105460593,71),
 "IS_XDATA_G0023":("13",108801501,111807358,14),
 "IS_XDATA_G0007":("4",78736549,81763187,10),
 "IS_XDATA_G0025":("17",5691879,8725409,92),
 "IS_XDATA_G0018":("11",377434,3387434,53),
}
DATASETS={"QTD000131","QTD000171","QTD000434"}
STUDY=["BBJ_IS","EAS_AIS"]

def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for part in iter(lambda:f.read(4*1024*1024),b""):h.update(part)
 return h.hexdigest()
def save(path,rows):
 if not rows:raise ValueError("Cannot write empty scientific output")
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)

def validated_selection(readiness,master,previous):
 if len(readiness)!=720 or len(master)!=720 or len(previous)!=393:
  raise ValueError("Frozen eQTL 720 gene-tissue and narrow GWAS overlap 393 drift")
 fdr=[r for r in readiness if r["BH05_ebi_reprocessed_cis_eGene"]=="True"]
 if len(fdr)!=ASSAYS or len({r["gene_id"] for r in fdr})!=EGENES:
  raise ValueError("Expected 98 eQTL BH q<=.05 gene-tissue and 66 genes")
 keys={(r["region_id"],r["QTL_dataset"],r["gene_id"]) for r in fdr}
 frozen={(r["region_id"],r["QTL_dataset"],r["gene_id"]):r for r in master}
 older={(r["region_id"],r["QTL_dataset"],r["gene_id"]):r for r in previous}
 if len(frozen)!=720 or len(older)!=393 or not keys.issubset(set(frozen)&set(older)):
  raise ValueError("Frozen 393 GWAS old-window reference incomplete")
 for r in fdr:
  key=(r["region_id"],r["QTL_dataset"],r["gene_id"])
  source=frozen[key]
  if source["restored_full_cis_n_variants_count"]!="True" or \
     source["permuted_official_lead_variant_in_extended_QTL"]!="True":
   raise ValueError("QTL source not certified cis-count-plus-lead "+str(key))
  if int(r["official_n_cis_tested_genetic_variants"])!=int(source["official_permuted_n_variants"]):
   raise ValueError("Official EBI tested n_variants mismatch")
  if r["QTL_dataset"] not in DATASETS or r["region_id"] not in REGIONS:
   raise ValueError("Not a predeclared pilot gene/dataset")
 return fdr,older

def file_for(rid,did,a):
 chrom,lo,hi,_=REGIONS[rid]
 if rid=="IS_XDATA_G0015":root=a.qtl_g0015
 elif rid in ("IS_XDATA_G0023","IS_XDATA_G0007"):root=a.qtl_t1
 elif rid=="IS_XDATA_G0018":root=a.qtl_g0018
 else:root=a.qtl_t4
 return root/f"{rid}__{did}__GRCh38_{chrom}_{lo}_{hi}.tsv.gz"

def insert_uniquely_mapped_variant(storage,ambiguous,key,item,counts):
 if key in ambiguous:
  counts["GRCh38_TO_GRCh37_COLLISION_OTHER_SOURCE_SNP"]+=1
  return "QUARANTINED_AFTER_COLLISION"
 previous=storage.get(key)
 if previous is not None:
  if previous!=item:
   storage.pop(key)
   ambiguous.add(key)
   counts["GRCh38_TO_GRCh37_DIFFERENT_SOURCE_SNP_COLLISIONS"]+=1
   return "QUARANTINED_TWO_DISTINCT_SOURCE_SNP"
  counts["IDENTICAL_GRCh37_MAPPING_DUPLICATE"]+=1
  return "SAME_SOURCE_VARIANT_AGAIN"
 storage[key]=item
 return "UNIQUE_GRCH37_MAPPED_SNP"

def read_exact_qtl(fdr,a,fa,lo,mod):
 bysource=defaultdict(dict)
 for r in fdr:
  bysource[(r["region_id"],r["QTL_dataset"])][r["gene_id"]]=r
 variants={}
 expected_byassay=defaultdict(set)
 counts=defaultdict(Counter)
 inps={}
 lifted_cache={}
 ambiguous=set()
 for (rid,did),targets in sorted(bysource.items()):
  p=file_for(rid,did,a)
  if not p.is_file():raise FileNotFoundError(p)
  inps[str(p)]=sha(p)
  chrom,minpos,maxpos,_=REGIONS[rid]
  n_all=0
  with gzip.open(p,"rt") as f:
   for line in f:
    if not line.strip():continue
    n_all+=1
    c=line.rstrip("\r\n").split("\t")
    if len(c)!=19:raise ValueError("Original EBI nominal QTL 19-column schema drift")
    gene=c[16].split(".")[0]
    if gene not in targets:continue
    if c[1]!=chrom or not minpos<=int(c[2])<=maxpos or \
       c[5]!=f"chr{c[1]}_{c[2]}_{c[3]}_{c[4]}":
      raise ValueError("Selected QTL row outside expected verified GRCh38 source window")
    assay=(rid,did,gene)
    key=(c[1],int(c[2]),c[3].upper(),c[4].upper())
    # Deterministically deduplicate rsID annotation replicas.
    source_sig=tuple(c[i] for i in (7,8,9,10,12,13))
    if key in expected_byassay[assay]:
     cnt=counts[assay]
     cnt["source_duplicate_rsID_rows"]+=1
     prev=variants.get((assay,key))
     if prev and prev["source_sig"]!=source_sig:
      raise ValueError("Same genomic QTL allele pair with different beta/p/ac/an")
     continue
    expected_byassay[assay].add(key)
    cnt=counts[assay];cnt["source_distinct_genomic_alleles"]+=1
    if c[11]!="SNP":cnt["source_non_SNP"]+=1;continue
    stat,var37,strand=mod.canonical_qtl(c[1],int(c[2]),c[3].upper(),
      c[4].upper(),lo,fa,lifted_cache)
    cnt[stat]+=1
    if stat!="GRCH37_FASTA_REF_VERIFIED":continue
    vals=[mod.number(c[i]) for i in (7,8,9,10,12,13)]
    if None in vals:cnt["QTL_NONFINITE_STATS"]+=1;continue
    maf,pval,beta,se,ac,an=vals
    if not (0<=maf<=.5 and 0<=pval<=1 and se>0 and an>0 and 0<=ac<=an):
     cnt["QTL_NUMERIC_RANGE_INVALID"]+=1;continue
    if abs(min(ac/an,1-ac/an)-maf)>0.0002:
     cnt["QTL_AC_AN_MAF_DISCORDANT"]+=1;continue
    cnt["QTL_SNP_VERIFIED_NUMERIC"]+=1
    # This 'effect ALT' follows official eQTL catalogue column semantics.
    pal=frozenset(var37.split(":")[2:]) in (
       frozenset(("A","T")),frozenset(("C","G")))
    variantkey=(assay,var37)
    item={"source_sig":source_sig,"source_variant_38":c[5],
      "mapped_variant_37":var37,"QTL_p":pval,"QTL_beta_ALT":beta,
      "QTL_se":se,"QTL_AF_ALT":ac/an,"QTL_maf":maf,
      "strand":strand,"palindrome":pal}
    insert_uniquely_mapped_variant(variants,ambiguous,variantkey,item,cnt)
  for gene,r in targets.items():
   key=(rid,did,gene)
   n=int(r["official_n_cis_tested_genetic_variants"])
   if counts[key]["source_distinct_genomic_alleles"]!=n:
    raise ValueError("Source cis tested count changed "+str(key)+" "+str(counts[key]))
  print("QTL_98_EGENE_FULL_CIS_SLICE",rid,did,
    "allrows",n_all,"target",len(targets),"verifiedQTL",
    sum(counts[(rid,did,g)]["QTL_SNP_VERIFIED_NUMERIC"] for g in targets),flush=True)
 if len(counts)!=ASSAYS:raise ValueError("Not all 98 FDR gene-tissue target assays in raw data")
 print("QTL_GRCH38_GRCH37_NONINJECTIVE_LIFTOVER_QUARANTINED",
  len(ambiguous),"ambiguous_geneSNP_mapping_keys",flush=True)
 return variants,dict(inps),counts

def checked_gigastroke_fasta_anchored(row,fa,mod):
 """A documented historical other:effect pair is NOT necessarily REF:ALT."""
 from gigastroke_fasta_variant_contract import normalize_record,VariantContractError
 try:
  if (row.get("dataset")!="GCST90104545" or row.get("phenotype")!="AIS" or
      row.get("ancestry")!="EAS" or row.get("build")!="GRCh37"):
   return "GIGA_SOURCE_TRAIT_BUILD_SCOPE_INVALID",None
  ref=fa.get(row["chr"].replace("chr",""),int(row["pos"]))
  if ref is None:return "GIGA_FASTA_REF_UNAVAILABLE",None
  v=normalize_record(row,ref)
  se=mod.number(row["se"]);pv=mod.number(row["p"])
  if se is None or pv is None or se<=0 or not 0<=pv<=1 or v.eaf_alt is None:
   return "GIGA_EFFECT_SE_P_EAF_INVALID",None
  if v.historical_id_status not in ("FASTA_REF_ALT","OTHER_EFFECT_ORDER"):
   return "GIGA_HISTORICAL_ID_ORDER_UNVERIFIED",None
  pal=frozenset((v.ref,v.alt)) in (
   frozenset(("A","T")),frozenset(("C","G")))
  code=("FASTA_ANCHORED_HISTORICAL_OTHER_EFFECT_ORDER" if
         v.historical_id_status=="OTHER_EFFECT_ORDER" else
         "FASTA_ANCHORED_HISTORICAL_REF_ALT_ORDER")
  return code,{
   "variant":v.variant_id,"beta_alt":v.beta_alt,"se":se,
   "p":pv,"af_alt":v.eaf_alt,"palindrome":pal,
   "historical_giga_id_order":v.historical_id_status,
   "effect_orientation":v.effect_orientation
  }
 except (KeyError,TypeError,ValueError,VariantContractError):
  return "GIGA_FASTA_ALLELE_PAIR_OR_SOURCE_ID_MISMATCH",None

def scan_gwas(path,trait,target_positions,mod,fa,allow_fasta_recoded_giga=False):
 idx={}
 stat=Counter()
 with gzip.open(path,"rt") as f:
  header=f.readline().rstrip("\r\n").split("\t")
  positions={name:i for i,name in enumerate(header)}
  required={"chr","pos","variant_id","effect_allele","other_allele","beta","se","p","eaf"}
  if not required.issubset(positions):raise ValueError("GWAS source schema changed")
  chrom_id=positions["chr"];pos_id=positions["pos"]
  # Avoid 20m allocations of DictReader dictionaries: decode only candidate
  # target coordinates, then use previously frozen FASTA/effect allele contract.
  for line in f:
   if not line.strip():continue
   stat["source_rows"]+=1
   c=line.rstrip("\r\n").split("\t")
   ch=c[chrom_id].replace("chr","")
   if ch not in target_positions:continue
   try:pos=int(c[pos_id])
   except ValueError:stat["source_position_invalid"]+=1;continue
   if pos not in target_positions[ch]:continue
   stat["GWAS_at_QTL_mapped_position"]+=1
   row={name:c[i] for name,i in positions.items()}
   if trait=="EAS_AIS" and allow_fasta_recoded_giga:
    qc,item=checked_gigastroke_fasta_anchored(row,fa,mod)
   else:qc,item=mod.checked_gwas_row(row,trait,fa)
   stat[qc]+=1
   if item:
    key=item["variant"]
    if key in idx:
     if idx[key]!=item:raise ValueError("GWAS exact-variant duplicate with incompatible data")
     stat["identical_GWAS_source_duplicate"]+=1
    else:idx[key]=item
 if stat["source_rows"]<5000000:raise ValueError("GWAS source unexpectedly incomplete")
 print("FULLCIS_GWAS_SCANNED",trait,"n",stat["source_rows"],
   "atQTLpositions",stat["GWAS_at_QTL_mapped_position"],
   "ref_correct_and_valid",len(idx),flush=True)
 return idx,dict(stat)

def integrate(fdr,older,variant_map,BBJ,AIS,counts,a):
 rows=[];genestats={};variant_rows=[]
 by_assay=defaultdict(list)
 for (gene_assay,mapped_id),qrecord in variant_map.items():
  by_assay[gene_assay].append((mapped_id,qrecord))
 for r in sorted(fdr,key=lambda z:(z["region_id"],z["gene_id"],z["QTL_dataset"])):
  assay=(r["region_id"],r["QTL_dataset"],r["gene_id"])
  old=older[assay]
  qs=counts[assay];tot=Counter()
  targetkeys=by_assay.get(assay,[])
  for v,q in targetkeys:
   b=BBJ.get(v);a1=AIS.get(v)
   if b:tot["BBJ_exact"]+=1
   if a1:tot["AIS_exact"]+=1
   if not (b or a1):continue
   tot["either_GWAS"]+=1
   if b and a1:tot["both_GWAS"]+=1
   if q["palindrome"]:tot["palindromic_seen"]+=1
   else:
    tot["either_nonpal"]+=1
    if b and a1:tot["both_nonpal"]+=1
   dif=abs(q["QTL_AF_ALT"]-a1["af_alt"]) if a1 else None
   if dif is not None:
    if dif>.1:tot["QTL_AIS_ALT_AF_DELTA_GT_0P1"]+=1
    if dif>.2:tot["QTL_AIS_ALT_AF_DELTA_GT_0P2"]+=1
   if b and a1 and not q["palindrome"]:
    if dif is not None and dif<=.1:tot["both_nonpal_AF_delta_le_0p1"]+=1
    if q["QTL_p"]<=1e-4 and a1["p"]<5e-8:
     tot["same_nonpal_SNP_QTL1e4_and_AIS5e8"]+=1
   variant_rows.append({
    "region_id":assay[0],"QTL_dataset":assay[1],"gene_id":assay[2],
    "original_QTL_variant_GRCh38":q["source_variant_38"],
    "mapped_genomic_variant_GRCh37":v,"QTL_liftover_strand":q["strand"],
    "QTL_beta_GRCh37_ALT":q["QTL_beta_ALT"],
    "QTL_p":q["QTL_p"],"QTL_se":q["QTL_se"],
    "QTL_alt_AF":q["QTL_AF_ALT"],"palindrome_review":q["palindrome"],
    "BBJ_allele_exact_Fasta37":bool(b),"BBJ_beta_ALT":b["beta_alt"] if b else "",
    "BBJ_p":b["p"] if b else "","BBJ_alt_AF":b["af_alt"] if b else "",
    "EAS_AIS_allele_exact_Fasta37":bool(a1),
    "EAS_AIS_beta_ALT":a1["beta_alt"] if a1 else "",
    "EAS_AIS_p":a1["p"] if a1 else "",
    "EAS_AIS_alt_AF":a1["af_alt"] if a1 else "",
    "EAS_AIS_historical_variant_ID_order":a1.get("historical_giga_id_order","") if a1 else "",
    "EAS_AIS_allele_orientation":a1.get("effect_orientation","") if a1 else "",
    "QTL_vs_AIS_ALT_AF_delta":dif if dif is not None else "",
    "same_SNP_joint_marginal_threshold_not_coloc":(
      bool(b and a1 and not q["palindrome"] and
       q["QTL_p"]<=1e-4 and a1["p"]<5e-8)),
    "scientific_status":"ASSOCIATION_ONLY_NO_SHARED_CAUSAL_VARIANT"
   })
  if tot["either_GWAS"]==0:raise ValueError("Full cis positive eGene with no valid GWAS SNP overlap")
  if tot["both_nonpal"]<int(old["exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS"]):
   raise ValueError("Expanded full-cis failed to retain narrow GWAS exact nonpal SNP subset")
  if tot["both_GWAS"]<int(old["exact_GRCh37_SNP_in_BOTH_BBJ_AIS"]):
   raise ValueError("Expanded full-cis GWAS exact overlap count decreased")
  if tot["either_GWAS"]<int(old["unique_lifted_valid_QTL_SNP_with_any_GWAS"]):
   raise ValueError("Expanded full-cis lost previous narrow GWAS overlap")
  rows.append({
   "region_id":assay[0],"QTL_dataset":assay[1],"gene_id":assay[2],
   "gene_symbol":r["gene_symbol"],"archived_ABF_assays":r["archived_646_gene_tissue_ABF_count"],
   "official_permuted_p_beta":r["official_p_beta"],
   "official_studywide_BH_q":r["BH_qvalue_across_ALL_STUDY_tested_genes"],
   "official_full_cis_tested_n_variants":r["official_n_cis_tested_genetic_variants"],
   "QTL_source_unique_biallelic_SNP_liftover_verified_n":qs["QTL_SNP_VERIFIED_NUMERIC"],
   "QTL_source_non_SNP_n":qs["source_non_SNP"],
   "QTL_liftover_many_to_one_distinct_GRCh38_mapping_collisions":qs["GRCh38_TO_GRCh37_DIFFERENT_SOURCE_SNP_COLLISIONS"],
   "QTL_liftover_many_to_one_quarantined_extra_source_SNP":qs["GRCh38_TO_GRCh37_COLLISION_OTHER_SOURCE_SNP"],
   "QTL_source_grch37_ref_conflict_n":qs["GRCH37_FASTA_REF_DISAGREES"],
   "OLD_2MB_BBJ_AIS_exact_nonpal_n":old["exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS"],
   "NEW_fullcis_BBJ_exact_n":tot["BBJ_exact"],
   "NEW_fullcis_EAS_AIS_exact_n":tot["AIS_exact"],
   "NEW_fullcis_either_GWAS_exact_n":tot["either_GWAS"],
   "NEW_fullcis_BOTH_GWAS_exact_n":tot["both_GWAS"],
   "NEW_fullcis_BOTH_nonpal_SNP_n":tot["both_nonpal"],
   "NEW_fullcis_BOTH_nonpal_QTL_AIS_AF_delta_le_0p1_n":tot["both_nonpal_AF_delta_le_0p1"],
   "NEW_fullcis_same_nonpal_SNP_QTL1e4_AIS5e8_n":tot["same_nonpal_SNP_QTL1e4_and_AIS5e8"],
   "any_eGene_to_AIS_coloc_computed":False,
   "source_QTL_genotype_matched_LD_verified":False,
   "BBJ_and_GIGA_independent_stroke_replication":False,
   "scientific_status":"PUBLIC_EGENE_BH05_EXACT_GWAS_SNP_HARMONIZATION_ONLY"
  })
  print("FULLCIS_98_EGENE_GWAS",assay[0],r["gene_symbol"],assay[1],
        "oldBOTHnonpal",old["exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS"],
        "newBOTHnonpal",tot["both_nonpal"],
        "sameSNPjoint",tot["same_nonpal_SNP_QTL1e4_and_AIS5e8"],flush=True)
 if len(rows)!=ASSAYS:raise ValueError("98 studywide BH eGene assays changed")
 genes={}
 for r in rows:
  key=(r["region_id"],r["gene_id"])
  entry=genes.setdefault(key,{
   "region_id":r["region_id"],"gene_id":r["gene_id"],"gene_symbol":r["gene_symbol"],
   "archived_ABF_assays":int(r["archived_ABF_assays"]),
   "QTL_FDR05_studies":0,"fullcis_BBJ_AIS_both_nonpal_source_geneSNP":0,
   "fullcis_geneQTL_source_same_nonpal_marginal_joint_SNP":0,
   "new_H4_coloc_or_causal_genes":0
  })
  entry["QTL_FDR05_studies"]+=1
  entry["fullcis_BBJ_AIS_both_nonpal_source_geneSNP"]+=int(r["NEW_fullcis_BOTH_nonpal_SNP_n"])
  entry["fullcis_geneQTL_source_same_nonpal_marginal_joint_SNP"]+=int(r["NEW_fullcis_same_nonpal_SNP_QTL1e4_AIS5e8_n"])
 if len(genes)!=EGENES:raise ValueError("66 eGenes count changed")
 summary={
  "eGene_BH05_assays":ASSAYS,"distinct_source_eGenes":len(genes),
  "GWAS_SNP_match_total_record_unit_gene_tissue_SNP":sum(int(r["NEW_fullcis_BOTH_nonpal_SNP_n"]) for r in rows),
  "OLD_GWAS_SNP_match_record_unit_gene_tissue_SNP":sum(int(r["OLD_2MB_BBJ_AIS_exact_nonpal_n"]) for r in rows),
  "sameSNP_marginal_threshold_joint_assays":sum(int(r["NEW_fullcis_same_nonpal_SNP_QTL1e4_AIS5e8_n"])>0 for r in rows),
  "sameSNP_marginal_threshold_joint_distinct_genes":sum(int(r["fullcis_geneQTL_source_same_nonpal_marginal_joint_SNP"])>0 for r in genes.values()),
  "all_66_are_EBI_reprocessed_eGenes_NOT_stroke_causal_genes":True,
  "full_cis_variant_count_lead_ID_independent_study_LD_not_assessed":True
 }
 return rows,list(sorted(genes.values(),key=lambda x:(x["region_id"],x["gene_symbol"]))),variant_rows,summary

def main():
 p=argparse.ArgumentParser()
 for opt in ("fdr","master","older","bbj","ais","chain","fasta",
             "vendor","qtl-g0015","qtl-t1","qtl-t4","qtl-g0018","out-dir"):
  p.add_argument("--"+opt,type=Path,required=True)
 p.add_argument("--giga-fasta-anchored",action="store_true",default=False)
 p.add_argument("--giga-original-source",type=Path)
 p.add_argument("--giga-metadata",type=Path)
 a=p.parse_args()
 if a.giga_fasta_anchored:
  import hashlib
  if not a.giga_original_source or not a.giga_metadata:
   raise ValueError("FASTA-anchored GIGA requires exact raw source and metadata")
  if not a.giga_original_source.is_file() or not a.giga_metadata.is_file():
   raise FileNotFoundError("Required original GIGA evidence missing")
  digest=hashlib.md5()
  with a.giga_original_source.open("rb") as source:
   for part in iter(lambda:source.read(4*1024*1024),b""):digest.update(part)
  if digest.hexdigest().lower()!="904e12e33ec3c2eeefe4dc93bec4995c":
   raise ValueError("GIGASTROKE EAS AIS original raw source MD5 drift")
  metadata=a.giga_metadata.read_text()
  if ("genome_assembly: GRCh37" not in metadata or
      "sample_size: 256274" not in metadata or
      "is_harmonised: false" not in metadata or
      "data_file_md5sum: 904e12e33ec3c2eeefe4dc93bec4995c" not in metadata):
   raise ValueError("GIGASTROKE pre-GWAS-SSF source metadata drift")
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No overwrite")
 for fp in (a.fdr,a.master,a.older,a.bbj,a.ais,a.chain,a.fasta):
  if not fp.is_file():raise FileNotFoundError(fp)
 sys.path.insert(0,str(Path(__file__).resolve().parent))
 import audit_is_eas5_public_qtl_gwas_snp_overlap as mod
 mod.validate_reverse_chain(a.chain)
 sys.path.insert(0,str(a.vendor))
 from pyliftover import LiftOver
 fdr,old=validated_selection(read(a.fdr),read(a.master),read(a.older))
 fa=mod.Fasta37(a.fasta)
 try:
  print("FULLCIS_98_EGENE_START",len(fdr),len({r["gene_id"] for r in fdr}),flush=True)
  qtl,sourcechecks,counts=read_exact_qtl(fdr,a,fa,LiftOver(str(a.chain)),mod)
  if len(qtl)<50000:raise ValueError("Unexpectedly too few full-cis selected QTL SNPs")
  targetpos=defaultdict(set)
  for assay,variant in qtl:
   chrom,position,ref,alt=variant.split(":")
   targetpos[chrom].add(int(position))
  print("FULLCIS_QTL_ALL_INPUT_READY",len(qtl),
    "DISTINCT_MAPPED_GWAS_POSITIONS",sum(len(x) for x in targetpos.values()),flush=True)
  bbj,bstats=scan_gwas(a.bbj,"BBJ_IS",targetpos,mod,fa)
  ais,astats=scan_gwas(a.ais,"EAS_AIS",targetpos,mod,fa,
     allow_fasta_recoded_giga=a.giga_fasta_anchored)
  out,genes,variantrows,summary=integrate(fdr,old,qtl,bbj,ais,counts,a)
  a.out_dir.mkdir(parents=True,exist_ok=True)
  save(a.out_dir/"IS_EAS5_98_BH05_FDR_EGENE_FULLCIS_GWAS_ALLELE_READINESS.tsv",out)
  save(a.out_dir/"IS_EAS5_66_BH05_FDR_EGENE_FULLCIS_SNP_READINESS.tsv",genes)
  with gzip.open(a.out_dir/"IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz","wt",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(variantrows[0]),delimiter="\t")
   w.writeheader();w.writerows(variantrows)
  m={
   "status":"98_STUDYWIDE_EBI_EGENE_FULL_CIS_GRCH37_FASTA_ALLELE_EXACT_GWAS_OVERLAP_ONLY",
   **summary,"n_fullcis_QTL_verified_geneSNP_records":len(qtl),
   "GWAS_source_scans":{"BBJ_IS":bstats,"GIGASTROKE_EAS_AIS":astats},
   "raw_eQTL_source_sha256":sourcechecks,
   "sha256_additional_inputs":{str(fp):sha(fp) for fp in (
       a.fdr,a.master,a.older,a.chain,a.bbj,a.ais)},
   "any_colocalization_calculated":False,
   "EAS_QTL_genotype_LD_verified":False,
   "GIGASTROKE_EAS_BBJ_independent_replication":False,
   "GIGASTROKE_pre_GWAS_SSF_historical_ID_is_other_effect_order_in_some_SNPs":a.giga_fasta_anchored,
   "GIGASTROKE_original_MD5_and_metadata_verified":a.giga_fasta_anchored,
   "GIGASTROKE_alleles_for_recoded_variants_anchored_to_GRCh37_FASTA":a.giga_fasta_anchored,
   "genetic_causal_gene_promotion":False
  }
  (a.out_dir/"IS_EAS5_98_FDR_EGENE_FULLCIS_GWAS_ALLELE_MANIFEST.json").write_text(json.dumps(m,indent=2)+"\n")
  print("IS_EAS5_98_EGENE_FULL_CIS_BBJ_GIGA_AIS_MATCH_PASS",json.dumps({
    k:m[k] for k in ("n_fullcis_QTL_verified_geneSNP_records",
    "GWAS_SNP_match_total_record_unit_gene_tissue_SNP",
    "OLD_GWAS_SNP_match_record_unit_gene_tissue_SNP",
    "sameSNP_marginal_threshold_joint_assays",
    "sameSNP_marginal_threshold_joint_distinct_genes")}),flush=True)
 finally:fa.close()
if __name__=="__main__":main()
