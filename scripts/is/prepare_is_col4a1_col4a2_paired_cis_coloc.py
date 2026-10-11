#!/usr/bin/env python3
"""Read-only COL4A1 vs COL4A2 public 3-tissue full-cis harmonized inputs.

Only original eQTL Catalogue 3-study genomic rows, 2 original GRCh37 stroke
GWAS datasets and independent Fasta37 ancestry-unmatched reference. Retains
the same SNP universe between the two genes for matched comparisons.
No absent sample imputation, no causal inference and no single-signal claims.
"""
import argparse,csv,gzip,hashlib,json,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path

GENES={"COL4A1":"ENSG00000187498","COL4A2":"ENSG00000134871"}
DATASETS={"QTD000131":387,"QTD000171":205,"QTD000434":560}
SOURCE_COUNTS={
 ("COL4A1","QTD000131"):9025,("COL4A2","QTD000131"):9026,
 ("COL4A1","QTD000171"):8218,("COL4A2","QTD000171"):8218,
 ("COL4A1","QTD000434"):8763,("COL4A2","QTD000434"):8764
}
BOUND=("13",108801501,111807358)
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
 return h.hexdigest()
def writetsv(p,rows):
 if not rows:raise ValueError("Empty result")
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def snp_base_count_only(f):
 c=Counter()
 for line in f:
  s=line.rstrip("\r\n").split("\t")
  if len(s)!=19:raise ValueError("eQTL Catalogue original 19-col nominal QTL changed")
  if s[16] in GENES.values():
   c[s[16]]+=1
   yield s
def qtl_for_dataset(path,did,lo,fa,mod):
 seen={};ambiguous=set();counts=defaultdict(Counter);distinct=defaultdict(set)
 N_by_gene=defaultdict(list)
 qtl_variant_to_gene={}
 with gzip.open(path,"rt") as f:
  for c in snp_base_count_only(f):
   gene=c[16];gn=next(k for k,v in GENES.items() if v==gene)
   if c[1]!=BOUND[0] or not BOUND[1]<=int(c[2])<=BOUND[2]:
    raise ValueError("Wrong GRCh38 region source boundaries")
   if c[5]!=f"chr{c[1]}_{c[2]}_{c[3]}_{c[4]}":
    raise ValueError("Official source ID not REF_ALT in GRCh38")
   key38=(c[1],c[2],c[3].upper(),c[4].upper())
   distinct[gn].add(key38)
   if c[11]!="SNP":
    counts[gn]["NONSNP_excluded"]+=1;continue
   st,vid,strand=mod.canonical_qtl(c[1],int(c[2]),c[3].upper(),c[4].upper(),lo,fa,seen)
   counts[gn][st]+=1
   if st!="GRCH37_FASTA_REF_VERIFIED":continue
   x=[mod.number(c[i]) for i in (7,8,9,10,12,13)]
   if None in x:counts[gn]["MISSING_QTL_VALUES"]+=1;continue
   maf,pval,beta,se,ac,an=x
   if not (0<maf<=.5 and 0<=pval<=1 and se>0 and an>0 and 0<=ac<=an):
    counts[gn]["QTL_NUMERIC_INVALID"]+=1;continue
   if abs(min(ac/an,1-ac/an)-maf)>2e-4:counts[gn]["QTL_AC_AN_AF_CONFLICT"]+=1;continue
   N_by_gene[gn].append(an/2)
   alt=vid.split(":")[-1]
   item={
    "variant_id":vid,"source_GRCh38":c[5],
    "QTL_beta_ALT":beta,"QTL_se":se,"QTL_p":pval,
    "QTL_alt_AF":ac/an,"QTL_n_variant_AN_div_2":an/2,
    "QTL_liftover_strand":strand,
    "QTL_alt_is_official_source_effect_allele":True,
    "QTL_source_sig":tuple(c[i] for i in (7,8,9,10,12,13))}
   key=(gn,vid)
   if key in ambiguous:
    counts[gn]["GRCh38_to_GRCh37_collision_quarantined_other"]+=1;continue
   prev=qtl_variant_to_gene.get(key)
   if prev is not None:
    if prev!=item and prev["source_GRCh38"]!=item["source_GRCh38"]:
     ambiguous.add(key);qtl_variant_to_gene.pop(key)
     counts[gn]["GRCh38_to_GRCh37_collision"]+=1
    elif prev!=item and prev["source_GRCh38"]==item["source_GRCh38"]:
     if prev["QTL_source_sig"]!=item["QTL_source_sig"]:
      raise ValueError("Same-source SNP inconsistent duplicated rsid record")
     counts[gn]["duplicate_annotation_same_genomic_SNP"]+=1
    else:counts[gn]["duplicate_exact_record"]+=1
    continue
   qtl_variant_to_gene[key]=item
 for gn in GENES:
  found=len(distinct[gn]);expected=SOURCE_COUNTS[(gn,did)]
  if found!=expected:raise ValueError("Original cis-tested unique GRCh38 SNP/INDEL count mismatch "+gn+" "+did+":"+str((found,expected)))
  if not N_by_gene[gn]:raise ValueError("No derived N from original source")
  if statistics.median(N_by_gene[gn])!=DATASETS[did]:
   raise ValueError("Sample count inferred from QTL source AN/2 inconsistent "+gn+" "+did)
 return qtl_variant_to_gene,counts,N_by_gene
def gwas_index(path,trait,positions,fa,mod,helper):
 ret={};tot=Counter()
 for h in (f"chr{z}" for z in []):pass
 with gzip.open(path,"rt") as f:
  reader=csv.DictReader(f,delimiter="\t")
  for row in reader:
   tot["scanned_total_GWAS"]+=1
   chr_=row["chr"].removeprefix("chr")
   if chr_!="13":continue
   try:pos=int(row["pos"])
   except ValueError:tot["bad_position"]+=1;continue
   if pos not in positions:continue
   tot["GWAS_positions_in_6_sources"]+=1
   if trait=="EAS_AIS":state,v=helper.checked_gigastroke_fasta_anchored(row,fa,mod)
   else:state,v=mod.checked_gwas_row(row,trait,fa)
   tot[state]+=1
   if v is None:continue
   key=v["variant"]
   if key in ret and ret[key]!=v:raise ValueError("Duplicate GWAS exact biological allele mismatched")
   ret[key]=v
 if tot["scanned_total_GWAS"]<5000000:raise ValueError("Expected complete GWAS source")
 return ret,dict(tot)
def compare(qtls,both,bbj,ais,prior):
 pairs=[]
 output=[]
 for did in DATASETS:
  ga={v:x for (g,d,v),x in qtls.items() if g=="COL4A1" and d==did}
  gb={v:x for (g,d,v),x in qtls.items() if g=="COL4A2" and d==did}
  both_study={v for v in set(ga)&set(gb) if v in both}
  base={v for v in both_study if v in bbj and v in ais and frozenset(v.split(":")[-2:]) not in (frozenset(("A","T")),frozenset(("C","G")))}
  for gn,q in (("COL4A1",ga),("COL4A2",gb)):
   for v in sorted(base,key=lambda a:int(a.split(":")[1])):
    x=q[v];b=bbj[v];a=ais[v]
    out={
     "QTL_dataset":did,"gene_symbol":gn,"gene_id":GENES[gn],
     "GRCh37_fasta_REF_ALT_variant":v,
     "GRCh38_original_variant_id":x["source_GRCh38"],
     "QTL_beta_alt":x["QTL_beta_ALT"],"QTL_se":x["QTL_se"],
     "QTL_p":x["QTL_p"],"QTL_ALT_AF":x["QTL_alt_AF"],
     "QTL_n_variant_AN_div_2":x["QTL_n_variant_AN_div_2"],
     "BBJ_IS_beta_alt":b["beta_alt"],"BBJ_IS_se":b["se"],
     "BBJ_IS_p":b["p"],"BBJ_IS_ALT_AF":b["af_alt"],
     "EAS_AIS_beta_alt":a["beta_alt"],"EAS_AIS_se":a["se"],
     "EAS_AIS_p":a["p"],"EAS_AIS_ALT_AF":a["af_alt"],
     "EAS_AIS_hist_ID_order":a.get("historical_giga_id_order",""),
     "QTL_EAS_AIS_AF_DELTA":abs(x["QTL_alt_AF"]-a["af_alt"]),
     "QTL_BBJ_IS_AF_DELTA":abs(x["QTL_alt_AF"]-b["af_alt"]),
     "QTL_COL4A1_AND_COL4A2_IDENTICAL_SNP_TEST_UNIVERSE":True,
     "palindrome_excluded":True,
     "causal_gene_label":"NOT_ESTABLISHED",
    }
    output.append(out)
   pairs.append({
    "QTL_dataset":did,"gene_symbol":gn,
    "official_tested_GRCh38_variant_count":SOURCE_COUNTS[(gn,did)],
    "QTL_fasta_valid_in_GWAS_allele_indices":len(base),
    "QTL_strong_p_le_1e_4_on_paired_SNPs":sum(q[z]["QTL_p"]<=1e-4 for z in base),
    "EAS_AIS_GWS_p_lt_5e_8_on_paired_SNPs":sum(ais[z]["p"]<5e-8 for z in base),
    "both_marginal_sameSNP_QTL1e4_AIS5e8":sum(q[z]["QTL_p"]<=1e-4 and ais[z]["p"]<5e-8 for z in base),
    "median_QTL_AN_over_2":DATASETS[did],
    "paired_nonpal_SNP_count":len(base),
    "paired_nonpal_AF_delta_le_0p1":sum(abs(q[z]["QTL_alt_AF"]-ais[z]["af_alt"])<=.1 for z in base),
    "paired_nonpal_AF_delta_le_0p2":sum(abs(q[z]["QTL_alt_AF"]-ais[z]["af_alt"])<=.2 for z in base)
   })
 if len(pairs)!=6:raise ValueError("All three tissues and both COL4A genes required")
 for did in DATASETS:
  a=[x for x in pairs if x["QTL_dataset"]==did]
  if len(a)!=2 or a[0]["paired_nonpal_SNP_count"]!=a[1]["paired_nonpal_SNP_count"]:
   raise ValueError("Paired same-SNP QTL comparison missing")
  for g in GENES:
   k=(did,g);r=[x for x in output if (x["QTL_dataset"],x["gene_symbol"])==k]
   if len(r)!=a[0]["paired_nonpal_SNP_count"]:raise ValueError("Source gene matched count unexpected")
 return output,pairs
def main():
 ap=argparse.ArgumentParser()
 for k in ("qtl-root","bbj","ais","fasta","chain","vendor","prior-98-file","out-dir"):
  ap.add_argument("--"+k,type=Path,required=True)
 a=ap.parse_args()
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("Fail-closed existing output")
 sys.path.insert(0,str(Path(__file__).resolve().parent))
 sys.path.insert(0,str(a.vendor))
 import audit_is_eas5_public_qtl_gwas_snp_overlap as mod
 import audit_is_eas5_98_egene_fullcis_gwas_exact as helper
 from pyliftover import LiftOver
 mod.validate_reverse_chain(a.chain)
 with a.prior_98_file.open() as f:
  archive=[r for r in csv.DictReader(f,delimiter="\t") if r["region_id"]=="IS_XDATA_G0023" and r["gene_symbol"]=="COL4A2" and r["QTL_dataset"]=="QTD000434"]
 if len(archive)!=1 or int(archive[0]["NEW_fullcis_BOTH_nonpal_SNP_n"])!=3947:
  raise ValueError("Previous COL4A2 DLPFC EAS+BBJ source count changed")
 fa=mod.Fasta37(a.fasta)
 try:
  lo=LiftOver(str(a.chain))
  data={};nstatus={};inputhash={}
  for did in DATASETS:
   p=a.qtl_root/f"IS_XDATA_G0023__{did}__GRCh38_13_108801501_111807358.tsv.gz"
   if not p.is_file():raise FileNotFoundError(p)
   inputhash[str(p)]=sha(p)
   ds,st,n=qtl_for_dataset(p,did,lo,fa,mod)
   for (g,v),x in ds.items():data[(g,did,v)]=x
   nstatus[did]={
    "median_AN_over_2":{g:statistics.median(n[g]) for g in GENES},
    "per_gene_QTL_counts":{g:dict(st[g]) for g in GENES},
    "original_nominal_cis_SNP_maps_count":{g:len([k for k in ds if k[0]==g]) for g in GENES}
   }
   print("COL4A1_COL4A2_SOURCE_QTL",did,nstatus[did],flush=True)
  positions={int(v.split(":")[1]) for _,_,v in data}
  both={v for _,_,v in data}
  bbj,bbjs=gwas_index(a.bbj,"BBJ_IS",positions,fa,mod,helper)
  ais,aiss=gwas_index(a.ais,"EAS_AIS",positions,fa,mod,helper)
  if len(bbj)<5000 or len(ais)<5000:
   raise ValueError("Unexpectedly few validated chr13 AIS/BJJ source SNPs")
  rows,pairs=compare(data,both,bbj,ais,archive)
  if sum(p["QTL_fasta_valid_in_GWAS_allele_indices"] for p in pairs)<14000:
   raise ValueError("Insufficient both-gene matched SNPs")
  a.out_dir.mkdir(parents=True)
  writetsv(a.out_dir/"IS_COL4A1_COL4A2_TISSUE_PAIRED_COMMON_SNP_INPUT.tsv",rows)
  writetsv(a.out_dir/"IS_COL4A1_COL4A2_THREE_TISSUE_SNP_QC.tsv",pairs)
  summary={
    "status":"PAIRED_GENES_FULL_CIS_NONPALIN_GWAS_SNP_MATCH_ANALYTIC_INPUT_NOT_CAUSAL",
    "genes":GENES,"datasets_sample_size_from_QTL_AN2":DATASETS,
    "input_SNP_n":len(rows),"n_unique_matched_gene_tissue":6,
    "pairwise_equal_SNP_universe_by_QTL_study":True,
    "no_new_coloc_or_shared_causal_variant_claim":True,
    "GIGASTROKE_historical_ID_normalized_with_GRCh37_FASTA":True,
    "COL4A2_DLPFC_prior_98_source_match_count":3947,
    "before_antisense_or_gene_structure_annotation":True,
    "two_stroke_GWAS_cohorts_independence_not_established":True,
    "QTL_STUDY_LD_GENETIC_ANCESTRY_NOT_MATCHED":True,
    "source_QTL_median_AN_div_2_and_QC":nstatus,
    "gwas_scan_source": {"BBJ_IS":bbjs,"EAS_AIS":aiss},
    "input_sha256":{**inputhash,**{str(p):sha(p) for p in (a.bbj,a.ais,a.chain,a.prior_98_file)}}
  }
  (a.out_dir/"IS_COL4A1_COL4A2_THREE_TISSUE_INPUT_SCIENCE_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
  print("IS_COL4A1_COL4A2_PAIRED_CIS_GWAS_QC_PASS",json.dumps({
   "geneTissueRows":len(rows),"matched_SNPs_each_tissue":{k:p["paired_nonpal_SNP_count"] for p in pairs for k in [p["QTL_dataset"]]},
   "sameSNP_QTL_GWAS_joint":{f"{p['QTL_dataset']}_{p['gene_symbol']}":p['both_marginal_sameSNP_QTL1e4_AIS5e8'] for p in pairs}
  }),flush=True)
 finally:fa.close()
if __name__=="__main__":main()
