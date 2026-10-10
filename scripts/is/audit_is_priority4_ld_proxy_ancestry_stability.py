#!/usr/bin/env python3
"""Noncausal ancestry-dependent GWAS proxy neighborhoods vs GTEx nominal eQTLs.

All variants remain in the universe. Does NOT recompute coloc PP.H4: coloc.abf
has no genotype LD matrix argument. No signal independence or allele effect
direction inference is made.
"""
import argparse,csv,hashlib,json,math,statistics
from collections import Counter
from pathlib import Path

LOCUS_GENES={
 "BBJ_IS_L001":(("FGF5","ENSG00000138675"),("C4orf22","ENSG00000197826")),
 "BBJ_IS_L002":(("CALHM2","ENSG00000138172"),("NEURL","ENSG00000107954")),
}
TISSUE="GTEx_V8__Brain_Cerebellar_Hemisphere"
R2_HIGH=0.8
R2_MATERIAL_DELTA=0.2
QTL_P_NOMINAL=1e-4

def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))

def number(v):
 try:
  x=float(v)
 except (ValueError,TypeError):
  return None
 return x if math.isfinite(x) else None

def classify(eas,eur):
 """Four-state proxy classification; missing LD is distinct from LD<0.8."""
 a=number(eas);b=number(eur)
 if a is None or b is None:return "UNINFORMATIVE_BOTH_ANCESTRIES"
 if not (0<=a<=1.0000001 and 0<=b<=1.0000001):
  raise ValueError("invalid r² outside 0..1")
 if a>=R2_HIGH and b>=R2_HIGH:return "HIGH_LD_BOTH"
 if a>=R2_HIGH:return "HIGH_LD_EAS_ONLY"
 if b>=R2_HIGH:return "HIGH_LD_EUR_ONLY"
 return "BELOW_HIGH_LD_BOTH"

def write(path,records):
 if not records:raise ValueError("empty audit records")
 with path.open("w",newline="") as f:
  writer=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
  writer.writeheader();writer.writerows(records)

def run(snp,source,out):
 data=read(snp)
 if len(data)!=3277:raise ValueError(f"Expected 3277 cross-ancestry SNP rows, found {len(data)}")
 loc={}
 for row in data:
  locus=row["locus"]
  if locus not in LOCUS_GENES:raise ValueError("unexpected locus")
  var=row["variant_b37"]
  loc.setdefault(locus,{})
  if var in loc[locus]:raise ValueError("duplicate reference SNP")
  loc[locus][var]=row
 if {k:len(v) for k,v in loc.items()}!={"BBJ_IS_L001":1694,"BBJ_IS_L002":1583}:
  raise ValueError("locus SNP denominator has changed")
 genes=[];panels=[];locus_details=[]
 for locus,markers in LOCUS_GENES.items():
  region=loc[locus]
  classes=Counter(classify(r["EAS_r2_with_GWAS_lead"],r["EUR_r2_with_GWAS_lead"]) for r in region.values())
  informative=[]
  reversals=0
  for row in region.values():
   ea=number(row["EAS_r2_with_GWAS_lead"]);eu=number(row["EUR_r2_with_GWAS_lead"])
   ra=number(row["EAS_signed_r_with_GWAS_lead"]);ru=number(row["EUR_signed_r_with_GWAS_lead"])
   if ea is None or eu is None:continue
   informative.append(abs(ea-eu))
   if ra is not None and ru is not None and abs(ra)>=.2 and abs(ru)>=.2 and ra*ru<0:reversals+=1
  if not informative:raise ValueError("no informative reference LD")
  panel=dict(locus=locus,region_snps=len(region),informative_both=len(informative),
    noninformative=len(region)-len(informative),
    high_ld_both=classes["HIGH_LD_BOTH"],high_ld_eas_only=classes["HIGH_LD_EAS_ONLY"],
    high_ld_eur_only=classes["HIGH_LD_EUR_ONLY"],
    below_high_ld_both=classes["BELOW_HIGH_LD_BOTH"],
    material_abs_delta_r2_gt_0p2=sum(x>R2_MATERIAL_DELTA for x in informative),
    median_abs_delta_r2=statistics.median(informative),
    mean_abs_delta_r2=statistics.mean(informative),
    max_abs_delta_r2=max(informative),
    confident_signed_r_direction_reversal=reversals,
    interpretation="REFERENCE_LD_TAGGING_ONLY_NOT_COLOC_POSTERIOR")
  panels.append(panel)
  for symbol,gene in markers:
   original_file=source/f"{locus}__{TISSUE}__{gene}.tsv"
   values=read(original_file)
   if len(values)!=len(region):raise ValueError("original input count mismatch")
   assay={}
   for r in values:
    k=r["variant_id"]
    if k in assay:raise ValueError("repeated source molecular trait SNP")
    if k not in region:raise ValueError("source SNP not reference aligned")
    p=number(r["eqtl_p"]);g=number(r["gwas_p"])
    if p is None or g is None or p<0 or p>1 or g<0 or g>1:
     raise ValueError("not a proper study association p")
    assay[k]=(p,g)
   if set(assay)!=set(region):raise ValueError("incomplete source overlap")
   qtl_top=sorted(assay,key=lambda snp:(assay[snp][0],snp))
   gwas_top=sorted(assay,key=lambda snp:(assay[snp][1],snp))
   lead=gwas_top[0]
   if any(r["lead_GWAS_SNP_b37"]!=lead for r in region.values()):
    raise ValueError("GWAS lead from source differs from LD table")
   lowp=[k for k in assay if assay[k][0]<QTL_P_NOMINAL]
   top20=qtl_top[:20]
   def counts(keys):
    cs=Counter(classify(region[k]["EAS_r2_with_GWAS_lead"],region[k]["EUR_r2_with_GWAS_lead"]) for k in keys)
    return cs
   cs=counts(lowp);t20=counts(top20)
   top=region[qtl_top[0]]
   a=number(top["EAS_r2_with_GWAS_lead"]);b=number(top["EUR_r2_with_GWAS_lead"])
   # Per-gene counts of nominal eQTL support in GWAS-lead LD neighborhoods,
   # NOT coloc posterior, qtl credible sets or independent QTL signals.
   genes.append(dict(locus=locus,symbol=symbol,gene_id=gene,tissue=TISSUE,
     n_original_snps=len(assay),gwas_lead_b37=lead,qtl_top_variant_b37=qtl_top[0],
     qtl_top_nominal_p=assay[qtl_top[0]][0],qtl_top_r2_eas=a,qtl_top_r2_eur=b,
     qtl_top_proxy_category=classify(top["EAS_r2_with_GWAS_lead"],top["EUR_r2_with_GWAS_lead"]),
     n_nominal_QTL_p_lt_1e4=len(lowp),
     nominal_QTL_high_ld_both=cs["HIGH_LD_BOTH"],
     nominal_QTL_high_ld_eas_only=cs["HIGH_LD_EAS_ONLY"],
     nominal_QTL_high_ld_eur_only=cs["HIGH_LD_EUR_ONLY"],
     nominal_QTL_below_high_ld_both=cs["BELOW_HIGH_LD_BOTH"],
     top20_QTL_high_ld_both=t20["HIGH_LD_BOTH"],
     top20_QTL_high_ld_eas_only=t20["HIGH_LD_EAS_ONLY"],
     top20_QTL_high_ld_eur_only=t20["HIGH_LD_EUR_ONLY"],
     top20_QTL_below_high_ld_both=t20["BELOW_HIGH_LD_BOTH"],
     coloc_abf_PP_H4_recomputed=False,valid_gtEx_in_study_LD=False,
     interpretation="EXPLORATORY_GWAS_PROXY_EQTL_OVERLAP_NOT_COLOCALIZATION"))
   for v in qtl_top[:20]:
    rr=region[v]
    locus_details.append(dict(locus=locus,symbol=symbol,
      rank_nominal_QTL_p=1+qtl_top.index(v),variant_b37=v,
      nominal_QTL_p=assay[v][0],GWAS_p=assay[v][1],
      EAS_reference_r2_to_GWAS_lead=rr["EAS_r2_with_GWAS_lead"],
      EUR_reference_r2_to_GWAS_lead=rr["EUR_r2_with_GWAS_lead"],
      proxy_category=classify(rr["EAS_r2_with_GWAS_lead"],rr["EUR_r2_with_GWAS_lead"]),
      status="NONCAUSAL_TAGGING_CONTEXT_ONLY"))
 out.mkdir(parents=True,exist_ok=True)
 write(out/"IS_PRIORITY4_LOCUS_EAS_EUR_PROXY_DISCORDANCE.tsv",panels)
 write(out/"IS_PRIORITY4_GENE_NOMINAL_EQTL_TAGGING.tsv",genes)
 write(out/"IS_PRIORITY4_TOP20_QTL_SNP_ANCESTRY_TAGGING.tsv",locus_details)
 j=dict(schema="IS_PRIORITY4_CROSS_ANCESTRY_TAGGING_V1",
   reference_snp_rows=len(data),unique_loci=len(panels),gene_tissue_assays=len(genes),
   r2_threshold_for_high_tag=R2_HIGH,qtl_p_threshold=QTL_P_NOMINAL,
   locus_summaries=panels,gene_summaries=genes,
   LD_matrix_directly_used_by_coloc_abf=False,
   PP_H4_recomputed=False,causal_gene_established=0,
   pvalue_trait_tested_universe_verified=False,
   original_candidate_gene_universe=2225,
   source_snp_table_sha256=hashlib.sha256(snp.read_bytes()).hexdigest(),
   warning="The original coloc.abf calculation uses SNP-level beta/SE, MAF and sample N but no LD matrix. Cross-ancestry reference r2 changes proxy relationships; it cannot numerically revise archived PP.H4, nor establish independent signals. Shared cis rows/tissue are an ascertained subset.")
 (out/"IS_PRIORITY4_LD_PROXY_STABILITY_SUMMARY.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("snp","source","out"):p.add_argument("--"+k,required=True,type=Path)
 a=p.parse_args()
 r=run(a.snp,a.source,a.out)
 for v in r["locus_summaries"]:print("LOCUS",v)
 for v in r["gene_summaries"]:print("GENE",v)
