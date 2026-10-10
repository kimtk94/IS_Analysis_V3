#!/usr/bin/env python3
"""Raw-SNP validation of 30 numerically replayed IS coloc ABF assays; noncausal."""
import argparse,csv,hashlib,json,math,statistics
from collections import Counter
from pathlib import Path

def run(replay,source,out):
 with replay.open(newline="") as f:rows=[r for r in csv.DictReader(f,delimiter="\t") if r["status"]=="PASS"]
 if len(rows)!=30:raise ValueError("Expected exactly thirty PASS assays")
 detail=[]
 for assay in rows:
  path=source/assay["filename"]
  if not path.exists():raise FileNotFoundError(path)
  with path.open(newline="") as f:
   vals=list(csv.DictReader(f,delimiter="\t"))
  if len(vals)!=int(assay["actual_nsnps"]):raise ValueError("SNP input count drift")
  maf=[];qtlN=[];pvals=[]
  for x in vals:
   if x["harmonization"]!="EXACT_REF_ALT_GRCH38":raise ValueError("alleles not exact")
   g=float(x["gwas_maf"]);q=float(x["eqtl_maf"])
   n=float(x["eqtl_n"]);qp=float(x["eqtl_p"])
   if not all(map(math.isfinite,(g,q,n,qp))) or not (0<=g<=1 and 0<=q<=1 and n>0 and 0<=qp<=1):raise ValueError("invalid source data")
   maf.append(abs(g-q));qtlN.append(n);pvals.append(qp)
  large=sum(x>.1 for x in maf)
  if large!=int(assay["maf_diff_gt_0_1"]):raise ValueError("Source AF discordance not reproduced")
  if round(max(qtlN)-min(qtlN))!=int(assay["observed_qtl_n_range"]):raise ValueError("QTL-N range drift")
  detail.append(dict(gene_id=assay["gene_base"],dataset_key=assay["dataset_key"],locus=assay["locus"],
   source_nsnps=len(vals),prior_H4=float(assay["PP_H4"]),prior_H3=float(assay["PP_H3"]),
   AF_abs_gt_0_1_count=large,AF_abs_gt_0_1_fraction=large/len(vals),
   AF_abs_median=statistics.median(maf),AF_abs_max=max(maf),
   QTL_N_min=min(qtlN),QTL_N_max=max(qtlN),QTL_N_unique=len(set(qtlN)),
   QTL_N_first=qtlN[0],QTL_N_does_vary=len(set(qtlN))>1,
   nominal_qtl_p_lt_5e_8=sum(x<5e-8 for x in pvals),
   source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
   ABF_reproduced_numerically=True,
   source_AF_allele_convention_verified=False,
   scientific_gate="SOURCE_INPUT_RISK_REQUIRES_SENSITIVITY_NOT_CAUSAL"))
 out.mkdir(parents=True,exist_ok=True)
 detail.sort(key=lambda x:-x["prior_H4"])
 with (out/"IS_30_ABF_RAW_SNP_INPUT_RISKS.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(detail[0]),delimiter="\t");w.writeheader();w.writerows(detail)
 report=dict(status="PASS_30_ORIGINAL_SNP_INPUT_AUDIT_ONLY",assays=30,
   genes=len({x["gene_id"] for x in detail}),
   assays_AF_discordance_gt_25pct=sum(x["AF_abs_gt_0_1_fraction"]>.25 for x in detail),
   assays_AF_discordance_gt_50pct=sum(x["AF_abs_gt_0_1_fraction"]>.5 for x in detail),
   assays_variable_QTL_N=sum(x["QTL_N_does_vary"] for x in detail),
   max_original_PP_H4=max(x["prior_H4"] for x in detail),
   top5_assays=[{k:x[k] for k in ("gene_id","dataset_key","prior_H4","AF_abs_gt_0_1_fraction","QTL_N_unique","QTL_N_min","QTL_N_max")} for x in detail[:5]],
   no_new_valid_coloc=True,no_causal_gene_confirmed=True)
 (out/"IS_30_ABF_RAW_INPUT_RISK_SUMMARY.json").write_text(json.dumps(report,indent=2)+"\n")
 return report
if __name__=="__main__":
 ap=argparse.ArgumentParser();ap.add_argument("--replay",type=Path,required=True);ap.add_argument("--source",type=Path,required=True);ap.add_argument("--out",type=Path,required=True);a=ap.parse_args()
 print(json.dumps(run(a.replay,a.source,a.out),indent=2))
