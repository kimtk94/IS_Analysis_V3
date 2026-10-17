#!/usr/bin/env python3
"""Compare CALHM2, NEURL, FGF5 and C4orf22 on the same evidence criteria.
Descriptive variant counts only. Gene rows sharing a locus are not independent.
"""
import argparse,csv,json,hashlib
from pathlib import Path

GENES=("CALHM2","NEURL","FGF5","C4orf22")
CUTOFF=1e-4
def digest(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
 return h.hexdigest()
def evaluate(source,output):
 output=Path(output)
 if output.exists():raise FileExistsError("Immutable output must not exist")
 with Path(source).open(newline="") as f: data=list(csv.DictReader(f,delimiter="\t"))
 if len(data)!=20689:raise ValueError("Frozen QTL/GWAS input count drift")
 result=[];variants=[]
 for gene in GENES:
  raw=[r for r in data if r["symbol"]==gene]
  if not raw:raise ValueError("Missing comparator "+gene)
  good=[r for r in raw if r["BBJ_exact_Fasta37"]=="True" and
        r["GIGA_EAS_AIS_exact_Fasta37"]=="True" and
        r["palindromic_review"]=="False"]
  unique={}
  for r in good:
   k=r["SNP_GRCh37"]
   if k in unique and unique[k]!=r:raise ValueError("Conflicting QTL row "+k)
   unique[k]=r
  if len(unique)<1000:raise ValueError("Unexpected loss of matched variants")
  q=[r for r in unique.values() if float(r["QTD000161_nominal_p"])<CUTOFF]
  both=[r for r in unique.values() if float(r["BBJ_p"])<5e-8 and float(r["GIGA_EAS_AIS_p"])<5e-8]
  common=[r for r in q if float(r["BBJ_p"])<5e-8 and float(r["GIGA_EAS_AIS_p"])<5e-8]
  direction=sum(float(r["BBJ_beta_ALT"])*float(r["GIGA_EAS_AIS_beta_ALT"])>0 for r in unique.values())
  result.append(dict(gene=gene,locus=raw[0]["locus"],
    molecular_assay_rows=len(raw),both_GWAS_allele_exact_nonpal_unique=len(unique),
    QTD_nominal_p_lt_1e4=len(q),BBJ_and_GIGA_GWS_same_SNP=len(both),
    triple_QTD_BBJ_GIGA_threshold=len(common),
    BBJ_GIGA_same_direction_fraction=direction/len(unique),
    independent_stroke_cohort_replication_established=False,
    causal_gene_evidence_established=False))
  for r in sorted(q,key=lambda r:int(r["SNP_GRCh37"].split(":")[1])):
   variants.append(dict(gene=gene,variant_GRCh37=r["SNP_GRCh37"],
      variant_GRCh38=r["SNP_GRCh38"],QTD_p=r["QTD000161_nominal_p"],
      QTD_beta=r["QTD000161_beta"],QTD_se=r["QTD000161_se"],
      BBJ_p=r["BBJ_p"],BBJ_beta_ALT=r["BBJ_beta_ALT"],
      GIGA_AIS_EAS_p=r["GIGA_EAS_AIS_p"],GIGA_beta_ALT=r["GIGA_EAS_AIS_beta_ALT"],
      two_GWAS_effect_directions_match=float(r["BBJ_beta_ALT"])*float(r["GIGA_EAS_AIS_beta_ALT"])>0))
 if {x["gene"]:x["QTD_nominal_p_lt_1e4"] for x in result}!={"CALHM2":3,"NEURL":0,"FGF5":5,"C4orf22":0}:
  raise ValueError("Expected molecular signal count differs from established source audit")
 output.mkdir(parents=True)
 for name,rows in (("IS_QTL_GENE_COMPARATOR_SUMMARY.tsv",result),
                   ("IS_QTL_GENE_COMPARATOR_SIGNIFICANT_VARIANTS.tsv",variants)):
  if not rows:raise ValueError("No comparator variant rows")
  with (output/name).open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
   w.writeheader();w.writerows(rows)
 meta=dict(status="SOURCE_MARGINAL_COMPARATOR_ONLY",source=str(source),
  source_sha256=digest(source),source_rows=len(data),qtl_cutoff=CUTOFF,
  qtl_sample_N_unverified=True,not_colocalization=True,
  correlated_variants_not_independent=True,original_2225_IS_genes_unchanged=True)
 (output/"IS_QTL_GENE_COMPARATOR_MANIFEST.json").write_text(json.dumps(meta,indent=2)+"\n")
 return result
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--source",required=True,type=Path)
 p.add_argument("--out",required=True,type=Path)
 a=p.parse_args()
 print(json.dumps(evaluate(a.source,a.out),indent=2))
