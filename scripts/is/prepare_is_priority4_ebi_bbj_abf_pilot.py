#!/usr/bin/env python3
"""Prepare separately-versioned EBI reprocessed QTD000161 × BBJ IS exploratory ABF.
Fail-closed exact REF/ALT; never mutate inherited GTEx_v8 historical evidence.
An inferred genotype AN ceiling is *not* experimentally attested effective QTL N.
"""
import argparse,csv,json,hashlib,math
from collections import defaultdict,Counter
from pathlib import Path

IDS={"FGF5":"ENSG00000138675","C4orf22":"ENSG00000197826",
     "CALHM2":"ENSG00000138172","NEURL":"ENSG00000107954"}
LOCUS={"FGF5":"BBJ_IS_L001","C4orf22":"BBJ_IS_L001",
       "CALHM2":"BBJ_IS_L002","NEURL":"BBJ_IS_L002"}
FIELDS=["gene_symbol","gene_id","locus","variant_GRCh37_REF_ALT","variant_GRCh38_REF_ALT",
        "BBJ_beta_ALT","BBJ_se","BBJ_p","BBJ_EAF_ALT","BBJ_MAF","EBI_QTD_beta_ALT",
        "EBI_QTD_se","EBI_QTD_p","EBI_QTD_maf","EBI_QTD_AC","EBI_QTD_AN",
        "ebi_genotype_AN_over_2_not_QTL_effective_N","harmonization",
        "QTD_nominal_full_cis_tested","BBJ_GWAS_exact_FASTA_ref"]
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def numeric(v):
 x=float(v)
 if not math.isfinite(x):raise ValueError("Nonfinite source statistic")
 return x
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b""):h.update(chunk)
 return h.hexdigest()
def run(matched,qdir,out):
 if out.exists() and any(out.iterdir()):raise FileExistsError("Refuse overwriting pilot evidence")
 merged=read(matched)
 if len(merged)!=20689:raise ValueError("Full 4-gene cis EBI source changed")
 qtl={}
 for locus in ("BBJ_IS_L001","BBJ_IS_L002"):
  for r in read(qdir/f"{locus}.QTD000161.wide.priority4.tsv"):
   gene=r["gene_id"]
   if gene not in IDS.values():continue
   key=(gene,f"{r['chromosome']}:{r['position']}:{r['ref'].upper()}:{r['alt'].upper()}")
   if key in qtl:
    if any(qtl[key][field]!=r[field] for field in ("maf","an","ac","beta","se","pvalue")):
     raise ValueError("Mismatched QTL duplicate rsIDs")
   else:qtl[key]=r
 allout=[];counts=defaultdict(Counter);seen=set()
 for r in merged:
  sym=r["symbol"];gene=r["gene_id"]
  if sym not in IDS or IDS[sym]!=gene or r["locus"]!=LOCUS[sym]:
   raise ValueError("Unrecognized molecular/gene/locus identity")
  counts[sym]["original_full_cis"]+=1
  if r["liftover_status"]!="GRCH37_FASTA_REF_VERIFIED":
   counts[sym]["excluded_not_Fasta37_valid"]+=1;continue
  if r["BBJ_exact_Fasta37"]!="True":
   counts[sym]["excluded_not_BBJ_exact"]+=1;continue
  if r["palindromic_review"]=="True":
   counts[sym]["excluded_palindromic"]+=1;continue
  k=(gene,r["SNP_GRCh38"])
  if k not in qtl:raise ValueError("Missing QTD original full cis nominal record")
  q=qtl[k]
  alt=q["alt"].upper()
  if q["gene_id"]!=gene or q["variant"]!="chr"+q["chromosome"]+"_"+q["position"]+"_"+q["ref"]+"_"+q["alt"]:
   raise ValueError("QTL source variant ID is not matched to REF/ALT")
  if r["SNP_GRCh38"]!=f"{q['chromosome']}:{q['position']}:{q['ref']}:{alt}":
   raise ValueError("Full cis evidence and QTL genomic SNP differ")
  if r["SNP_GRCh37"].split(":")[2:] not in (
       [q["ref"],q["alt"]],
       [q["ref"].translate(str.maketrans("ACGT","TGCA")),q["alt"].translate(str.maketrans("ACGT","TGCA"))]):
   raise ValueError("Invalid GRCh38 to 37 strand mapping")
  beta=numeric(q["beta"]);se=numeric(q["se"]);p=numeric(q["pvalue"])
  bb=numeric(r["BBJ_beta_ALT"]);bp=numeric(r["BBJ_p"])
  an=numeric(q["an"]);ac=numeric(q["ac"]);maf=numeric(q["maf"])
  eaf=numeric(r["BBJ_AF_ALT"])
  # BBJ s.e. must come from the frozen canonical source. We avoid reconstructing
  # it using beta/p because p-values can be rounded to zero at very strong loci.
  # A companion BBJ GWAS exact SNP keyed SE table is joined separately below.
  if not (0<=p<=1 and 0<=bp<=1 and 0<se and 0<an and 0<=ac<=an and 0<maf<=.5 and 0<eaf<1):
   counts[sym]["invalid_numeric"]+=1;continue
  if abs(min(ac/an,1-ac/an)-maf)>2e-5:raise ValueError("QTL MAF does not match genotype AC/AN")
  key=(sym,r["SNP_GRCh37"])
  if key in seen:raise ValueError("Duplicate harmonized gene SNP")
  seen.add(key)
  allout.append(dict(gene_symbol=sym,gene_id=gene,locus=r["locus"],
   variant_GRCh37_REF_ALT=r["SNP_GRCh37"],variant_GRCh38_REF_ALT=r["SNP_GRCh38"],
   BBJ_beta_ALT=bb,BBJ_se="",BBJ_p=bp,BBJ_EAF_ALT=eaf,BBJ_MAF=min(eaf,1-eaf),
   EBI_QTD_beta_ALT=beta,EBI_QTD_se=se,EBI_QTD_p=p,EBI_QTD_maf=maf,
   EBI_QTD_AC=ac,EBI_QTD_AN=an,
   ebi_genotype_AN_over_2_not_QTL_effective_N=an/2,
   harmonization="ALT_AFTER_CHAIN_AND_FASTA_REF_EQUIVALENCE",
   QTD_nominal_full_cis_tested=True,BBJ_GWAS_exact_FASTA_ref=True))
  counts[sym]["preliminary_usable"]+=1
 out.mkdir(parents=True,exist_ok=True)
 path=out/"IS_PRIORITY4_EBI_BBJ_ABF_INPUT_PRESE.tsv"
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(allout)
 summary=dict(input_SNP_gene_rows=len(merged),candidate_not_independent_studies=4,
   preliminary_excluding_palindrome=len(allout),
   QTL_effective_sample_size_proven=False,
   BBJ_effect_se_must_be_exactly_joined=True,
   scientific_status="PREP_ONLY_REQUIRE_CANONICAL_GWAS_SE_AND_QTL_N",
   per_gene={k:dict(v) for k,v in counts.items()},
   source_hash={"GWAS_annotated":sha(matched),
   **{f.name:sha(f) for f in sorted(qdir.glob("BBJ_IS_L00*.QTD000161.wide.priority4.tsv"))}})
 (out/"IS_PRIORITY4_EBI_BBJ_ABF_PILOT_PRESE.json").write_text(json.dumps(summary,indent=2)+"\n")
 return summary

if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("matched","qdir","out"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(run(a.matched,a.qdir,a.out),indent=2))
