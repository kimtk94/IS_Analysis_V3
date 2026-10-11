#!/usr/bin/env python3
"""Four priority IS genes: exact public EBI full-cis QTL × BBJ/GIGASTROKE GWAS.

Reuses the frozen, independently FASTA-anchored GWAS allele/effect contract.
No mixed-source coloc, no independent replication, no causal-gene promotion.
"""
import argparse,csv,hashlib,json,sys
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED={"FGF5":5511,"C4orf22":5490,"CALHM2":4868,"NEURL":4820}
INPUT_GENES={"FGF5":"ENSG00000138675","C4orf22":"ENSG00000197826",
             "CALHM2":"ENSG00000138172","NEURL":"ENSG00000107954"}
THRESH_QTL=1e-4
THRESH_GWAS=5e-8

def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save(path,records):
 if not records:raise ValueError("Refuse empty science output")
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
  w.writeheader();w.writerows(records)
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()

def prepare(source,lo,fa,canon):
 inp=read(source)
 if len(inp)!=sum(EXPECTED.values()):raise ValueError("EBI source should be 20689 records")
 sets=defaultdict(set);mapping={};qc=defaultdict(Counter);seen={}
 for x in inp:
  symbol=x["symbol"]
  if symbol not in EXPECTED or x["gene_id"]!=INPUT_GENES[symbol]:
   raise ValueError("Source gene identity drift")
  if x["status"]!="EBI_REPROCESSED_CIS_TESTED_VARIANT_ONLY":
   raise ValueError("Not an independently reviewed EBI reprocessed cis source")
  try:
   chrom=x["chromosome"];pos=int(x["position_GRCh38"]);ref=x["REF_GRCh38"];alt=x["ALT_GRCh38"]
   p=float(x["QTD000161_nominal_QTL_p"])
  except ValueError as e:raise ValueError("Incompatible source schema") from e
  if not 0<=p<=1:raise ValueError("Invalid QTL nominal p")
  key=(symbol,chrom,pos,ref,alt)
  if key in seen:raise ValueError("Duplicate gene-SNP full cis record")
  seen[key]=True
  sets[symbol].add(key)
  stat,build37,strand=canon.canonical_qtl(chrom,pos,ref,alt,lo,fa,{})
  qc[symbol][stat]+=1
  record=dict(locus=x["locus"],symbol=symbol,gene_id=x["gene_id"],
    SNP_GRCh38=x["variant_id_GRCh38"],SNP_GRCh37=build37,
    liftover_status=stat,liftover_strand=strand or "",
    in_archived_original_QTL_BBJ_intersection=x["in_original_archived_GWAS_QTL_intersection"],
    QTD000161_nominal_p=p,QTD000161_beta=x["QTD000161_beta"],QTD000161_se=x["QTD000161_se"],
    QTD000161_source="EBI_REPROCESSED_NOT_HISTORIC_ORIGINAL_FASTQTL",
    BBJ_exact_Fasta37=False,BBJ_p="",BBJ_beta_ALT="",BBJ_AF_ALT="",
    GIGA_EAS_AIS_exact_Fasta37=False,GIGA_EAS_AIS_p="",GIGA_EAS_AIS_beta_ALT="",GIGA_EAS_AIS_AF_ALT="",
    GIGA_EAS_AIS_allele_order_provenance="",palindromic_review=False,
    same_variant_QTL_BBJ_joint_marginal_threshold=False,
    same_variant_QTL_GIGA_joint_marginal_threshold=False,
    scientific_gate="INPUT_READY_DIAGNOSTIC_NOT_COLOCALIZATION")
  if stat=="GRCH37_FASTA_REF_VERIFIED":
   if build37 in mapping and mapping[build37]["SNP_GRCh38"]!=record["SNP_GRCh38"]:
    # Fail closed rather than silently choosing a many-to-one chain mapping.
    raise ValueError("Distinct original QTL SNPs collapse onto identical GRCh37 variant")
   mapping.setdefault(build37,{"SNP_GRCh38":record["SNP_GRCh38"],"records":[]})["records"].append(record)
 return inp,mapping,qc,sets

def run(args):
 root=Path(__file__).resolve().parent
 sys.path.insert(0,str(root))
 sys.path.insert(0,str(args.vendor))
 from pyliftover import LiftOver
 import audit_is_eas5_public_qtl_gwas_snp_overlap as canon
 import audit_is_eas5_98_egene_fullcis_gwas_exact as frozen
 canon.validate_reverse_chain(args.chain)
 if args.out.exists() and any(args.out.iterdir()):
  raise FileExistsError("Do not overwrite science output "+str(args.out))
 fasta=canon.Fasta37(args.fasta)
 try:
  lo=LiftOver(str(args.chain))
  input_records,mapping,qc,sets=prepare(args.full_cis,lo,fasta,canon)
  positions=defaultdict(set)
  for v in mapping:
   ch,pos,*_=v.split(":")
   positions[ch].add(int(pos))
  bbj,bbjstats=frozen.scan_gwas(args.bbj,"BBJ_IS",positions,canon,fasta)
  giga,gigastats=frozen.scan_gwas(args.giga,"EAS_AIS",positions,canon,fasta,allow_fasta_recoded_giga=True)
 finally:fasta.close()
 for variant,group in mapping.items():
  b=bbj.get(variant)
  a=giga.get(variant)
  chrom,pos,ref,alt=variant.split(":")
  pal=frozenset((ref,alt)) in (frozenset(("A","T")),frozenset(("C","G")))
  for x in group["records"]:
   x["palindromic_review"]=pal
   if b:
    x["BBJ_exact_Fasta37"]=True
    x["BBJ_p"]=b["p"];x["BBJ_beta_ALT"]=b["beta_alt"];x["BBJ_AF_ALT"]=b["af_alt"]
    x["same_variant_QTL_BBJ_joint_marginal_threshold"]=x["QTD000161_nominal_p"]<THRESH_QTL and float(b["p"])<THRESH_GWAS and not pal
   if a:
    x["GIGA_EAS_AIS_exact_Fasta37"]=True
    x["GIGA_EAS_AIS_p"]=a["p"];x["GIGA_EAS_AIS_beta_ALT"]=a["beta_alt"];x["GIGA_EAS_AIS_AF_ALT"]=a["af_alt"]
    x["GIGA_EAS_AIS_allele_order_provenance"]=a["historical_giga_id_order"]
    x["same_variant_QTL_GIGA_joint_marginal_threshold"]=x["QTD000161_nominal_p"]<THRESH_QTL and float(a["p"])<THRESH_GWAS and not pal
 records=[]
 for group in mapping.values():records.extend(group["records"])
 # Include all source records that failed SNP liftover rather than hiding excluded rows.
 successful={(r["symbol"],r["SNP_GRCh38"]) for r in records}
 for x in input_records:
  k=(x["symbol"],x["variant_id_GRCh38"])
  if k in successful:continue
  records.append(dict(locus=x["locus"],symbol=x["symbol"],gene_id=x["gene_id"],
    SNP_GRCh38=x["variant_id_GRCh38"],SNP_GRCh37="",liftover_status="NO_VALID_BIALLELIC_FASTA_SNP",
    liftover_strand="",in_archived_original_QTL_BBJ_intersection=x["in_original_archived_GWAS_QTL_intersection"],
    QTD000161_nominal_p=float(x["QTD000161_nominal_QTL_p"]),
    QTD000161_beta=x["QTD000161_beta"],QTD000161_se=x["QTD000161_se"],
    QTD000161_source="EBI_REPROCESSED_NOT_HISTORIC_ORIGINAL_FASTQTL",
    BBJ_exact_Fasta37=False,BBJ_p="",BBJ_beta_ALT="",BBJ_AF_ALT="",
    GIGA_EAS_AIS_exact_Fasta37=False,GIGA_EAS_AIS_p="",GIGA_EAS_AIS_beta_ALT="",GIGA_EAS_AIS_AF_ALT="",
    GIGA_EAS_AIS_allele_order_provenance="",palindromic_review=False,
    same_variant_QTL_BBJ_joint_marginal_threshold=False,
    same_variant_QTL_GIGA_joint_marginal_threshold=False,
    scientific_gate="INPUT_READY_DIAGNOSTIC_NOT_COLOCALIZATION"))
 if len(records)!=sum(EXPECTED.values()):raise ValueError("Input gene-SNP universe not fully accounted")
 records.sort(key=lambda r:(r["locus"],r["symbol"],r["SNP_GRCh38"]))
 summaries=[]
 for sym in EXPECTED:
  rr=[r for r in records if r["symbol"]==sym]
  if len(rr)!=EXPECTED[sym]:raise ValueError("Full-cis per-gene denominator drift")
  ok=[r for r in rr if r["liftover_status"]=="GRCH37_FASTA_REF_VERIFIED"]
  common=[r for r in rr if str(r["in_archived_original_QTL_BBJ_intersection"])=="True"]
  fresh=[r for r in rr if str(r["in_archived_original_QTL_BBJ_intersection"])=="False"]
  summaries.append(dict(symbol=sym,full_ebi_reprocessed_cis_snps=len(rr),
    biallelic_snp_FASTA37_verified=len(ok),
    SNP_not_valid_mappable=len(rr)-len(ok),
    exact_BBJ_GWAS_overlap=sum(bool(r["BBJ_exact_Fasta37"]) for r in ok),
    exact_GIGASTROKE_EAS_AIS_overlap=sum(bool(r["GIGA_EAS_AIS_exact_Fasta37"]) for r in ok),
    exact_BOTH_BBJ_GIGA_overlap=sum(bool(r["BBJ_exact_Fasta37"] and r["GIGA_EAS_AIS_exact_Fasta37"]) for r in ok),
    exact_BOTH_nonpal_overlap=sum(bool(r["BBJ_exact_Fasta37"] and r["GIGA_EAS_AIS_exact_Fasta37"] and not r["palindromic_review"]) for r in ok),
    archived_original_shared_count=len(common),
    EBI_only_vs_original_QTL_BBJ_intersection_count=len(fresh),
    EBI_only_with_BBJ_GWAS=sum(bool(r["BBJ_exact_Fasta37"]) for r in fresh),
    EBI_only_with_GIGA_EAS_AIS=sum(bool(r["GIGA_EAS_AIS_exact_Fasta37"]) for r in fresh),
    EBI_only_with_BOTH_GWAS=sum(bool(r["BBJ_exact_Fasta37"] and r["GIGA_EAS_AIS_exact_Fasta37"]) for r in fresh),
    same_QTL_BBJ_nominal_joint_threshold_not_coloc=sum(bool(r["same_variant_QTL_BBJ_joint_marginal_threshold"]) for r in ok),
    same_QTL_GIGA_AIS_nominal_joint_threshold_not_coloc=sum(bool(r["same_variant_QTL_GIGA_joint_marginal_threshold"]) for r in ok),
    genotype_level_GTEx_QTL_LD_available=False,
    shared_causal_variant_proven=False))
 args.out.mkdir(parents=True,exist_ok=True)
 save(args.out/"IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_GENE_READINESS.tsv",summaries)
 save(args.out/"IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv",records)
 j=dict(status="COMPLETE_REPROCESSED_QTD000161_CIS_BY_BBJ_GIGA_EXACT_FASTA_GRCh37",
   source_is_reprocessed_ebi_gtEx_NOT_original_coloc_QTL=True,
   gene_source_records=len(records),
   distinct_QTL_GRCh37_SNP=len(mapping),
   studywide_original_bbj_snp_records=bbjstats["source_rows"],
   studywide_eas_ais_snp_records=gigastats["source_rows"],
   source_GWAS_QC={"BBJ":bbjstats,"GIGA_EAS_AIS":gigastats},
   full_cis_ebi_study_universe_tested_complete=all(r["full_ebi_reprocessed_cis_snps"]==EXPECTED[r["symbol"]] for r in summaries),
   original_GTEx_V8_full_cis_tested_variants_authenticated=False,
   GTEx_in_study_molecular_LD_verified=False,
   original_2225_positional_candidate_genes_unchanged=True,
   new_coloc_performed=False,
   BBJ_GIGA_independent_replication_claimed=False,
   causal_gene_count=0,
   summaries=summaries,
   source_sha256={p.name:sha(p) for p in (args.full_cis,args.chain)})
 (args.out/"IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGA_MANIFEST.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 ap=argparse.ArgumentParser()
 for k in ("full-cis","bbj","giga","chain","fasta","vendor","out"):
  ap.add_argument("--"+k,type=Path,required=True)
 a=ap.parse_args()
 r=run(a)
 print("OUTPUT_SUMMARY",json.dumps({k:v for k,v in r.items() if k not in ("summaries","source_GWAS_QC","source_sha256")},indent=2),flush=True)
 for x in r["summaries"]:print("GENE",x,flush=True)
