#!/usr/bin/env python3
"""Read-only 30 ABF assay multi-signal ancestry readiness, no false coloc."""
import csv,json,argparse,hashlib
from pathlib import Path
from collections import Counter

def records(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def run(replay,genes,ldfile,out):
 rr=[x for x in records(replay) if x["status"]=="PASS"]
 gg={x["gene_id_stable"]:x["canonical_GENCODE_v19_gene_symbol"] for x in records(genes)}
 ld=records(ldfile)
 if len(rr)!=30 or len(ld)!=8:raise ValueError("Prior source cardinality changed")
 seen=set();results=[]
 for x in rr:
  key=(x["locus"],x["dataset_key"],x["gene_base"])
  if key in seen:raise ValueError("Duplicate assay")
  seen.add(key)
  symbol=gg.get(x["gene_base"])
  aliases={symbol}
  if symbol=="NEURL":aliases.add("NEURL1")
  attached=[v for v in ld if v["locus"]==x["locus"] and v["assay"]==x["dataset_key"] and v["gene"] in aliases]
  if len(attached)>1:raise ValueError("Ambiguous LD assay join")
  match=attached[0] if attached else None
  if match and int(match["n_input"])!=int(x["actual_nsnps"]):raise ValueError("LD and coloc SNP denominator mismatch")
  results.append(dict(gene_id=x["gene_base"],symbol=symbol,locus=x["locus"],dataset_key=x["dataset_key"],
   original_assay_SNPs=x["actual_nsnps"],archived_PP_H4=x["PP_H4"],
   available_EAS_reference_LD_alignment=bool(match),
   EAS_ref_SNP_coverage_fraction=(int(match["n_EAS_LD_match_GRCh37"])/int(match["n_input"]) if match else None),
   EUR_GTEx_in_study_QTL_signed_LD_obtained=False,
   independent_BBJ_GWAS_in_study_signed_LD_obtained=False,
   multi_signal_SuSiE_coloc_performed=False,
   study_matched_full_QTL_cis_set_verified=False,
   readiness="REFERENCE_GWAS_LD_ONLY_QTL_LD_MISSING" if match else "NO_ASSAY_SPECIFIC_SIGNED_LD_ATTESTATION",
   inferred_causal_gene=False))
 out.mkdir(parents=True,exist_ok=True)
 results.sort(key=lambda x:-float(x["archived_PP_H4"]))
 with (out/"IS_30_MULTISIGNAL_ANCESTRY_GATES.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter="\t");w.writeheader();w.writerows(results)
 j=dict(assays=30,genes=len({r["gene_id"] for r in results}),
    assays_EAS_ref_LD_SNP_covered=sum(x["available_EAS_reference_LD_alignment"] for x in results),
    assays_QTL_cohort_LD_available=0,
    assays_ready_for_validated_multisignal_coloc=0,
    source_sha256={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in [("replay",replay),("genes",genes),("reference_ld",ldfile)]},
    warnings=["EAS reference LD matches BBJ ancestry approximately but not GTEx QTL cohort","EAS and GTEx EUR LD must not be silently substituted","Missing LD is not a negative molecular genetic result"])
 (out/"IS_30_MULTISIGNAL_ANCESTRY_GATES.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 a=argparse.ArgumentParser()
 for n in ["replay","genes","ldfile","out"]:a.add_argument("--"+n,type=Path,required=True)
 p=a.parse_args()
 print(json.dumps(run(p.replay,p.genes,p.ldfile,p.out),indent=2))
