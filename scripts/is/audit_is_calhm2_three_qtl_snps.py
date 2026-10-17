#!/usr/bin/env python3
"""CALHM2 lead-linked QTL SNP evidence ledger (read-only, not colocalization)."""
import argparse,csv,json
from pathlib import Path

LEADS=("BBJ","GIGASTROKE")
POPS=("EAS504","JPT104","otherEAS400")
FIELDS=("SNP_GRCh37","SNP_GRCh38","QTD000161_nominal_p","QTD000161_beta",
        "QTD000161_se","BBJ_p","BBJ_beta_ALT","GIGA_EAS_AIS_p",
        "GIGA_EAS_AIS_beta_ALT")
def rows(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def create(evidence,ld,out,cutoff=1e-4):
 e=rows(evidence)
 q={r["SNP_GRCh37"]:r for r in e if r["locus"]=="BBJ_IS_L002"
    and r["symbol"]=="CALHM2" and r["BBJ_exact_Fasta37"]=="True"
    and r["GIGA_EAS_AIS_exact_Fasta37"]=="True"
    and r["palindromic_review"]=="False"
    and float(r["QTD000161_nominal_p"])<cutoff}
 if not q:raise ValueError("No harmonized significant CALHM2 QTL variants")
 if len(q)!=3:raise ValueError(f"Unexpected CALHM2 QTL count {len(q)}; examine universe before analysis")
 L={(r["population"],r["lead"],r["variant_b37"]):r for r in rows(ld)}
 ledger=[]
 for key in sorted(q,key=lambda k:int(k.split(":")[1])):
  r=q[key]
  rec={"gene":"CALHM2","SNP_GRCh37":key}
  for col in FIELDS:
   rec[col]=r[col]
  for pop in POPS:
   for lead in LEADS:
    tag=f"{pop}_{lead}"
    item=L.get((pop,lead,key))
    rec[tag+"_signed_r"]=item["signed_r"] if item else ""
    rec[tag+"_r2"]=item["r2"] if item else ""
    rec[tag+"_tag_r2_ge_0p8"]=item["r2_ge_0p8"] if item else ""
  rec["statistical_evidence"]="NOMINAL_QTL_AND_PUBLIC_REFERENCE_LD_ONLY"
  ledger.append(rec)
 for pop in ("EAS504","JPT104"):
  for lead in LEADS:
   n=sum(x[f"{pop}_{lead}_tag_r2_ge_0p8"]=="True" for x in ledger)
   if n!=(3 if lead=="BBJ" else 0):raise ValueError("Prior lead tag count not reproduced")
 out.mkdir(parents=True,exist_ok=False)
 with (out/"CALHM2_THREE_QTL_SNP_LEDGER.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(ledger[0]),delimiter="\t");w.writeheader();w.writerows(ledger)
 meta={"status":"EXPLORATORY_VARIANT_LEDGER_NOT_CREDIBLE_SET",
       "n_CALHM2_qtl_variants":len(ledger),"qtl_p_cutoff":cutoff,
       "no_conditioned_gwas":True,"no_multisignal_coloc":True,
       "source_evidence":str(evidence),"source_reference_ld":str(ld),
       "discovery_universe_2225_unchanged":True}
 (out/"CALHM2_THREE_QTL_SNP_MANIFEST.json").write_text(json.dumps(meta,indent=2)+"\n")
 return ledger
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--evidence",type=Path,required=True)
 p.add_argument("--ld",type=Path,required=True)
 p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(create(a.evidence,a.ld,a.out),indent=2))
