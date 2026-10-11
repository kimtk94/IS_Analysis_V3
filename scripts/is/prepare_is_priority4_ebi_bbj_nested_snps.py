#!/usr/bin/env python3
"""Controlled SNP-universe sensitivity: fixed EBI-QTD × BBJ stats, nested old-overlap SNPs."""
import argparse,csv,json,hashlib
from pathlib import Path
EXPECTED={"FGF5":1286,"C4orf22":1286,"CALHM2":1342,"NEURL":1342}
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def run(full,marker,out):
 if out.exists() and any(out.iterdir()):raise FileExistsError("Refuse overwrite")
 rows=read(full);r=read(marker)
 if len(rows)!=11050 or len(r)!=20689:raise ValueError("Input SNP denominator changed")
 old={}
 for x in r:
  if x["liftover_status"]!="GRCH37_FASTA_REF_VERIFIED":continue
  key=(x["symbol"],x["SNP_GRCh37"])
  if key in old:raise ValueError("Multiple original variants map to one genomic key")
  old[key]=x["in_archived_original_QTL_BBJ_intersection"]=="True"
 nested=[]
 for row in rows:
  key=(row["gene_symbol"],row["variant_GRCh37_REF_ALT"])
  if key not in old:raise ValueError("Full source variant missing provenance")
  if old[key]:nested.append(row)
 by={gene:sum(x["gene_symbol"]==gene for x in nested) for gene in EXPECTED}
 if by!=EXPECTED:raise ValueError("Expected old-source nested allele-aligned counts drift "+str(by))
 out.mkdir(parents=True,exist_ok=True)
 dest=out/"IS_PRIORITY4_SAME_EBI_BBJ_ARCHIVED_SNP_NESTED_INPUT.tsv"
 with dest.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(nested)
 j=dict(full_compatible_geneSNP=11050,nested_compatible_geneSNP=len(nested),
  per_gene_nested=by,source_statistics_identical_between_nested_and_full=True,
  only_SNP_universe_changed=True, original_historical_GTEx_v8_FastQTL_analysis_reproduced=False,
  interpretation="Controlled same-source nested SNP coverage sensitivity, not causal replication")
 (out/"IS_PRIORITY4_SAME_SOURCE_NESTED_SNP_MANIFEST.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("full","marker","out"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args();print(json.dumps(run(a.full,a.marker,a.out),indent=2))
