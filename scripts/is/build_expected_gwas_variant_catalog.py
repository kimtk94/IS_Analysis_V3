#!/usr/bin/env python3
"""Build expected per-locus GWAS SNP catalog from source-verifiable cached 30 inputs.

Comparisons across 646 will be performed in Drive one-cell Colab; G1 native
BBJ/chain proof for the cached union remains the authoritative baseline.
"""
import csv,hashlib,collections
from pathlib import Path
source=Path("/srv/is-analysis/results/is/audits/legacy_coloc_snp_replay_20261010_v1/input")
ref=Path("/srv/is-analysis/results/is/audits/is_g1_native_bbj_chain_20261011_v1/IS_G1_NATIVE_BBJ_ALLELE_CHAIN_SNP.tsv")
a={}
for p in sorted(source.glob("BBJ_IS_L00*__GTEx_V8__*.tsv")):
 with p.open() as f:
  for z in csv.DictReader(f,delimiter="\t"):
   k=(z["locus"],z["variant_id"])
   value=(z["match_key"],z["gwas_beta"],z["gwas_eaf"])
   if k in a and a[k]!=value:raise ValueError("SOURCE_SNP_ASSOCIATION_INCONSISTENT")
   a[k]=value
with ref.open() as f:
 r={(z["locus"],z["variant_id"]):z for z in csv.DictReader(f,delimiter="\t")}
assert len(a)==8309 and len(r)==8309 and set(a)==set(r)
for k,(mk,beta,eaf) in a.items():
 z=r[k]
 assert mk==z["match_key"] and beta==z["source_gwas_beta"] and eaf==z["source_gwas_eaf"]
out=Path("/srv/is-analysis/results/is/audits/is_g1_native_bbj_chain_20261011_v1/IS_BBJ_FOUR_LOCUS_EXPECTED_GWAS_CATALOG.tsv")
assert not out.exists()
with out.open("w") as h:
 w=csv.writer(h,delimiter="\t",lineterminator="\n")
 w.writerow(["locus","variant_id","match_key","gwas_beta","gwas_eaf"])
 for (locus,variant),v in sorted(a.items()):
  w.writerow([locus,variant,*v])
print("EXPECTED_FULL_646_GWAS_CATALOG_UNIQUE_VARIANTS",len(a))
print("EXPECTED_FULL_646_GWAS_CATALOG_SHA256",hashlib.sha256(out.read_bytes()).hexdigest())
print("PER_LOCUS",dict(collections.Counter(k[0] for k in a)))
