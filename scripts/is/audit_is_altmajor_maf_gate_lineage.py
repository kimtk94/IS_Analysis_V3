#!/usr/bin/env python3
"""Trace original GTEx ALT-major variant omissions through BBJ+1KG EAS QC.

This is a read-only forensic audit of *existing* native GWAS, genotype QC log,
1KG EAS reference .pvar, post-QC .pvar, stage3 SuSiE input and the prior
10-site diagnostic ledger. It does NOT retroactively change coloc SNP sets.
"""
import argparse,csv,json,re,hashlib
from collections import Counter
from pathlib import Path

def read_tsv(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
 return h.hexdigest()
def locate(pvar,variant):
 # Both source and post-QC pvar have columns #CHROM POS ID REF ALT ...
 records=[]
 with pvar.open() as f:
  rd=csv.DictReader((line for line in f if not line.startswith("##")),delimiter="\t")
  if rd.fieldnames is None or not {"#CHROM","POS","ID","REF","ALT"}.issubset(rd.fieldnames):
   raise ValueError("Invalid 1KG EAS PVAR schema "+str(rd.fieldnames))
  for r in rd:
   if r["ID"]==variant:
    records.append(r)
 return records
def audit(ten,ref,filtered,log_text,susie):
 if len(ten)!=10 or len({x["site"] for x in ten})!=10:
  raise ValueError("Exactly ten unique raw GTEx ALT-major source positions required")
 retained=[x for x in ten if x["native_status"]=="BBJ_ORIGINAL_SAME_REF_ALT_PAIR"]
 missing=[x for x in ten if x["native_status"]=="NOT_REPORTED_AT_GRCH37_POSITION"]
 if len(retained)!=1 or len(missing)!=9:
  raise ValueError("Native 10-site original GWAS coverage evidence changed")
 row=retained[0]
 vid=row["native_variant_candidates"]
 if vid!="12:112090022:C:A" or row["site"]!="12:111652218:C:A":
  raise ValueError("The native BBJ ALT-major followup marker identity changed")
 src=locate(ref,vid)
 out=locate(filtered,vid)
 if len(src)!=1 or len(out)!=0:
  raise ValueError("1KG source variant must be present once and absent after QC")
 r=src[0]
 if (r["#CHROM"],r["POS"],r["REF"],r["ALT"])!=("12","112090022","C","A"):
  raise ValueError("PVAR allele content drift")
 raw_info=r.get("INFO","")
 af_match=re.search(r"(?:^|;)EAS_AF=([0-9.]+)(?:;|$)",raw_info)
 if not af_match:raise ValueError("1KG source EAS AF missing")
 af=float(af_match.group(1));maf=min(af,1-af)
 if abs(af-0.999)>0.0001 or maf>=0.01:
  raise ValueError("The EAS SNP is not MAF < 0.01; cannot attribute exclusion")
 if not re.search(r"(?m)^\s*--maf 0\.01\s*$",log_text) or not re.search(r"(?m)^\s*--geno 0\.05\s*$",log_text):
  raise ValueError("Expected historical PLINK QC command unavailable")
 if "0 variants removed due to missing genotype data" not in log_text:
  raise ValueError("Cannot attribute site exclusion exclusively to the MAF gate")
 counts=re.search(r"(\d+) variants removed due to allele frequency threshold",log_text)
 if not counts:raise ValueError("MAF filtered count not in original QC log")
 if not any(x["variant_id"]==vid for x in susie):
  stage3_status="NOT_IN_STAGE3_SUSIE_V3"
 else:
  raise ValueError("Variant unexpectedly in Stage3 despite QC exclusion")
 minimum=min(float(x["maf"]) for x in susie)
 if minimum<0.01:
  raise ValueError("Stage3 finemap input contains MAF<0.01 unexpectedly")
 return {
    "status":"ONE_NATIVE_BBJ_SITE_1KG_EAS_MAF_FILTER_LINEAGE_CONFIRMED",
    "n_original_alt_major_sites":len(ten),
    "n_native_BBJ_GWAS_absent":len(missing),
    "n_native_BBJ_GWAS_present":len(retained),
    "native_variant_GRCh37":vid,
    "QTL_variant_GRCh38":row["site"],
    "1kg_EAS_ref_present":True,
    "1kg_EAS_alt_AF_recorded":af,
    "1kg_EAS_MAF_derived":maf,
    "PLINK_filter_maf_threshold":0.01,
    "PLINK_filter_geno_threshold":0.05,
    "PLINK_variants_removed_by_maf":int(counts.group(1)),
    "1kg_EAS_post_QC_variant_absent":True,
    "stage3_BBJ_susie_input_variant_status":stage3_status,
    "stage3_BBJ_susie_input_min_maf":minimum,
    "BBJ_native_effect_ALT_af":float(row["native_allele2_AF"]),
    "BBJ_native_p":float(row["native_p"]),
    "saved_GTEx_ALT_AF_range":[float(row["raw_QTD_ALT_AF_min"]),float(row["raw_QTD_ALT_AF_max"])],
    "primary_absence_explanation":"REFERENCE_1KG_EAS_MAF_FILTER_AND_UPSTREAM_STAGE3_SNP_UNIVERSE",
    "full_coloc_variant_universe_not_modified":True,
    "LD_ancestry_matched_to_BBJ_study":False,
    "posthoc_SNP_reinclusion":False,
    "causal_gene_established":False,
    "caveat":"PLINK QC with --maf 0.01 demonstrably removed this EAS rare variant. Original Stage3 SNP provenance is limited to the available files; no retroactive re-inclusion or causal signal inference."
 }
def main():
 p=argparse.ArgumentParser()
 for name in ("ten-sites","reference-pvar","post-qc-pvar","plink-qc-log","stage3-bbj-input","out-dir"):
  p.add_argument("--"+name,type=Path,required=True)
 a=p.parse_args()
 paths={"ten":a.ten_sites,"reference":a.reference_pvar,"post_QC":a.post_qc_pvar,
        "PLINK_QC":a.plink_qc_log,"stage3":a.stage3_bbj_input}
 for f in paths.values():
  if not f.is_file():raise FileNotFoundError(f)
 result=audit(read_tsv(a.ten_sites),a.reference_pvar,a.post_qc_pvar,
              a.plink_qc_log.read_text(),read_tsv(a.stage3_bbj_input))
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("Refuse overwrite")
 a.out_dir.mkdir(parents=True,exist_ok=True)
 result["sources_sha256"]={k:sha(f) for k,f in paths.items()}
 result["source_files"]={k:str(f) for k,f in paths.items()}
 (a.out_dir/"IS_GTEX_ALTMAJOR_ORIGINAL_BBJ_LD_MAF_GATE.json").write_text(json.dumps(result,indent=2)+"\n")
 print("IS_GTEX_BBJ_LD_MAF_GATE_CONFIRMED",json.dumps({k:result[k] for k in (
  "n_native_BBJ_GWAS_absent","n_native_BBJ_GWAS_present","native_variant_GRCh37",
  "1kg_EAS_MAF_derived","PLINK_filter_maf_threshold","stage3_BBJ_susie_input_variant_status")}))

if __name__=="__main__":main()
