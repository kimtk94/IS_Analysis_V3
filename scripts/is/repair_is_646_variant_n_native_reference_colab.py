#!/usr/bin/env python3
"""Repair V1's unjustified 8309 == full 646-source variant-union assertion.

Uses a source-native original BBJ ZIP/chain derived fingerprint reference for
all source variants, and the 30 frozen source file checksums for provenance.
No raw GWAS summary statistics shipped, just cryptographic fingerprints.
"""
import argparse,ast,hashlib,json,subprocess
from pathlib import Path

def build(source,ref_gz,cache30_json,out):
    if out.exists():raise FileExistsError("DO_NOT_OVERWRITE_EXISTING_NOTEBOOK")
    original=json.loads(source.read_text())
    if len(original["cells"])!=1:raise ValueError("ONECELL_SOURCE_REQUIRED")
    py="".join(original["cells"][0]["source"])
    start="# Exact independent native BBJ chain source matching across ALL 646 QTL inputs."
    stop='if shutil.which("Rscript") is None:'
    if py.count(start)!=1 or py.count(stop)!=1:raise ValueError("PREFLIGHT_ANCHORS_MISSING")
    expected_sha=hashlib.sha256(ref_gz.read_bytes()).hexdigest()
    with open(cache30_json) as f: cache30=json.load(f)
    if len(cache30)!=30 or len(set(cache30))!=30:raise ValueError("FROZEN_CACHE30_INVALID")
    url="https://raw.githubusercontent.com/kimtk94/IS_Analysis_V3/feat/is-p0-evidence-20261010/data/is_reference/"+ref_gz.name
    preflight = '''# FIX: 8309 unique variants covered only the earlier 30-source local cache.
# The original 646 source assays contain 8956 distinct GWAS variants.
# Audit extra variants against native Sakaue/Kanai GWAS Allele2/BETA/EAF
# + independent hg19ToHg38 chain by using server-generated SHA fingerprints.
# The reference contains NO original native GWAS summary-stat beta/EAF values.
import gzip
import io
import urllib.request
import math
NATIVE_REFERENCE_URL = '''+repr(url)+'''
NATIVE_REFERENCE_SHA256 = '''+repr(expected_sha)+'''
CACHED_30_SHA256 = '''+repr(cache30)+'''
with urllib.request.urlopen(NATIVE_REFERENCE_URL,timeout=120) as response:
    native_reference_gzip=response.read()
assert hashlib.sha256(native_reference_gzip).hexdigest()==NATIVE_REFERENCE_SHA256, "NATIVE_REFERENCE_TAMPER_OR_URL_MISMATCH"
native_reference={}
with gzip.open(io.BytesIO(native_reference_gzip),"rt",encoding="utf8",newline="") as nh:
    for r in csv.DictReader(nh,delimiter="\\t"):
        key=(r["locus"],r["variant_id"])
        if key in native_reference:
            raise ValueError("DUPLICATE_NATIVE_GWAS_REFERENCE_VARIANT")
        native_reference[key]=r["native_variant_fingerprint_sha256_160bits"]
assert len(native_reference)==56189, "NATIVE_REFERENCE_4LOCUS_INPUT_COUNT_CHANGED"
print("ORIGINAL_NATIVE_BBJ_CHAIN_FINGERPRINT_REFERENCE_COUNT=",len(native_reference))
def native_fingerprint(locus,variant,match_key,beta,eaf):
    canonical="IS_BBJ_NATIVE_4LOCI_V1\\t"+"\\t".join((
         locus,variant,match_key,format(float(beta),".12f"),format(float(eaf),".12f")))
    return hashlib.sha256(canonical.encode()).hexdigest()[:40]
for filename,expected in CACHED_30_SHA256.items():
    fn=INPUT/filename
    if not fn.is_file() or sha(fn)!=expected:
        raise RuntimeError("ORIGINAL_CACHED30_FROZEN_SOURCE_SHA256_MISMATCH: "+filename)
print("ORIGINAL_CACHED_30_FILE_SHA256_VERIFIED=",len(CACHED_30_SHA256))
SNP_CATALOG=WORK/"IS_FULL_646_GWAS_VARIANT_NATIVE_STATUS.tsv"
combined={}
old_keys=set()
for fn in input_files:
    with fn.open(encoding="utf8",newline="") as infile:
        for rec in csv.DictReader(infile,delimiter="\\t"):
            key=(rec["locus"],rec["variant_id"])
            sample=(rec["match_key"],rec["gwas_beta"],rec["gwas_eaf"])
            if key in combined and combined[key]!=sample:
                raise RuntimeError("GWAS_SOURCE_CONFLICT_ACROSS_GENE_TISSUE_ASSAYS: "+str(key))
            combined[key]=sample
            if fn.name in CACHED_30_SHA256:
                old_keys.add(key)
print("ALL_ORIGINAL_646_UNIQUE_BBJ_GWAS_VARIANTS=",len(combined))
print("PREVIOUS_NATIVE_VERIFIED_CACHED30_SNP_UNION=",len(old_keys))
print("NEW_VARIANTS_BETWEEN_646_AND_CACHED30=",len(set(combined)-old_keys))
assert len(old_keys)==8309, "ORIGINAL_NATIVE_VERIFIED_CACHE30_SET_CHANGED"
assert len(combined)==8956, "FULL_ORIGINAL_646_UNIQUE_VARIANT_UNIVERSE_CHANGED"
assert old_keys.issubset(combined), "OLD_VARIANTS_MISSING_FROM_ALL_646"
assert len(set(combined)-old_keys)==647, "EXPECTED_647_ADDITIONAL_VARIANTS"
native_status=[]
for (locus,variant),(mk,beta,eaf) in sorted(combined.items()):
    source_sha=native_fingerprint(locus,variant,mk,beta,eaf)
    known=native_reference.get((locus,variant))
    status=("PASS_NATIVE_BBJ_ALLELE2_BETA_EAF_GRCH38_CHAIN"
            if known==source_sha else
            ("NOT_FOUND_IN_EXTENDED_NATIVE_REFERENCE" if known is None else "NATIVE_EFFECT_OR_CHAIN_MISMATCH"))
    native_status.append((locus,variant,mk,
       "CACHED30" if (locus,variant) in old_keys else "NEW647",status,source_sha))
with SNP_CATALOG.open("w",encoding="utf8",newline="") as handle:
    w=csv.writer(handle,delimiter="\\t",lineterminator="\\n")
    w.writerow(["locus","variant_id","match_key","original_subset",
                "native_BBJ_effect_chain_validation","source_fingerprint_sha256_160bits"])
    w.writerows(native_status)
native_counts={k:sum(x[4]==k for x in native_status) for k in
  ("PASS_NATIVE_BBJ_ALLELE2_BETA_EAF_GRCH38_CHAIN",
   "NOT_FOUND_IN_EXTENDED_NATIVE_REFERENCE","NATIVE_EFFECT_OR_CHAIN_MISMATCH")}
print("NATIVE_BBJ_PRIMARY_SOURCE_AND_LIFTOVER_COUNTS=",native_counts)
print("NATIVE_GWAS_SNP_CATALOG_SHA256=",sha(SNP_CATALOG))
if native_counts["PASS_NATIVE_BBJ_ALLELE2_BETA_EAF_GRCH38_CHAIN"]!=8956:
    failure_output=GDRIVE/"results/IS_646_GWAS_NATIVE_PREFLIGHT_FAILURE_V2"
    if failure_output.exists():
        raise FileExistsError("NATIVE_PREFLIGHT_OUTPUT_EXISTS_REQUIRES_VERSIONING")
    failure_output.mkdir(parents=True)
    shutil.copyfile(SNP_CATALOG,failure_output/SNP_CATALOG.name)
    raise RuntimeError("FULL_646_NATIVE_BBJ_OR_CHAIN_INCOMPLETE; audit: "+str(failure_output))
print("FULL_8956_NATIVE_BBJ_ORIGINAL_EFFECT_ALLELE_AND_LIFTOVER_VERIFIED")
'''
    py=py.replace(py[py.index(start):py.index(stop)],preflight)
    # Preserve completed original 646 G0/G1 and extreme MIN-N sensitivity untouched.
    py=py.replace('is_phase10_11_646_variant_n_sdy','is_phase10_11_646_variant_n_sdy_v2')
    old='IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V1'
    assert py.count(old)==1
    py=py.replace(old,'IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V2')
    old='shutil.copyfile(SNP_CATALOG,OUTPUT/"IS_FOUR_BBJ_LOCI_NATIVE_VERIFIED_GWAS_CATALOG.tsv")'
    assert py.count(old)==1
    py=py.replace(old,'shutil.copyfile(SNP_CATALOG,OUTPUT/SNP_CATALOG.name)')
    old='"native_BBJ_8309_catalog_sha256":sha(SNP_CATALOG),'
    assert py.count(old)==1
    py=py.replace(old,'"native_BBJ_8956_validation_status":"PASS_NATIVE_BBJ_ALLELE2_BETA_EAF_GRCH38_CHAIN",\n    "native_BBJ_reference_fingerprint_SHA256":NATIVE_REFERENCE_SHA256,\n    "native_BBJ_full_8956_status_SHA256":sha(SNP_CATALOG),\n    "native_BBJ_30_previous_union":len(old_keys),\n    "native_BBJ_new_variants_audited":len(set(combined)-old_keys),')
    old='"native_BBJ_8309_unique_variants_covered_by_all_646":len(unique_gwas),'
    assert py.count(old)==1
    py=py.replace(old,'"native_BBJ_unique_variants_full_646":len(combined),')
    py=py.replace('One-cell 646 historical BBJ x GTEx variant-specific QTL-N assisted common-sdY calibration.',
      'One-cell 646 original coloc QTL-N common-sdY and 8956 original BBJ native source-chain verification.')
    old_schema='"schema":"IS_PHASE10_11_646_VARIANT_N_SDY_COLAB_V1",'
    assert py.count(old_schema)==1, "V2_SCHEMA_EXPECTED_FROM_TEMPLATE"
    py=py.replace(old_schema,'"schema":"IS_PHASE10_11_646_VARIANT_N_SDY_COLAB_V2",')
    ast.parse(py)
    original["cells"][0]["source"]=py.splitlines(keepends=True)
    original["cells"][0]["outputs"]=[]
    original["cells"][0]["execution_count"]=None
    original["metadata"]["colab"]["name"]=out.name
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(original,ensure_ascii=False,indent=1)+"\n")
    assert "".join(json.loads(out.read_text())["cells"][0]["source"])==py
    print("REPAIR_NOTEBOOK",out.name,"EXPECTED_NATIVE_REFERENCE_SHA256",expected_sha,
        "SOURCE_NOTEBOOK_SHA256",hashlib.sha256(source.read_bytes()).hexdigest(),
        "OUTPUT_SHA256",hashlib.sha256(out.read_bytes()).hexdigest())
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--ref",type=Path,required=True)
    p.add_argument("--cache30",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();build(a.source,a.ref,a.cache30,a.out)
