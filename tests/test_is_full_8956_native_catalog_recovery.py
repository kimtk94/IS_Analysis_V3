"""Regression of 8309 local-cache native variants vs 8956 full historical universe.

Protect source QTL inputs and audit lineage. The extra 647 SNPs are not proven
until full 646 source Colab validates each against the BBJ original ZIP chain
fingerprint reference. Test local native fingerprints against cached 8309 only.
"""
import ast,csv,gzip,hashlib,importlib.util,json,tempfile,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/"data/is_reference/native_bbj_is_4loci_extended_20261011_v2.fingerprints.tsv.gz"
META=ROOT/"data/is_reference/native_bbj_is_4loci_extended_20261011_v2.manifest.json"
FROZEN=ROOT/"data/is_reference/native_verified_original_30_input_SHA256.json"
BUILDER=ROOT/"scripts/is/repair_is_646_variant_n_native_reference_colab.py"
SOURCE=ROOT/"notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_ONECELL.ipynb"
TARGET=ROOT/"notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_NATIVE_8956_REPAIRED_ONECELL.ipynb"
DATA=Path("/srv/is-analysis/results/is/audits/legacy_coloc_snp_replay_20261010_v1/input")

def module():
 spec=importlib.util.spec_from_file_location("is_g1_native_8956_builder",BUILDER)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 return mod

def fp(z):
 normal="IS_BBJ_NATIVE_4LOCI_V1\t"+"\t".join((z["locus"],z["variant_id"],
      z["match_key"],format(float(z["gwas_beta"]),".12f"),
      format(float(z["gwas_eaf"]),".12f")))
 return hashlib.sha256(normal.encode()).hexdigest()[:40]

class TestFullNativeAudit(unittest.TestCase):
 def test_reference_manifest_integrity_and_no_raw_GWAS_stat_columns(self):
  manifest=json.loads(META.read_text())
  self.assertEqual(manifest["reference_variants"],56189)
  self.assertEqual(hashlib.sha256(REF.read_bytes()).hexdigest(),manifest["gz_sha256"])
  with gzip.open(REF,"rt") as h:
   reader=csv.DictReader(h,delimiter="\t")
   self.assertEqual(reader.fieldnames,
       ["locus","variant_id","native_variant_fingerprint_sha256_160bits"])
   n=0
   for z in reader:
    self.assertEqual(len(z["native_variant_fingerprint_sha256_160bits"]),40)
    n+=1
   self.assertEqual(n,56189)
  self.assertNotIn("gwas_beta",gzip.decompress(REF.read_bytes())[:120].decode())
  self.assertNotIn("gwas_eaf",gzip.decompress(REF.read_bytes())[:120].decode())

 def test_all_8309_prior_source_SNPs_have_original_BBJ_native_fingerprint(self):
  with gzip.open(REF,"rt") as h:
   ref={(z["locus"],z["variant_id"]):z["native_variant_fingerprint_sha256_160bits"]
      for z in csv.DictReader(h,delimiter="\t")}
  self.assertEqual(len(ref),56189)
  cached=json.loads(FROZEN.read_text())
  self.assertEqual(len(cached),30)
  original={}
  for filename,digest in cached.items():
   path=DATA/filename
   self.assertTrue(path.is_file())
   self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),digest)
   with path.open(encoding="utf8",newline="") as handle:
    for z in csv.DictReader(handle,delimiter="\t"):
     key=(z["locus"],z["variant_id"])
     value=fp(z)
     self.assertEqual(value,ref.get(key))
     self.assertEqual(original.get(key,value),value)
     original[key]=value
  self.assertEqual(len(original),8309)

 def test_repaired_onecell_has_no_8309_full_universe_constraint(self):
  nb=json.loads(TARGET.read_text())
  self.assertEqual(len(nb["cells"]),1)
  py="".join(nb["cells"][0]["source"])
  ast.parse(py)
  self.assertIn("assert len(combined)==8956",py)
  self.assertIn("assert len(old_keys)==8309",py)
  self.assertIn("assert len(set(combined)-old_keys)==647",py)
  self.assertIn("FULL_8956_NATIVE_BBJ_ORIGINAL_EFFECT_ALLELE_AND_LIFTOVER_VERIFIED",py)
  self.assertIn("ORIGINAL_CACHED_30_FILE_SHA256_VERIFIED",py)
  self.assertIn("NATIVE_EFFECT_OR_CHAIN_MISMATCH",py)
  self.assertIn("NATIVE_PREFLIGHT_OUTPUT_EXISTS_REQUIRES_VERSIONING",py)
  self.assertIn("IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V2",py)
  self.assertIn("FULL_646_VARIANT_N_SDY_CALIBRATION_COMPLETE",py)
  self.assertIn("original_index",py)
  self.assertIn(json.loads(META.read_text())["gz_sha256"],py)
  self.assertNotIn('assert len(unique_gwas)==8309',py)
  self.assertIn("BASELINE_SDY_EXPLICIT_REPRODUCTION_FAILURE",py)
  self.assertIn("if OUTPUT.exists()",py)

 def test_builder_reproducible_source_and_refuse_overwrite(self):
  with tempfile.TemporaryDirectory() as t:
   dst=Path(t)/TARGET.name
   module().build(SOURCE,REF,FROZEN,dst)
   self.assertEqual(dst.read_bytes(),TARGET.read_bytes())
   with self.assertRaises(FileExistsError):
    module().build(SOURCE,REF,FROZEN,dst)
   self.assertEqual(dst.read_bytes(),TARGET.read_bytes())

if __name__=="__main__":
 unittest.main()
