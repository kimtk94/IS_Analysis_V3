"""Regression tests for BBJ native Allele2, beta/EAF and chain audit guards."""
import importlib.util,sys,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/is/audit_all_8956_bbj_native_allele_chain.py"
spec=importlib.util.spec_from_file_location("bbj8956",P)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

class LiftPlus:
    def convert_coordinate(self,chrom,pos):
        return [(chrom,pos-437804,"+",1)]
class LiftWrong:
    def convert_coordinate(self,chrom,pos):
        return [(chrom,pos-437804,"-",1)]

class TestFullNativeAlleleAudit(unittest.TestCase):
    def fixture(self):
        r={
          "locus":"BBJ_IS_L003",
          "variant_id_GRCh37":"12:112241766:G:A",
          "match_key_GRCh38":"12:111803962:G:A",
          "gwas_beta_cached":-0.153,
          "gwas_eaf_cached":0.225,
          "ref_GRCh38":"G","alt_GRCh38":"A",
          "count_assay_rows":16,
          "per_assay_tissues":{"Artery_Aorta","Brain_Cortex"},
          "cache_consistency_status":"PASS",
        }
        n={"source_allele1":"G","source_allele2":"A","source_beta":-0.153,
           "source_AF_Allele2":0.225,
           "source_v":r["variant_id_GRCh37"],"source_position":"112241766",
           "source_chromosome":"12"}
        return {("BBJ_IS_L003",r["match_key_GRCh38"]):r},{r["variant_id_GRCh37"]:n}

    def test_native_allele2_and_chain_pass(self):
        sites,n=self.fixture()
        out=mod.audit(sites,n,set(),LiftPlus())[0]
        self.assertEqual(out["status"],"PASS_NATIVE_EFFECT_BETA_AF_AND_CHAIN")
        self.assertEqual(out["native_effect_vs_GRCh37_alt"],"ALT")
        self.assertTrue(out["GTEx_effect_allele_not_independently_checked"])

    def test_altered_native_gwas_beta_fails(self):
        sites,n=self.fixture()
        n["12:112241766:G:A"]["source_beta"]=-0.2
        out=mod.audit(sites,n,set(),LiftPlus())[0]
        self.assertEqual(out["status"],"REVIEW")
        self.assertIn("CACHED_GWAS_BETA_NOT_NATIVE",out["issues"])

    def test_missing_and_wrong_orientation_not_auto_flipped(self):
        sites,n=self.fixture()
        out=mod.audit(sites,{},set(),LiftPlus())[0]
        self.assertEqual(out["status"],"REVIEW")
        self.assertIn("NATIVE_VARIANT_NOT_FOUND",out["issues"])
        n["12:112241766:G:A"]["source_allele2"]="G"
        out=mod.audit(sites,n,set(),LiftPlus())[0]
        self.assertIn("NATIVE_ALLELE_REFALT_MISMATCH",out["issues"])

    def test_mismatched_chain_strand_rejected(self):
        sites,n=self.fixture()
        out=mod.audit(sites,n,set(),LiftWrong())[0]
        self.assertEqual(out["status"],"REVIEW")
        self.assertIn("CHAIN_BUILD_STRAND_OR_POS_MISMATCH",out["issues"])

if __name__=="__main__":
    unittest.main()
