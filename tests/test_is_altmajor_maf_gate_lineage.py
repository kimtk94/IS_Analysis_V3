"""Enforce native-to-reference genotype --maf QC causal chronology (not gene causality)."""
import importlib.util,sys,tempfile,unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_altmajor_maf_gate_lineage.py"
spec=importlib.util.spec_from_file_location("maf_lineage",P)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

class TestMaFGateLineage(unittest.TestCase):
    def setup_fixture(self,td):
        base=Path(td)
        raw=base/"1kg.pvar"
        qc=base/"1kg.QC.pvar"
        header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        native="12\t112090022\t12:112090022:C:A\tC\tA\t100\tPASS\tEAS_AF=0.999\n"
        other="12\t112094444\t12:112094444:G:A\tG\tA\t100\tPASS\tEAS_AF=0.12\n"
        raw.write_text(header+native+other)
        qc.write_text(header+other)
        ten=[{
            "site":f"12:{111652218+i}:C:A",
            "native_status":"NOT_REPORTED_AT_GRCH37_POSITION",
            "native_variant_candidates":"",
        } for i in range(9)]
        ten.append({
          "site":"12:111652218:C:A",
          "native_status":"BBJ_ORIGINAL_SAME_REF_ALT_PAIR",
          "native_variant_candidates":"12:112090022:C:A",
          "native_allele2_AF":"0.997129291706149","native_p":"0.161785997450046",
          "raw_QTD_ALT_AF_min":"0.826873","raw_QTD_ALT_AF_max":"0.864726"
        })
        ten[:9]=[{"site":f"12:{111739971+i}:T:C",**{k:v for k,v in r.items() if k!="site"}} for i,r in enumerate(ten[:9])]
        log=("Options in effect:\n  --geno 0.05\n  --maf 0.01\n"
             "0 variants removed due to missing genotype data.\n"
             "64428 variants removed due to allele frequency threshold(s)\n")
        stage3=[{"variant_id":"12:112090111:G:C","maf":"0.0100043"}]
        return raw,qc,ten,log,stage3

    def test_exact_maf_qc_lineage(self):
        with tempfile.TemporaryDirectory() as tmp:
            src,qc,ten,log,susie=self.setup_fixture(tmp)
            out=mod.audit(ten,src,qc,log,susie)
            self.assertTrue(out["1kg_EAS_post_QC_variant_absent"])
            self.assertEqual(out["n_native_BBJ_GWAS_absent"],9)
            self.assertAlmostEqual(out["1kg_EAS_MAF_derived"],0.001)
            self.assertEqual(out["PLINK_filter_maf_threshold"],0.01)

    def test_maf_invalidates_causal_exclusion_inference(self):
        with tempfile.TemporaryDirectory() as tmp:
            src,qc,ten,log,susie=self.setup_fixture(tmp)
            src.write_text(src.read_text().replace("EAS_AF=0.999","EAS_AF=0.96"))
            with self.assertRaises(ValueError):mod.audit(ten,src,qc,log,susie)

    def test_variant_after_qc_invalidates_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            src,qc,ten,log,susie=self.setup_fixture(tmp)
            qc.write_text(src.read_text())
            with self.assertRaises(ValueError):mod.audit(ten,src,qc,log,susie)

    def test_no_explicit_qc_maf_threshold_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            src,qc,ten,log,susie=self.setup_fixture(tmp)
            with self.assertRaises(ValueError):
                mod.audit(ten,src,qc,log.replace("--maf 0.01","--maf 0.001"),susie)
if __name__=="__main__":
    unittest.main()
