"""Native BBJ 8956 and 646 sdY Colab V2 output release gate tests."""
import csv,gzip,hashlib,importlib.util,json,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"scripts/is/audit_is_646_variant_n_sdy_v2_results.py"
NOTE=ROOT/"notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_NATIVE_8956_REPAIRED_ONECELL.ipynb"
REF=ROOT/"data/is_reference/native_bbj_is_4loci_extended_20261011_v2.fingerprints.tsv.gz"
G0=Path("/srv/is-analysis/results/is/audits/legacy_coloc_646_full_colab_20261011_v1/source")
PILOT=Path("/srv/is-analysis/results/is/audits/is_646_variant_n_sdy_cache30_20261011_v2")
INDEX=Path("/srv/is-analysis/results/is/audits/legacy_coloc_snp_replay_20261010_v1/source_reference/COLOC_ABF_INPUT_INDEX_V2.tsv")
MASTER=Path("/srv/is-analysis/results/is/audits/legacy_coloc_snp_replay_20261010_v1/source_reference/COLOC_ABF_MASTER_V2.tsv")
RUNNER=ROOT/"scripts/is/is_646_qtl_variant_n_sdy_calibration.R"

def audit():
    spec=importlib.util.spec_from_file_location("native_v2_auditor",SRC)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

class TestV2Release(unittest.TestCase):
    def test_onecell_v2_pins_full_native_and_source(self):
        nb=json.loads(NOTE.read_text())
        self.assertEqual(len(nb["cells"]),1)
        code="".join(nb["cells"][0]["source"])
        self.assertIn('"schema":"IS_PHASE10_11_646_VARIANT_N_SDY_COLAB_V2"',code)
        for marker in ("NATIVE_REFERENCE_TAMPER_OR_URL_MISMATCH",
                       "GWAS_SOURCE_CONFLICT_ACROSS_GENE_TISSUE_ASSAYS",
                       "NATIVE_EFFECT_OR_CHAIN_MISMATCH",
                       "assert len(combined)==8956","if OUTPUT.exists():"):
            self.assertIn(marker,code)
        self.assertIn(hashlib.sha256(REF.read_bytes()).hexdigest(),code)

    def test_audit_does_not_overwrite_existing(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/"audit";out.mkdir()
            sentinel=out/"before.txt";sentinel.write_text("KEEP")
            with self.assertRaises(FileExistsError):
                audit().check(Path(t),G0,REF,INDEX,MASTER,RUNNER,out)
            self.assertEqual(sentinel.read_text(),"KEEP")

    def test_missing_full_source_fails(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/"new"
            with self.assertRaisesRegex(FileNotFoundError,"MISSING_V2_SOURCE"):
                audit().check(Path(t),G0,REF,INDEX,MASTER,RUNNER,out)
            self.assertFalse(out.exists())

    def test_existing_30_real_sdy_posterior_matches_original(self):
        mod=audit()
        g0=mod.keys(mod.tsv(G0/"IS_LEGACY_646_SNP_REPLAY_STATUS.tsv"),"g0")
        local=[z for z in mod.tsv(PILOT/"IS_VARIANT_N_SDY_CALIBRATION_STATUS_646.tsv") if z["status"]=="PASS"]
        self.assertEqual(len(local),30)
        for z in local:
            k=tuple(z[x] for x in mod.KEY)
            self.assertAlmostEqual(float(z["H3_original"]),float(g0[k]["PP_H3"]),delta=1e-8)
            self.assertAlmostEqual(float(z["H4_original"]),float(g0[k]["PP_H4"]),delta=1e-8)
            self.assertAlmostEqual(float(z["H4_delta"]),float(z["H4_variantN_sdY"])-float(z["H4_original"]),delta=1e-10)

    def test_native_reference_non_raw_56189(self):
        with gzip.open(REF,"rt") as f:
            reader=csv.DictReader(f,delimiter="\t")
            self.assertEqual(reader.fieldnames,["locus","variant_id","native_variant_fingerprint_sha256_160bits"])
            vals=[x["native_variant_fingerprint_sha256_160bits"] for x in reader]
        self.assertEqual(len(vals),56189)
        self.assertTrue(all(len(x)==40 for x in vals))

if __name__=="__main__":
    unittest.main()
