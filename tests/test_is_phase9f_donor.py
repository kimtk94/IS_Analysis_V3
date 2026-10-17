"""Regression guards for Phase9F R3_7 donor-aware pseudobulk and notebook generator."""
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "scripts/is/build_phase9f_r3_7_donor_colab.py"
spec = importlib.util.spec_from_file_location("build_phase9f_r3_7",CODE)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


class Phase9FDonorTests(unittest.TestCase):
    def test_signature_fail_closed(self):
        with self.assertRaises(ValueError):
            builder.replace_once("xx aa xx aa", "aa", "bb")
        with self.assertRaises(ValueError):
            builder.replace_once("nothing", "aa", "bb")
        self.assertEqual(builder.replace_once("one aa two","aa","bb"),"one bb two")

    def test_extract_r_source_and_python_parse(self):
        self.assertEqual(builder.extract_embedded_r("r_code=r'''x<-1\\n'''\n"),"x<-1\\n")
        with self.assertRaises(ValueError):
            builder.extract_embedded_r("something = 1\n")

    def test_missing_base_fails_without_creating_output(self):
        with tempfile.TemporaryDirectory() as t:
            out = Path(t) / "x.ipynb"
            with self.assertRaises(FileNotFoundError):
                builder.build_notebook(Path(t)/"absent.ipynb",
                                       ROOT/"scripts/is/phase9f_donor_pseudobulk.R",
                                       out,rscript="")
            self.assertFalse(out.exists())

    @unittest.skipUnless(shutil.which("Rscript"),"Rscript not installed in test environment")
    def test_donor_helper_synthetic_r_fixture(self):
        helper=ROOT/"scripts/is/phase9f_donor_pseudobulk.R"
        fixture=ROOT/"tests/fixtures/is_phase9f_donor_pseudobulk_selftest.R"
        proc=subprocess.run(["Rscript",str(fixture),str(helper)],
                            capture_output=True,text=True,timeout=30)
        if proc.returncode and "there is no package called" in proc.stderr:
            self.skipTest("R Matrix package not installed")
        self.assertEqual(proc.returncode,0,msg=proc.stdout+"\n"+proc.stderr)
        self.assertIn("IS_PHASE9F_DONOR_PSEUDOBULK_SYNTHETIC_TEST_PASS",proc.stdout)


if __name__=="__main__":
    unittest.main()
