import csv, json, math, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).parents[1]

class WorkflowTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
    def test_ukb_ppp_exposure_conversion_and_qc(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            command=["Rscript", str(ROOT/"scripts/01_prepare_exposure_fast.R"),
                "--archive", str(ROOT/"tests/fixtures/ukb_ppp_summary.tsv"),
                "--gene", "IDO1", "--ancestry", "EUR",
                "--coordinates", str(ROOT/"tests/fixtures/gene_coordinates_hg38.tsv"),
                "--standardized-dir", str(root/"standardized"),
                "--instrument-dir", str(root/"instruments"),
                "--legacy-dir", str(root/"legacy"), "--limit-type", "lines",
                "--limit", "20", "--p-value-threshold", "1e-7",
                "--f-statistic-threshold", "10", "--cis-window-bp", "1000000"]
            subprocess.run(command, cwd=ROOT, check=True)

            output=root/"standardized/EUR/batch_001/IDO1.tsv"
            with output.open() as handle:
                rows=list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(["rs_effect", "rs_zero"], [row["rsid"] for row in rows])
            self.assertEqual(("G", "A"), (rows[0]["effect_allele"], rows[0]["other_allele"]))
            self.assertEqual(("8", "39910000"), (rows[0]["chr"], rows[0]["pos"]))
            self.assertTrue(math.isclose(float(rows[0]["eaf"]), 0.25))
            self.assertTrue(math.isclose(float(rows[0]["p_value"]), 1e-8, rel_tol=1e-12))
            self.assertTrue(math.isclose(float(rows[0]["f_statistic"]), 16.0))
            self.assertTrue(math.isclose(float(rows[1]["p_value"]), 1.0))
            self.assertTrue(math.isclose(float(rows[1]["f_statistic"]), 9.0))

            with (root/"standardized/EUR/batch_001/IDO1.qc.tsv").open() as handle:
                qc=next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual("LOG10P", qc["p_value_source"])
            self.assertEqual("7", qc["input_rows"]); self.assertEqual("2", qc["output_rows"])
            self.assertEqual("5", qc["invalid_any"]); self.assertEqual("3", qc["invalid_p_value"])
            self.assertEqual("1", qc["invalid_eaf"]); self.assertEqual("1", qc["invalid_se"])
            self.assertEqual("1", qc["log10p_zero"]); self.assertEqual("1", qc["log10p_negative"])
            self.assertEqual("2", qc["log10p_nonfinite_or_missing"])

            with (root/"instruments/IDO1.tsv").open() as handle:
                instruments=list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(["rs_effect"], [row["rsid"] for row in instruments])

    def test_gigastroke_alias_adapter(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"out.tsv"
            subprocess.run([sys.executable,str(ROOT/"workflow/gigastroke_outcome_adapter.py"),"--input",str(ROOT/"tests/fixtures/gigastroke.tsv"),"--output",str(out)],cwd=ROOT,check=True)
            with out.open() as f: rows=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual("39910000",rows[0]["pos"]); self.assertEqual("0.12",rows[0]["beta"])

    def test_checkpoint_contains_input_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"checkpoint.json"; fixture=ROOT/"tests/fixtures/gigastroke.tsv"
            subprocess.run([sys.executable,str(ROOT/"workflow/causal_checkpoint_analysis.py"),"--exposure",str(fixture),"--outcome",str(fixture),"--ancestry","EUR","--output",str(out)],cwd=ROOT,check=True)
            payload=json.loads(out.read_text()); self.assertEqual("validated",payload["status"])
            self.assertEqual(64,len(payload["inputs"][0]["sha256"]))

if __name__ == "__main__": unittest.main()
