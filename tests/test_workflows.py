import csv, json, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).parents[1]

class WorkflowTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
    def test_ukb_ppp_delimiters_produce_same_canonical_schema(self):
        outputs = []
        with tempfile.TemporaryDirectory() as tmp:
            for fixture in ("ukb_ppp_exposure_tab.tsv", "ukb_ppp_exposure_spaces.txt"):
                root = Path(tmp) / fixture
                command = [
                    "Rscript", str(ROOT / "scripts/01_prepare_exposure_fast.R"),
                    "--archive", str(ROOT / "tests/fixtures" / fixture),
                    "--gene", "IDO1", "--ancestry", "EUR",
                    "--coordinates", str(ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"),
                    "--standardized-dir", str(root / "standardized"),
                    "--instrument-dir", str(root / "instruments"),
                    "--legacy-dir", str(root / "legacy"),
                    "--limit-type", "lines", "--limit", "100",
                ]
                subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
                output = root / "standardized/EUR/batch_001/IDO1.tsv"
                with output.open() as stream:
                    outputs.append(list(csv.DictReader(stream, delimiter="\t")))

        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(
            ["gene", "ancestry", "chr", "pos", "rsid", "effect_allele",
             "other_allele", "beta", "se", "p_value", "eaf", "f_statistic"],
            list(outputs[0][0]),
        )

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
