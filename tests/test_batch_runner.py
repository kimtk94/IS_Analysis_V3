import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location("runner", ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py")
runner = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(runner)

class BatchRunnerTests(unittest.TestCase):
    def test_paired_batches_exclude_unpaired_and_are_stable(self):
        rows = runner.read_tsv(ROOT / "tests/fixtures/ukb_ppp_download_manifest.tsv")
        batches = runner.paired_batches(rows, 15)
        self.assertEqual([15, 1], list(map(len, batches)))
        self.assertNotIn("UNPAIRED", sum(batches, []))
        self.assertEqual(sorted(sum(batches, [])), sum(batches, []))

    def test_focus_dry_run_writes_limits_and_no_runtime_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); qc = root / "qc"; raw = root / "raw"
            command = [sys.executable, str(ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py"),
                "--base", str(raw), "--qc-dir", str(qc), "--outdir", str(root/"out"),
                "--standardized-dir", str(root/"std"), "--instrument-dir", str(root/"inst"),
                "--download-manifest", str(ROOT/"tests/fixtures/ukb_ppp_download_manifest.tsv"),
                "--gene-coordinate-file", str(ROOT/"tests/fixtures/gene_coordinates_hg38.tsv"),
                "--batch-size", "15", "--focus-gene", "IDO1", "--focus-max-bytes", "20000000",
                "--other-max-file-lines", "1000"]
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertFalse(raw.exists())
            with (qc/"execution_plan.tsv").open() as f: plan=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(30, len(plan))
            ido1=[r for r in plan if r["gene"] == "IDO1"]
            self.assertEqual({("bytes","20000000")}, {(r["limit_type"],r["limit"]) for r in ido1})
            self.assertTrue(all(int(r["limit"]) <= 1000 for r in plan if r["gene"] != "IDO1"))

    def test_download_only_downloads_focused_batch_without_running_r(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["gene", "ancestry", "source_url", "sha256", "size_bytes"], delimiter="\t")
                writer.writeheader()
                fixtures = {"ACE": ROOT / "tests/fixtures/gigastroke.tsv",
                            "IDO1": ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"}
                for gene, source in fixtures.items():
                    for ancestry in ("EUR", "EAS"):
                        writer.writerow({"gene": gene, "ancestry": ancestry,
                                         "source_url": str(source), "sha256": "",
                                         "size_bytes": source.stat().st_size})

            raw = root / "raw"
            qc = root / "qc"
            command = [sys.executable, str(ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py"),
                "--base", str(raw), "--qc-dir", str(qc), "--outdir", str(root/"out"),
                "--standardized-dir", str(root/"std"), "--instrument-dir", str(root/"inst"),
                "--download-manifest", str(manifest),
                "--gene-coordinate-file", str(ROOT/"tests/fixtures/gene_coordinates_hg38.tsv"),
                "--batch-size", "2", "--focus-gene", "IDO1", "--download-only"]
            subprocess.run(command, check=True, capture_output=True, text=True)

            self.assertEqual(4, len(list(raw.glob("*/*.tsv"))))
            self.assertFalse((root / "std").exists())
            with (qc / "batch_progress.tsv").open() as handle:
                progress = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({"downloaded"}, {row["status"] for row in progress})

if __name__ == "__main__": unittest.main()
