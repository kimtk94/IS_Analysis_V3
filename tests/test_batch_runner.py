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
    def test_read_tsv_normalizes_production_synapse_schema(self):
        rows = runner.read_tsv(ROOT / "tests/fixtures/ukb_ppp_download_manifest_synapse.tsv")
        self.assertEqual("IDO1", rows[0]["gene"])
        self.assertEqual("https://www.synapse.org/Synapse:syn00000001", rows[0]["source_url"])
        self.assertEqual("12", rows[0]["size_bytes"])
        self.assertEqual("IDO1_EUR.tar", rows[0]["source_file"])

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

    def test_duplicate_gene_ancestry_sources_are_each_staged_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            qc = root / "qc"
            command = [sys.executable, str(ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py"),
                "--base", str(raw), "--qc-dir", str(qc), "--outdir", str(root / "out"),
                "--standardized-dir", str(root / "std"), "--instrument-dir", str(root / "inst"),
                "--download-manifest", str(ROOT / "tests/fixtures/ukb_ppp_duplicate_assays_manifest.tsv"),
                "--gene-coordinate-file", str(ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"),
                "--batch-size", "1", "--focus-gene", "IDO1", "--download-only"]
            subprocess.run(command, check=True, capture_output=True, text=True, cwd=ROOT)

            staged = sorted(path.name for path in raw.glob("*/*"))
            self.assertEqual([
                "IDO1_EAS_assay_a.tsv", "IDO1_EAS_assay_b.tsv",
                "IDO1_EUR_assay_a.tsv", "IDO1_EUR_assay_b.tsv"], staged)
            with (qc / "batch_progress.tsv").open() as handle:
                progress = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(4, len(progress))
            self.assertEqual(4, len({row["source_key"] for row in progress}))
            self.assertEqual(set(staged), {row["source_file"] for row in progress})
            self.assertEqual({"downloaded"}, {row["status"] for row in progress})

    def test_test_data_mode_applies_focus_bytes_and_other_line_limits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"
            target = root / "target.tsv"
            other = root / "other.tsv"

            runner.write_test_sample(source, target, "bytes", 80)
            runner.write_test_sample(source, other, "lines", 2)

            self.assertLessEqual(target.stat().st_size, 80)
            self.assertTrue(target.read_bytes().endswith(b"\n"))
            self.assertEqual(source.read_text().splitlines()[:2], other.read_text().splitlines())

    def test_test_data_mode_preserves_multiple_assays_for_gene_and_ancestry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            sources = [ROOT / "tests/fixtures/gene_coordinates_hg38.tsv",
                       ROOT / "tests/fixtures/gigastroke.tsv"]
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                fields = ["gene", "ancestry", "source_url", "source_file", "synapse_id"]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for ancestry in ("EUR", "EAS"):
                    for number, source in enumerate(sources, 1):
                        writer.writerow({"gene": "IDO1", "ancestry": ancestry,
                                         "source_url": str(source),
                                         "source_file": f"IDO1_{ancestry}_assay_{number}.tsv.tar.gz",
                                         "synapse_id": f"syn{ancestry}{number}"})

            command = [sys.executable, str(ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py"),
                "--base", str(root / "raw"), "--qc-dir", str(root / "qc"),
                "--outdir", str(root / "out"), "--standardized-dir", str(root / "std"),
                "--instrument-dir", str(root / "inst"), "--download-manifest", str(manifest),
                "--gene-coordinate-file", str(ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"),
                "--batch-size", "1", "--focus-gene", "IDO1", "--focus-max-bytes", "10000",
                "--test-data-only"]
            subprocess.run(command, check=True, capture_output=True, text=True, cwd=ROOT)

            samples = sorted((root / "out/test_data").glob("*/*.tsv"))
            self.assertEqual(4, len(samples))
            self.assertEqual({
                "IDO1_EAS_assay_1.tsv", "IDO1_EAS_assay_2.tsv",
                "IDO1_EUR_assay_1.tsv", "IDO1_EUR_assay_2.tsv",
            }, {sample.name for sample in samples})
            for sample in samples:
                expected = sources[int(sample.stem.rsplit("_", 1)[1]) - 1].read_bytes()
                self.assertEqual(expected, sample.read_bytes())

    def test_write_test_sample_refuses_to_overwrite_existing_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "sample.tsv"
            destination.write_text("existing\ncontent\n", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "refusing to overwrite"):
                runner.write_test_sample(ROOT / "tests/fixtures/gigastroke.tsv",
                                         destination, "lines", 2)
            self.assertEqual("existing\ncontent\n", destination.read_text())

if __name__ == "__main__": unittest.main()
