import csv
import gzip
import hashlib
import importlib.util
import io
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest import mock
from pathlib import Path

ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location("runner", ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py")
runner = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(runner)

class BatchRunnerTests(unittest.TestCase):
    def test_write_test_sample_uncompresses_plain_gzip_and_tar_member_gzip(self):
        content = b"variant\tbeta\nrs1\t0.1\nrs2\t0.2\n"
        expected_rows = [{"variant": "rs1", "beta": "0.1"},
                         {"variant": "rs2", "beta": "0.2"}]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plain = root / "input.tsv"
            standalone = root / "input.tsv.gz"
            archive = root / "input.tar"
            plain.write_bytes(content)
            with gzip.open(standalone, "wb") as handle:
                handle.write(content)
            compressed = gzip.compress(content)
            with tarfile.open(archive, "w") as handle:
                info = tarfile.TarInfo("nested/input.tsv.gz")
                info.size = len(compressed)
                handle.addfile(info, io.BytesIO(compressed))

            for number, source in enumerate((plain, standalone, archive)):
                with self.subTest(source=source.name):
                    destination = root / f"sample-{number}.tsv"
                    runner.write_test_sample(source, destination, "lines", 3)

                    self.assertEqual(".tsv", destination.suffix)
                    self.assertNotEqual(b"\x1f\x8b", destination.read_bytes()[:2])
                    with destination.open(newline="", encoding="utf-8") as handle:
                        reader = csv.DictReader(handle, delimiter="\t")
                        self.assertEqual(["variant", "beta"], reader.fieldnames)
                        self.assertEqual(expected_rows, list(reader))

    def test_stage_reuses_existing_file_when_all_manifest_values_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "download.tsv"
            destination.write_bytes(b"verified fixture\n")
            row = {"gene": "IDO1", "ancestry": "EUR", "source_url": "https://invalid.test/file",
                   "size_bytes": str(destination.stat().st_size),
                   "sha256": runner.sha256(destination), "md5": runner.md5(destination)}

            with mock.patch.object(runner.urllib.request, "urlretrieve") as download, \
                    mock.patch.object(runner, "file_digests",
                                      wraps=runner.file_digests) as calculate:
                self.assertEqual(destination, runner.stage(row, destination))

            download.assert_not_called()
            calculate.assert_called_once_with(destination, ["sha256", "md5"])

    def test_file_digests_calculates_sha256_and_md5_together(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.tsv"
            content = b"simultaneously verified fixture\n"
            source.write_bytes(content)

            self.assertEqual({
                "sha256": hashlib.sha256(content).hexdigest(),
                "md5": hashlib.md5(content, usedforsecurity=False).hexdigest(),
            }, runner.file_digests(source, ["sha256", "md5"]))

    def test_stage_rejects_sha256_mismatch_when_md5_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.tsv"
            destination = Path(tmp) / "download.tsv"
            content = b"dual digest fixture\n"
            source.write_bytes(content)
            row = {"gene": "IDO1", "ancestry": "EUR", "source_url": str(source),
                   "size_bytes": str(len(content)), "sha256": "0" * 64,
                   "md5": hashlib.md5(content, usedforsecurity=False).hexdigest()}

            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                runner.stage(row, destination)

            self.assertFalse(destination.exists())

    def test_stage_rejects_md5_mismatch_when_sha256_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.tsv"
            destination = Path(tmp) / "download.tsv"
            content = b"dual digest fixture\n"
            source.write_bytes(content)
            row = {"gene": "IDO1", "ancestry": "EAS", "source_url": str(source),
                   "size_bytes": str(len(content)),
                   "sha256": hashlib.sha256(content).hexdigest(), "md5": "0" * 32}

            with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
                runner.stage(row, destination)

            self.assertFalse(destination.exists())

    def test_stage_redownloads_checksum_mismatch_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "download.tsv"
            expected = b"replacement fixture\n"
            destination.write_bytes(b"x" * len(expected))
            expected_sha256 = hashlib.sha256(expected).hexdigest()
            row = {"gene": "IDO1", "ancestry": "EAS", "source_url": "https://example.test/file",
                   "size_bytes": str(len(expected)), "sha256": expected_sha256, "md5": ""}

            def download(_url, path):
                Path(path).write_bytes(expected)

            with mock.patch.object(runner.urllib.request, "urlretrieve",
                                   side_effect=download) as mocked_download:
                runner.stage(row, destination)

            mocked_download.assert_called_once_with(row["source_url"], destination)
            self.assertEqual(expected, destination.read_bytes())

    def test_stage_redownloads_unverified_file_unless_opted_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "download.tsv"
            destination.write_bytes(b"existing\n")
            row = {"gene": "IDO1", "ancestry": "EUR", "source_url": "https://example.test/file",
                   "size_bytes": "", "sha256": "", "md5": ""}

            def download(_url, path):
                Path(path).write_bytes(b"fresh\n")

            with mock.patch.object(runner.urllib.request, "urlretrieve",
                                   side_effect=download) as mocked_download:
                runner.stage(row, destination)
            mocked_download.assert_called_once()
            self.assertEqual(b"fresh\n", destination.read_bytes())

            with mock.patch.object(runner.urllib.request, "urlretrieve") as mocked_download:
                runner.stage(row, destination, reuse_unverified=True)
            mocked_download.assert_not_called()
            self.assertEqual(b"fresh\n", destination.read_bytes())

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
                "--batch-size", "15", "--focus-gene", "IDO1",
                "--other-max-file-lines", "1000"]
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertFalse(raw.exists())
            with (qc/"execution_plan.tsv").open() as f: plan=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(30, len(plan))
            ido1=[r for r in plan if r["gene"] == "IDO1"]
            self.assertEqual({("lines", "500000")},
                             {(r["limit_type"], r["limit"]) for r in ido1})
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

    def test_test_data_mode_applies_line_limits(self):
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

    def test_write_test_sample_byte_limit_discards_partial_final_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"
            destination = Path(tmp) / "sample.tsv"

            runner.write_test_sample(source, destination, "bytes", 80)

            sample = destination.read_bytes()
            self.assertLessEqual(len(sample), 80)
            self.assertTrue(sample.endswith(b"\n"))
            self.assertEqual(source.read_bytes().splitlines(keepends=True)[:2],
                             sample.splitlines(keepends=True))

    def test_write_test_sample_line_limit_stops_at_requested_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"
            destination = Path(tmp) / "sample.tsv"

            runner.write_test_sample(source, destination, "lines", 2)

            self.assertEqual(source.read_bytes().splitlines(keepends=True)[:2],
                             destination.read_bytes().splitlines(keepends=True))

    def test_write_test_sample_removes_partial_destination_on_no_data_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = ROOT / "tests/fixtures/gigastroke.tsv"
            destination = Path(tmp) / "sample.tsv"

            with self.assertRaisesRegex(ValueError, "sample contains no complete data rows"):
                runner.write_test_sample(source, destination, "bytes", 50)

            self.assertFalse(destination.exists())

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
                                         "source_file": f"IDO1_{ancestry}_assay_{number}.tsv",
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

    def test_test_data_mode_reuses_existing_unverified_raw_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            source = ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"
            raw = root / "raw"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                fields = ["gene", "ancestry", "source_url", "source_file"]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for ancestry in ("EUR", "EAS"):
                    source_file = f"IDO1_{ancestry}.tsv"
                    writer.writerow({"gene": "IDO1", "ancestry": ancestry,
                                     "source_url": "/missing/source.tsv",
                                     "source_file": source_file})
                    destination = raw / ancestry / source_file
                    destination.parent.mkdir(parents=True)
                    destination.write_bytes(source.read_bytes())

            command = [sys.executable, str(ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py"),
                "--base", str(raw), "--qc-dir", str(root / "qc"),
                "--outdir", str(root / "out"), "--standardized-dir", str(root / "std"),
                "--instrument-dir", str(root / "inst"), "--download-manifest", str(manifest),
                "--gene-coordinate-file", str(ROOT / "tests/fixtures/gene_coordinates_hg38.tsv"),
                "--batch-size", "1", "--focus-gene", "IDO1", "--focus-max-file-lines", "2",
                "--test-data-only"]
            subprocess.run(command, check=True, capture_output=True, text=True, cwd=ROOT)

            samples = sorted((root / "out/test_data").glob("*/*.tsv"))
            self.assertEqual(2, len(samples))
            expected = b"".join(source.read_bytes().splitlines(keepends=True)[:2])
            self.assertTrue(all(sample.read_bytes() == expected for sample in samples))
            with (root / "qc/batch_progress.tsv").open() as handle:
                progress = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({"sampled"}, {row["status"] for row in progress})

    def test_write_test_sample_refuses_to_overwrite_existing_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "sample.tsv"
            destination.write_text("existing\ncontent\n", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "refusing to overwrite"):
                runner.write_test_sample(ROOT / "tests/fixtures/gigastroke.tsv",
                                         destination, "lines", 2)
            self.assertEqual("existing\ncontent\n", destination.read_text())

if __name__ == "__main__": unittest.main()
