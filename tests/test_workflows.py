import csv, gzip, io, json, shutil, subprocess, sys, tarfile, tempfile, unittest
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

@unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
class ExposurePreparationTests(unittest.TestCase):
    FIXTURES = ROOT / "tests/fixtures"

    def prepare(self, source, root, *, ancestry="EUR", extra=(), check=True):
        command = ["Rscript", str(ROOT / "scripts/01_prepare_exposure_fast.R"),
                   "--archive", str(source), "--gene", "IDO1", "--ancestry", ancestry,
                   "--coordinates", str(self.FIXTURES / "gene_coordinates_hg38.tsv"),
                   "--standardized-dir", str(root / "std"),
                   "--instrument-dir", str(root / "inst"), "--legacy-dir", str(root / "legacy"),
                   "--limit-type", "lines", "--limit", "100", *extra]
        return subprocess.run(command, cwd=ROOT, check=check, capture_output=True, text=True)

    @staticmethod
    def rows(path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def archive(self, path, members, *, mode="w"):
        with tarfile.open(path, mode) as handle:
            for member_name, fixture_name in members:
                content = (self.FIXTURES / fixture_name).read_bytes()
                if member_name.endswith(".gz"):
                    content = gzip.compress(content)
                info = tarfile.TarInfo(member_name)
                info.size = len(content)
                handle.addfile(info, io.BytesIO(content))

    def test_original_headers_log10p_alleles_and_direct_ancestry_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for ancestry in ("EUR", "EAS"):
                self.prepare(self.FIXTURES / "ukb_ppp_raw_ido1_assay_a.txt", root,
                             ancestry=ancestry)
                canonical = self.rows(root / "std" / ancestry / "batch_001/IDO1.tsv")
                instruments = self.rows(root / "inst" / ancestry / "batch_001/IDO1.tsv")
                self.assertEqual(["rs_assay_a", "rs_nonsignificant"],
                                 [row["rsid"] for row in canonical])
                self.assertEqual(["rs_assay_a"], [row["rsid"] for row in instruments])
                row = instruments[0]
                self.assertEqual(("8", "39910000", "A", "G"),
                                 tuple(row[k] for k in ("chr", "pos", "effect_allele", "other_allele")))
                self.assertAlmostEqual(1e-10, float(row["p_value"]), delta=1e-15)
                self.assertAlmostEqual(0.4, float(row["eaf"]))
                self.assertAlmostEqual(100, float(row["f_statistic"]))
                self.assertEqual(ancestry, row["ancestry"])
            self.assertFalse((root / "inst/IDO1.tsv").exists())
            self.assertEqual(2, len(list((root / "inst").glob("*/batch_001/IDO1.tsv"))))

    def test_log10p_conversion_and_qc_retain_pr18_numeric_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.prepare(self.FIXTURES / "ukb_ppp_summary.tsv", root)
            canonical = self.rows(root / "std/EUR/batch_001/IDO1.tsv")
            self.assertEqual(["rs_effect", "rs_zero"], [row["rsid"] for row in canonical])
            self.assertEqual(("G", "A"),
                             (canonical[0]["effect_allele"], canonical[0]["other_allele"]))
            self.assertAlmostEqual(1e-8, float(canonical[0]["p_value"]), delta=1e-15)
            self.assertEqual(1.0, float(canonical[1]["p_value"]))
            self.assertEqual(16.0, float(canonical[0]["f_statistic"]))
            qc = self.rows(root / "std/EUR/batch_001/IDO1.qc.tsv")[0]
            expected = {"input_rows": "7", "evaluated_rows": "7", "output_rows": "2",
                        "invalid_any": "5", "invalid_p_value": "3", "invalid_eaf": "1",
                        "invalid_se": "1", "invalid_f_statistic": "1", "log10p_zero": "1",
                        "log10p_negative": "1", "log10p_nonfinite_or_missing": "2",
                        "p_value_source": "LOG10P", "status": "PREPARED"}
            for key, value in expected.items():
                self.assertEqual(value, qc[key], key)
            self.assertEqual(["rs_effect"], [row["rsid"] for row in
                             self.rows(root / "inst/EUR/batch_001/IDO1.tsv")])

    def test_numeric_qc_counts_additive_gene_rows_and_rejects_nonfinite_f(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.prepare(self.FIXTURES / "ukb_ppp_numeric_qc.tsv", root)
            rows = self.rows(root / "std/EUR/batch_001/IDO1.tsv")
            self.assertEqual(["rs_underflow", "rs_frequency_one"], [row["rsid"] for row in rows])
            self.assertEqual([0.0, 1.0], [float(row["eaf"]) for row in rows])
            self.assertEqual(0.0, float(rows[0]["p_value"]))
            qc = self.rows(root / "std/EUR/batch_001/IDO1.qc.tsv")[0]
            expected = {"input_rows": "10", "non_additive_rows": "1",
                        "off_gene_chromosome_rows": "1", "evaluated_rows": "8",
                        "invalid_eaf": "3", "invalid_beta": "1", "invalid_coordinate": "1",
                        "invalid_f_statistic": "2", "invalid_any": "6", "output_rows": "2",
                        "instrument_rows": "2", "log10p_underflow": "1"}
            for key, value in expected.items():
                self.assertEqual(value, qc[key], key)

    def test_direct_p_precedence_and_optional_eaf(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.prepare(self.FIXTURES / "ukb_ppp_direct_p_optional_eaf.tsv", root)
            rows = self.rows(root / "std/EUR/batch_001/IDO1.tsv")
            self.assertEqual(["rs_direct_valid"], [row["rsid"] for row in rows])
            self.assertEqual("NA", rows[0]["eaf"])
            self.assertAlmostEqual(1e-9, float(rows[0]["p_value"]), delta=1e-15)
            qc = self.rows(root / "std/EUR/batch_001/IDO1.qc.tsv")[0]
            self.assertEqual("direct", qc["p_value_source"])
            self.assertEqual("FALSE", qc["eaf_supplied"])
            self.assertEqual("0", qc["invalid_eaf"])
            self.assertEqual("1", qc["invalid_p_value"])
            self.assertEqual("0", qc["log10p_nonfinite_or_missing"])

    def test_all_invalid_rows_fail_with_source_qc_and_no_data_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self.FIXTURES / "ukb_ppp_all_invalid.tsv"
            result = self.prepare(source, root, check=False,
                                  extra=("--batch-id", "batch_007", "--source-id", "assay_qc"))
            self.assertNotEqual(0, result.returncode)
            self.assertIn("no valid additive variants", result.stderr)
            path = root / "std/EUR/batch_007/IDO1/assay_qc.failed.qc.tsv"
            qc = self.rows(path)[0]
            self.assertEqual(("batch_007", "assay_qc", "NO_VALID_ROWS", "0"),
                             tuple(qc[key] for key in ("batch_id", "source_id", "status", "output_rows")))
            self.assertEqual(str(source.resolve()), qc["source_path"])
            self.assertFalse(path.with_name("assay_qc.tsv").exists())
            self.assertFalse(path.with_name("assay_qc.qc.tsv").exists())
            self.assertFalse((root / "inst").exists())
            self.assertFalse((root / "legacy").exists())

    def test_all_invalid_rerun_preserves_the_previous_successful_bundle(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "same_source.tsv"
            source.write_bytes((self.FIXTURES / "ukb_ppp_raw_ido1_assay_a.txt").read_bytes())
            identity = ("--batch-id", "batch_007", "--source-id", "assay_qc")
            self.prepare(source, root, extra=identity)
            relative = Path("EUR/batch_007/IDO1/assay_qc.tsv")
            successful = [root / part / relative for part in ("std", "inst", "legacy")]
            successful.append((root / "std" / relative).with_suffix(".qc.tsv"))
            before = {path: path.read_bytes() for path in successful}

            source.write_bytes((self.FIXTURES / "ukb_ppp_all_invalid.tsv").read_bytes())
            result = self.prepare(source, root, extra=identity, check=False)
            self.assertNotEqual(0, result.returncode)
            self.assertIn("failed-attempt QC", result.stderr)
            self.assertEqual(before, {path: path.read_bytes() for path in successful})
            failed = (root / "std" / relative).with_suffix(".failed.qc.tsv")
            qc = self.rows(failed)[0]
            self.assertEqual(("assay_qc", "NO_VALID_ROWS", "1", "0"),
                             tuple(qc[key] for key in ("source_id", "status", "invalid_any", "output_rows")))
            self.assertEqual("PREPARED", self.rows(successful[-1])[0]["status"])

    def test_tar_selects_gene_chromosome_and_decompresses_inner_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for suffix, mode in ((".tar", "w"), (".tar.gz", "w:gz")):
                with self.subTest(archive=suffix):
                    archive = root / ("source" + suffix)
                    self.archive(archive, [
                        ("nested/discovery_chr17_other.tsv.gz", "ukb_ppp_raw_ace.txt"),
                        ("nested/discovery_chr8_IDO1.tsv.gz", "ukb_ppp_raw_ido1_assay_a.txt"),
                    ], mode=mode)
                    output = root / suffix.replace(".", "_")
                    self.prepare(archive, output)
                    rows = self.rows(output / "std/EUR/batch_001/IDO1.tsv")
                    self.assertEqual({"8"}, {row["chr"] for row in rows})
                    self.assertEqual("rs_assay_a", rows[0]["rsid"])
                    qc = self.rows(output / "std/EUR/batch_001/IDO1.qc.tsv")[0]
                    self.assertEqual("nested/discovery_chr8_IDO1.tsv.gz", qc["archive_member"])

    def test_standalone_gzip_and_complete_line_byte_limit(self):
        content = (self.FIXTURES / "ukb_ppp_raw_ido1_assay_a.txt").read_bytes()
        lines = content.splitlines(keepends=True)
        limit = len(lines[0]) + len(lines[1]) + len(lines[2]) // 2
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            standalone = root / "summary.txt.gz"
            standalone.write_bytes(gzip.compress(content))
            archive = root / "summary.tar"
            self.archive(archive, [("discovery_chr8_IDO1.txt.gz", "ukb_ppp_raw_ido1_assay_a.txt")])
            for source in (standalone, archive):
                with self.subTest(source=source.name):
                    output = root / source.stem
                    self.prepare(source, output, extra=("--limit-type", "bytes", "--limit", str(limit)))
                    rows = self.rows(output / "std/EUR/batch_001/IDO1.tsv")
                    self.assertEqual(["rs_assay_a"], [row["rsid"] for row in rows])

    def test_missing_or_ambiguous_chromosome_members_fail_before_outputs(self):
        cases = [
            ("missing", [("discovery_chr17_ACE.txt.gz", "ukb_ppp_raw_ace.txt")],
             "exactly one summary member for chromosome 8"),
            ("ambiguous", [("discovery_chr8_assay_a.txt.gz", "ukb_ppp_raw_ido1_assay_a.txt"),
                           ("discovery_chr8_assay_b.txt.gz", "ukb_ppp_raw_ido1_assay_b.txt")],
             "multiple members for chromosome 8"),
            ("mislabeled", [("discovery_chr8_IDO1.txt.gz", "ukb_ppp_raw_ace.txt")],
             "gene chromosome"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for label, members, error in cases:
                with self.subTest(case=label):
                    source = root / (label + ".tar")
                    self.archive(source, members)
                    output = root / label
                    completed = self.prepare(source, output, check=False)
                    self.assertNotEqual(0, completed.returncode)
                    self.assertIn(error, completed.stderr)
                    self.assertFalse((output / "std").exists())
                    self.assertFalse((output / "inst").exists())

    def test_batch_run_preserves_each_assay_ancestry_and_real_batch_number(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = [sys.executable, str(ROOT / "scripts/ukb_ppp_batch_manifest_runner_fast.py"),
                       "--base", str(root / "raw"), "--qc-dir", str(root / "qc"),
                       "--outdir", str(root / "out"), "--standardized-dir", str(root / "std"),
                       "--instrument-dir", str(root / "inst"), "--download-manifest",
                       str(self.FIXTURES / "ukb_ppp_preparation_manifest.tsv"),
                       "--gene-coordinate-file", str(self.FIXTURES / "gene_coordinates_hg38.tsv"),
                       "--batch-size", "1", "--focus-gene", "IDO1", "--run", "--stop-on-error"]
            subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
            progress = self.rows(root / "qc/batch_progress.tsv")
            self.assertEqual(4, len(progress))
            self.assertEqual({"complete"}, {row["status"] for row in progress})
            self.assertEqual({"batch_002"}, {row["batch_id"] for row in progress})
            for field in ("standardized_path", "instrument_path", "legacy_path", "qc_path"):
                self.assertEqual(4, len({row[field] for row in progress}))
                self.assertTrue(all(Path(row[field]).is_file() for row in progress))
            for task in progress:
                qc = self.rows(Path(task["qc_path"]))[0]
                self.assertEqual(tuple(task[key] for key in ("gene", "ancestry", "batch_id", "source_id")),
                                 tuple(qc[key] for key in ("gene", "ancestry", "batch_id", "source_id")))
                expected = "rs_assay_a" if "assay_a" in task["source_file"] else "rs_assay_b"
                for field in ("instrument_path", "legacy_path"):
                    rows = self.rows(Path(task[field]))
                    self.assertEqual([expected], [row["rsid"] for row in rows])
                    self.assertEqual(task["ancestry"], rows[0]["ancestry"])
                self.assertEqual(4, len(list((root / "inst").glob("*/batch_002/IDO1/*.tsv"))))
            self.assertFalse(any((root / "std").glob("*/batch_001")))


if __name__ == "__main__": unittest.main()
