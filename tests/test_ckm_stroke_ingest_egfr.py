import gzip
import io
import tempfile
import unittest
from pathlib import Path
from scripts.ckm_stroke import ingest_egfr, normalize_sumstats


class Response:
    status = 200
    def __init__(self, blob, url, size=None):
        self.blob, self.url = io.BytesIO(blob), url
        self.headers = {"Content-Length": str(len(blob) if size is None else size)}
    def __enter__(self): return self
    def __exit__(self, *args): return False
    def read(self, n): return self.blob.read(n)
    def geturl(self): return self.url


def mk_source():
    return gzip.compress(
        b"Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum\n"
        b"1 736289 rs79010578 a t 0.14 0.000315 0.0007081 0.6563 347599\n")


class TestEgfrSource(unittest.TestCase):
    def test_eur_full_read_is_not_official_md5(self):
        blob = mk_source()
        with tempfile.TemporaryDirectory() as temp:
            x = ingest_egfr.ingest(temp, max_mb=2, opener=lambda req, timeout:
                                   Response(blob, ingest_egfr.EGFR["url"]), retry_sleep=lambda _: None)
            self.assertEqual(x["status"], "PASS")
            self.assertEqual(x["scientific_status"], "HOLD")
            self.assertFalse(x["local_sha256_verified_against_official"])
            self.assertEqual(x["source_checksum_status"], "UNVERIFIED_OFFICIAL_CHECKSUM")
            self.assertEqual(x["input_rows"], 1)
            self.assertEqual(x["normalized_rows"], 1)

    def test_invalid_length_is_fail(self):
        blob = mk_source()
        with tempfile.TemporaryDirectory() as temp:
            x = ingest_egfr.ingest(temp, max_mb=2, opener=lambda req, timeout:
                                   Response(blob, ingest_egfr.EGFR["url"], len(blob)+5),
                                   retry_sleep=lambda _: None)
            self.assertEqual(x["status"], "FAIL")
            self.assertIn("length", x["reason"])

    def test_truncated_gzip_fails(self):
        blob = mk_source()[:-5]
        with tempfile.TemporaryDirectory() as temp:
            x = ingest_egfr.ingest(temp, max_mb=2, opener=lambda req, timeout:
                                   Response(blob, ingest_egfr.EGFR["url"]), retry_sleep=lambda _: None)
            self.assertEqual(x["status"], "FAIL")

    def test_exact_missing_value_classification(self):
        columns=normalize_sumstats.MAP["ckdgen2019"]
        row={"Chr":"1", "Pos_b37":"123", "Allele1":"A", "Allele2":"G",
             "Effect":"NA","StdErr":"0.1","P-value":"0.2"}
        reason, parsed=normalize_sumstats.classify(row, columns)
        self.assertEqual(reason, "missing_value")
        self.assertIsNone(parsed)
        row["Effect"]="nan"
        self.assertEqual(normalize_sumstats.classify(row, columns)[0], "non_finite")
        row["Effect"]="0.1"
        row["P-value"]="0"
        self.assertEqual(normalize_sumstats.classify(row, columns)[0], "p_underflow")

    def test_missing_and_malformed_counts(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            raw=root/"source.gz"
            with gzip.open(raw, "wt") as fh:
                fh.write("Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum\n")
                fh.write("1 123 rs1 A G 0.2 0.1 0.03 0.01 100\n")
                fh.write("1 124 rs2 A G 0.2 NA 0.03 0.01 100\n")
                fh.write("1 malformed\n")
            x=normalize_sumstats.normalize(raw, root/"out.gz", "ckdgen2019")
            self.assertEqual(x["status"], "HOLD")
            self.assertEqual(x["input_rows"], 3)
            self.assertEqual(x["excluded_total"], 2)
            self.assertEqual(x["excluded"]["missing_value"], 1)
            self.assertEqual(x["excluded"]["malformed_row"], 1)


    def test_recover_blank_rsid_losslessly(self):
        fields = ["Chr", "Pos_b37", "RSID", "Allele1", "Allele2",
                  "Freq1", "Effect", "StdErr", "P-value", "n_total_sum"]
        original_tokens = ["1", "2556125", "t", "c", "0.3260",
                           "-0.0035140", "0.0021112", "0.09603", "547349"]
        recovered = normalize_sumstats.recover_missing_rsid(fields, original_tokens)
        self.assertIsNotNone(recovered)
        self.assertEqual(recovered["RSID"], "")
        self.assertEqual(recovered["Allele1"], "t")
        self.assertEqual(recovered["Effect"], "-0.0035140")
        self.assertIsNone(normalize_sumstats.classify(recovered, normalize_sumstats.MAP["ckdgen2019"])[0])

    def test_do_not_recover_other_missing_columns(self):
        fields = ["Chr", "Pos_b37", "RSID", "Allele1", "Allele2",
                  "Freq1", "Effect", "StdErr", "P-value", "n_total_sum"]
        # Here the third token is an RSID, so some other field was lost.
        wrong = ["1", "2556125", "rs123", "t", "0.3260",
                 "-0.0035140", "0.0021112", "0.09603", "547349"]
        self.assertIsNone(normalize_sumstats.recover_missing_rsid(fields, wrong))

    def test_normalize_recovered_rsid_is_not_exclusion(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "src.gz"
            with gzip.open(raw, "wt") as fh:
                fh.write("Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum\n")
                fh.write("1 200 rs2 A G 0.2 0.1 0.03 0.01 100\n")
                fh.write("1 2556125 t c 0.3260 -0.0035140 0.0021112 0.09603 547349\n")
            result = normalize_sumstats.normalize(raw, Path(d) / "normalized.gz", "ckdgen2019")
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["output_rows"], 2)
            self.assertEqual(result["excluded_total"], 0)
            self.assertEqual(result["recovered_missing_rsid"], 1)


if __name__ == "__main__":
    unittest.main()
