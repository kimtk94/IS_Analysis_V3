"""Offline source-format regression tests for Japanese OASIS sc-eQTL nominal."""
import gzip
import io
import json
import sys
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_egead1054_nominal_archive import (
    REQUIRED_COLUMNS, RANGE_BYTES, build, safe_extract_prefix
)


def synthetic_zip(first_row_p=0.3488858046, wrong_schema=False, all_significant=False):
    lines = ["\t".join(REQUIRED_COLUMNS if not wrong_schema else ("gene", "phenotype_id"))]
    for i in range(120):
        p = 0.001 if all_significant else first_row_p if i == 0 else 0.6
        lines.append("\t".join(["ENSG00000237491", "LINC01409",
            f"chr1_{100000+i}_A_G", str(i-100), str(p), "0.1817",
            "0.1934", "0.0739", "29", "29"]))
    inner = gzip.compress(("\n".join(lines)+"\n").encode("utf8"))
    tarbytes = io.BytesIO()
    with tarfile.open(fileobj=tarbytes, mode="w") as tf:
        ti = tarfile.TarInfo("B_Activated_PC15_MAF0.05.cis_nominal.txt.gz")
        ti.size = len(inner)
        tf.addfile(ti, io.BytesIO(inner))
    zipbytes = io.BytesIO()
    with zipfile.ZipFile(zipbytes, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("eQTL_summary_statistics.tar", tarbytes.getvalue())
    buf = zipbytes.getvalue()
    if len(buf) > RANGE_BYTES:
        raise ValueError("fixture unexpectedly bigger than range")
    return buf + b"\0" * (RANGE_BYTES - len(buf))


class TestEGEAD1054NominalProbe(unittest.TestCase):
    def test_partial_zip_tar_gzip_nominal_with_non_significant(self):
        result = safe_extract_prefix(synthetic_zip())
        self.assertEqual(result["zip_member_name"], "eQTL_summary_statistics.tar")
        self.assertTrue(result["tar_first_member_name"].endswith("cis_nominal.txt.gz"))
        self.assertEqual(result["nominal_qtl_columns"], list(REQUIRED_COLUMNS))
        self.assertGreater(result["sample_p_greater_than_0_05"], 0)
        self.assertFalse(result["sample_contains_g0022_verified"])

    def test_filtered_only_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "non-significant"):
            safe_extract_prefix(synthetic_zip(all_significant=True))

    def test_invalid_schema_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "schema"):
            safe_extract_prefix(synthetic_zip(wrong_schema=True))

    def test_non_zip_header_fails_closed(self):
        corrupt = b"BAD!" + synthetic_zip()[4:]
        with self.assertRaisesRegex(ValueError, "ZIP"):
            safe_extract_prefix(corrupt)

    def test_offline_output_cannot_claim_ready(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            blob=root/"head.bin"
            blob.write_bytes(synthetic_zip())
            class Args:
                outdir=root/"output"
                local_prefix=blob
                probe_remote=False
            response=build(Args)
            self.assertFalse(response["full_cis_coloc_ready"])
            self.assertFalse(response["http_range_enforced"])
            self.assertIsNone(response["eas_AIS_four_cs_variant_coverage_verified"])
            self.assertEqual(response["eas_AIS_four_cs_variant_coverage_status"], "NOT_ASSESSED")
            self.assertTrue((Args.outdir/"G0022_EGEAD1054_SOURCE_NOMINAL_PROBE.json").is_file())


if __name__ == "__main__":
    unittest.main()
