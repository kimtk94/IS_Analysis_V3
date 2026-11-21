"""Bounded first-party stroke log(OR) scale fixtures and missing N detection."""
from __future__ import annotations
import gzip
import math
import unittest
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import COLUMNS
from scripts.ckm_stroke.audit_gigastroke_logor_prefix import (
    READ_BYTES, validate_prefix
)


def sample(*,beta="0.1", odds=None, schema=COLUMNS, rows=600):
    odds=odds or str(math.exp(float(beta)))
    values=["1","12345","0.25",beta,"0.05","0.04",
            odds,"1.01","1.3","A","C"]
    text="\t".join(schema)+"\n"+("\t".join(values)+"\n")*rows
    data=gzip.compress(text.encode("ascii"),mtime=0)
    if len(data)>READ_BYTES:
        raise ValueError("fixture gzip too large")
    return data.ljust(READ_BYTES,b"\x00")


class TestLogOR(unittest.TestCase):
    def test_correct_original_effect_scale_is_bounded_and_not_h2(self):
        report=validate_prefix("AIS","GCST90104540",sample())
        self.assertEqual(report["observed"]["parsed_effect_rows"],600)
        self.assertLess(report["maximum_abs_beta_minus_lnOR"],1e-10)
        self.assertFalse(report["source_per_variant_N_available"])
        self.assertTrue(report["bounded_sample_not_genomewide_proof"])

    def test_original_log_or_disagreement_fails(self):
        with self.assertRaisesRegex(ValueError,"beta != log"):
            validate_prefix("AIS","GCST90104540",sample(beta="0.7",odds=str(math.exp(.1))))

    def test_reversed_CE_LAA_fails(self):
        with self.assertRaisesRegex(ValueError,"accession"):
            validate_prefix("CE","GCST90104542",sample())

    def test_changed_header_fails(self):
        header=list(COLUMNS)
        header[-1]="N"
        with self.assertRaisesRegex(ValueError,"schema changed"):
            validate_prefix("AIS","GCST90104540",sample(schema=tuple(header)))

    def test_truncated_prefix_fails(self):
        with self.assertRaisesRegex(ValueError,"gzip prefix bytes"):
            validate_prefix("AIS","GCST90104540",gzip.compress(b"incomplete"))


if __name__=="__main__":
    unittest.main()
