#!/usr/bin/env python3
"""Strict synthetic signed-Z transformation tests; no simulated stroke evidence."""
import io
import math
from pathlib import Path
import sys
import tarfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_ImmuNexUT_eCAVIAR_signed_z import (
    signed_z_from_slope_p,first_gene_source_rows,convert_block
)
from audit_EGEAD420_nominal_zip64 import IMMUNEXUT_13_COL_HEADER


def fixture(n=7,complete=True):
    f=io.BytesIO()
    lines=["\t".join(IMMUNEXUT_13_COL_HEADER)]
    for i in range(n):
        lines.append("\t".join([
            "ENSG00000227232.5","WASH7P","chr1","29570",str(n),
            f"rs{i+1}", "TTTTC" if i==0 else "A",
            "T" if i==0 else "G",
            "chr1",str(54712+i),str(54716+i) if i==0 else str(54712+i),
            "0.05" if i==0 else "0.7",
            "-0.2" if i==0 else "0.3",
        ]))
    if complete:
        lines.append("\t".join([
            "ENSG0000011","TEST","chr1","80000","2","rs99",
            "G","C","chr1","80001","80001","0.1","0.15"]))
    body=("\n".join(lines)+"\n").encode()
    with tarfile.open(fileobj=f,mode="w:gz") as ar:
        i=tarfile.TarInfo("CD16p_Mono_nominal.txt")
        i.size=len(body)
        ar.addfile(i,io.BytesIO(body))
    return f.getvalue()


class SignedZInputTests(unittest.TestCase):
    def test_normal_2_sided_p_to_signed_z(self):
        self.assertAlmostEqual(signed_z_from_slope_p(0.2,0.05),1.9599639845,places=6)
        self.assertAlmostEqual(signed_z_from_slope_p(-0.2,0.05),-1.9599639845,places=6)
        self.assertAlmostEqual(signed_z_from_slope_p(0,1.0),0.0,places=10)

    def test_extreme_p_survives_float_without_rounding_one_minus_p(self):
        z=signed_z_from_slope_p(-0.1,1e-90)
        self.assertLess(z,-19)
        self.assertTrue(math.isfinite(z))

    def test_zero_or_negative_p_rejected(self):
        with self.assertRaisesRegex(ValueError,"finite signed Z"):
            signed_z_from_slope_p(0.2,0)
        with self.assertRaisesRegex(ValueError,"finite signed Z"):
            signed_z_from_slope_p(0.2,-0.1)
        with self.assertRaisesRegex(ValueError,"zero beta"):
            signed_z_from_slope_p(0,0.01)

    def test_first_gene_real_style_13col_source_cis_block(self):
        rows,meta=first_gene_source_rows(fixture(7))
        self.assertEqual(len(rows),7)
        self.assertEqual(meta["source_declared_cis_variants"],7)
        self.assertTrue(meta["source_next_gene_boundary_observed"])
        converted,status=convert_block(rows)
        self.assertEqual(status["real_original_gene_cis_rows"],7)
        self.assertEqual(status["indel_records"],1)
        self.assertEqual(status["positive_signed_z"],6)
        self.assertEqual(status["negative_signed_z"],1)
        self.assertFalse(status["GWAS_harmonized"])
        self.assertFalse(status["eCAVIAR_CLPP_calculated"])
        self.assertIn("TTTTC>T",converted[0]["original_grch38_variant_alleles"])

    def test_partial_original_first_gene_is_hold(self):
        with self.assertRaisesRegex(ValueError,"incomplete"):
            first_gene_source_rows(fixture(7,False))

    def test_only_published_alt_effect_direction_used(self):
        rows,_=first_gene_source_rows(fixture())
        rows[0][7]="TTTTC"
        with self.assertRaisesRegex(ValueError,"alleles"):
            convert_block(rows)


if __name__=="__main__":
    unittest.main()
