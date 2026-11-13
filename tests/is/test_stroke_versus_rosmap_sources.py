#!/usr/bin/env python3
"""Synthetic GEO SOFT tests; never infer eQTL from transcriptome only."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_stroke_versus_rosmap_sources import (  # noqa:E402
    SERIES, geo_soft_url, parse_soft, derive_group
)


class StrokeGeoMetadataTests(unittest.TestCase):
    def test_known_stroke_GEO_list(self):
        self.assertEqual(set(SERIES), {"GSE331114", "GSE267023", "GSE285659"})
        self.assertFalse(any(d["genotype_mapped_eqtl"] for d in SERIES.values()))

    def test_official_ncbi_ftp(self):
        self.assertEqual(geo_soft_url("GSE331114"),
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE331nnn/GSE331114/soft/GSE331114_family.soft.gz")
        with self.assertRaises(ValueError):
            geo_soft_url("GSE331114/../../")

    def test_parse_stroke_brain_sample_groups(self):
        content = """^SERIES = GSE331114
!Series_title = Human acute IS brain
^SAMPLE = GSM111
!Sample_title = AIS_brain1
!Sample_source_name_ch1 = brain
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = group: AIS
!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/GSM111_matrix.mtx.gz
^SAMPLE = GSM222
!Sample_title = Control_brain1
!Sample_source_name_ch1 = brain
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = group: Control
!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/GSM222_matrix.mtx.gz
"""
        info, samples = parse_soft(content)
        self.assertEqual(info["accession"], "GSE331114")
        self.assertEqual(len(samples), 2)
        self.assertEqual(derive_group(samples["GSM111"]), "AIS")
        self.assertEqual(derive_group(samples["GSM222"]), "Control")

    def test_duplicate_sample_rejected(self):
        with self.assertRaises(ValueError):
            parse_soft("^SAMPLE = GSM1\n^SAMPLE = GSM1\n")

    def test_unlabeled_sample_is_unknown_not_healthy(self):
        self.assertEqual(derive_group({"title": "baseline", "characteristics": []}), "UNKNOWN")


if __name__ == "__main__":
    unittest.main()
