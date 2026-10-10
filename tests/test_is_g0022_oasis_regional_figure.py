"""Generated descriptive figure must not imply coloc or causal mediation."""
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from render_g0022_oasis_regional_overlap_svg import format_svg


class TestOasisFigure(unittest.TestCase):
    def test_empty_mono_source_is_blocked(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError,"monocyte"):
                format_svg(Path(t)/"figure.svg",[])

    def test_svg_is_valid_and_contains_explicit_caveats(self):
        rows=[]
        for gene in ("ALDH2","BRAP","RPH3A"):
            rows.append({"gene":gene,"cell":"Mono-L1",
                         "mapped_grch38_pos":111730205,
                         "gwas_p":9e-18,"oasis_p":0.3})
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/"output.svg"
            result=format_svg(p,rows)
            ET.parse(p)
            svg=p.read_text()
            self.assertEqual(result["gene_mono_matched_counts"]["ALDH2"],1)
            self.assertIn("not evidence of shared causal variant",svg)
            self.assertIn("rs671 (NOT observed",svg)
            self.assertIn("RPH3A",svg)


if __name__=="__main__":
    unittest.main()
