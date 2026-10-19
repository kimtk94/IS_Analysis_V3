"""Non-destructive, fixture-level tests for independent R3_7 audit report."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"scripts/is/audit_phase9f_r3_7_donor_results.py"
spec = importlib.util.spec_from_file_location("audit_phase9f_donor",SOURCE)
auditor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(auditor)

class TestIndependentDonorAudit(unittest.TestCase):
    def test_nan_does_not_become_zero(self):
        import math
        self.assertTrue(math.isnan(auditor.num("NA")))
        self.assertTrue(math.isnan(auditor.num("")))
        self.assertEqual(auditor.num("0"),0)
        with self.assertRaises(ValueError):
            auditor.num("inf")

    def test_fails_closed_without_all_source_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/"no_data"
            source.mkdir()
            output=Path(tmp)/"output"
            with self.assertRaises(FileNotFoundError):
                auditor.run(source,output)
            self.assertFalse(output.exists())

    def test_fails_closed_if_output_nonempty(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/"source"
            source.mkdir()
            output=Path(tmp)/"output"
            output.mkdir()
            (output/"keep.txt").write_text("HISTORICAL")
            with self.assertRaises(FileExistsError):
                auditor.run(source,output)
            self.assertEqual((output/"keep.txt").read_text(),"HISTORICAL")

    def test_svg_figures_are_well_formed(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp=Path(tmp)/"comparisons.svg"
            dp=Path(tmp)/"donors.svg"
            rows=[dict(gene="COL4A2",celltype_a="Pericytes",
                       celltype_b="Endothelial cells",n=4,
                       n_positive=4,detect_positive=2,
                       median=1.5,individual_deltas=[1.0,1.2,1.8,2.0])]
            auditor.render_chart(rows,cp)
            auditor.render_balance({
                ("D1","Pericytes"):25,("D1","Endothelial cells"):30,
                ("D2","Pericytes"):40,("D2","Endothelial cells"):20
            },dp)
            for path in (cp,dp):
                ET.parse(path)
                self.assertGreater(path.stat().st_size,800)

if __name__=="__main__":
    unittest.main()
