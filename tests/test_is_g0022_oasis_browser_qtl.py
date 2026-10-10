import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_oasis_browser_qtl import parse

class TestOasisBrowser(unittest.TestCase):
 def test_missing_bokeh_must_fail_not_infer_negative(self):
  with self.assertRaisesRegex(ValueError,"Missing embedded"):
   parse("<html><body>No source</body></html>","ALDH2")
 def test_wrong_gene_or_column_validation(self):
  # Fails closed without Bokeh JSON, not equivalent to zero true QTL signal.
  with self.assertRaises(ValueError):
   parse('<script type="application/json">{"x":{}}</script>',"RPH3A")

if __name__=="__main__":unittest.main()
