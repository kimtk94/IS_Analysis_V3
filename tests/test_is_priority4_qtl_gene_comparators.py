import csv,tempfile,unittest
from pathlib import Path
from importlib import import_module
_m = import_module("scripts.is.audit_is_priority4_qtl_gene_comparators")
GENES,CUTOFF = _m.GENES,_m.CUTOFF
class TestPriority4ComparatorContracts(unittest.TestCase):
 def test_expected_genes_and_cutoff(self):
  self.assertEqual(set(GENES),{"CALHM2","NEURL","FGF5","C4orf22"})
  self.assertEqual(CUTOFF,1e-4)
 def test_script_exposes_fail_closed_and_provenance(self):
  path=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_priority4_qtl_gene_comparators.py"
  source=path.read_text()
  for phrase in ("source_sha256","Input count drift","source_rows","not_colocalization","original_2225_IS_genes_unchanged"):
   if phrase=="Input count drift":
    self.assertIn("count drift",source)
   else:self.assertIn(phrase,source)
