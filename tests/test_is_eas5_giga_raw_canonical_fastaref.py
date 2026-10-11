"""Original GIGASTROKE EAS AIS raw<->canonical FASTA source invariants."""
import importlib.util,sys,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas5_giga_raw_canonical_fastaref.py"
sp=importlib.util.spec_from_file_location("giga_raw_canonical_test",p)
m=importlib.util.module_from_spec(sp);sys.modules[sp.name]=m
sp.loader.exec_module(m)
sys.path.insert(0,str(p.parent))
from gigastroke_fasta_variant_contract import normalize_record
class Fasta:
 def get(self,c,p):
  if (c,p)!=("5",29439275):raise ValueError("Unexpected coordinate")
  return "T"
class Tests(unittest.TestCase):
 def test_numeric_input_equivalence_tight(self):
  self.assertTrue(m.eq_float("0.0106","0.0106"))
  self.assertFalse(m.eq_float("0.0106","-0.0106"))
  self.assertFalse(m.eq_float("nan","nan"))
 def test_raw_canonical_fasta_recoded_variant(self):
  key=("5",29439275,"T","C")
  row={"chrom":"5","pos":29439275,"effect_allele":"T","other_allele":"C",
    "beta":"0.0106","se":"0.0136","p":"0.4351","eaf":"0.7798"}
  can={**row,"build":"GRCh37","dataset":"GCST90104545",
    "ancestry":"EAS","historical_variant_id":"5:29439275:C:T"}
  rows,stat=m.audit({key:row},{key:can},{"5:29439275:T:C"},Fasta(),normalize_record)
  self.assertEqual(len(rows),1)
  self.assertEqual(rows[0]["historical_id_order"],"OTHER_EFFECT_ORDER")
  self.assertAlmostEqual(float(rows[0]["beta_effect_ALT"]),-.0106)
  self.assertEqual(stat["RAW_CANONICAL_ALLELE_EFFECT_FREQ_NUMERIC_AND_FASTA_VERIFIED"],1)
 def test_raw_canonical_numeric_divergence_fails_closed(self):
  key=("5",29439275,"T","C")
  row={"chrom":"5","pos":29439275,"effect_allele":"T","other_allele":"C",
    "beta":"0.0106","se":"0.0136","p":"0.4351","eaf":"0.7798"}
  can={**row,"build":"GRCh37","dataset":"GCST90104545",
    "ancestry":"EAS","historical_variant_id":"5:29439275:C:T",
    "beta":"0.5"}
  with self.assertRaises(ValueError):
   m.audit({key:row},{key:can},{"5:29439275:T:C"},Fasta(),normalize_record)
 def test_raw_source_allele_pair_mismatch_not_repaired(self):
  key=("5",29439275,"A","C")
  row={"chrom":"5","pos":29439275,"effect_allele":"A","other_allele":"C",
    "beta":"0.0106","se":"0.0136","p":"0.4351","eaf":"0.7798"}
  can={**row,"build":"GRCh37","dataset":"GCST90104545",
    "ancestry":"EAS","historical_variant_id":"5:29439275:C:A"}
  with self.assertRaises(ValueError):
   m.audit({key:row},{key:can},{"5:29439275:T:C"},Fasta(),normalize_record)
if __name__=="__main__":unittest.main()
