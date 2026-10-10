"""Guard 30 groups / allele-aware six-outcome GWAS coverage as noncausal audit."""
import csv,gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_30_eas_regions_six_gwas_sources.py"
spec=importlib.util.spec_from_file_location("is30gwas",P)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)
def fixture_regions():
 return [{
  "group_id":f"IS_XDATA_G{i:04d}","chr":"1","start":str(i*1000000),
  "end":str(i*1000000),"lead_variant":f"1:{i*1000000}:A:G",
  "phenotypes":"AIS","has_gws":"0",
  "window_start_GRCh37":i*1000000-500000,
  "window_end_GRCh37":i*1000000+500000
 } for i in range(1,31)]
class TestSourceScan(unittest.TestCase):
 def test_no_false_replication_or_allele_flip(self):
  groups=fixture_regions();bychr={"1":groups}
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"test.tsv.gz"
   h=["dataset","phenotype","ancestry","build","chr","pos","effect_allele",
      "other_allele","beta","se","p","eaf","or","variant_id"]
   data=[
      ["GWS","AIS","EAS","GRCh37","1","1000000","A","G","0.2",".03","1e-9",".7","1","1:1000000:A:G"],
      ["GWS","AIS","EAS","GRCh37","1","1000001","X","G","0.2",".03","1e-8",".7","1","1:1000001:A:G"],
      ["GWS","AIS","EAS","GRCh37","1","700000","G","A","0.1",".1",".04",".2","1","1:700000:A:G"]
   ]
   with gzip.open(p,"wt") as f:
    f.write("\t".join(h)+"\n")
    for row in data:f.write("\t".join(row)+"\n")
   records,meta=m.scan(p,"TEST",groups,bychr)
   self.assertEqual(len(records),30)
   first=records[0]
   self.assertEqual(first["n_valid_variant_allele_records"],2)
   self.assertEqual(first["n_invalid_or_ambiguous_record"],1)
   self.assertTrue(first["group_source_lead_present"])
   self.assertEqual(first["n_genomewide_p_lt_5e8"],1)
   self.assertAlmostEqual(float(first["group_source_lead_beta_ALT"]),-.2)
   self.assertAlmostEqual(float(first["group_source_lead_ALT_AF"]),.3)
   self.assertEqual(meta["total_native_rows_scanned"],3)
 def test_allele_set_invariant(self):
  self.assertTrue(m.variant_same("1:123:G:A","1:123:A:G"))
  self.assertFalse(m.variant_same("1:123:G:A","1:124:A:G"))
  self.assertFalse(m.variant_same("1:123:G:A","1:123:A:C"))
 def test_wrong_region_count_rejected(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"regions.tsv"
   with p.open("w") as f:
    f.write("group_id\tchr\tstart\tend\nIS_XDATA_G0001\t1\t100\t100\n")
   with self.assertRaises(ValueError):m.regions_read(p,500000)
if __name__=="__main__":unittest.main()
