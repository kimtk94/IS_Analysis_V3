"""PLINK2 clump screen only selects study-group GWS and never claims causality."""
import csv,gzip,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from run_gws_ld_clump_screen import prepare
class ClumpScreenTest(unittest.TestCase):
 def test_gws_only_and_no_cross_study_pooling(self):
  with tempfile.TemporaryDirectory() as temp:
   root=Path(temp)
   readiness=[{"dataset":f"S{i//7}","group_id":f"G{i%7}"} for i in range(42)]
   with (root/"IS_GWS_STUDY_INPUT_READINESS.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(readiness[0]),delimiter="\t")
    w.writeheader();w.writerows(readiness)
   records=[{"dataset":"S0","group_id":"G0","qc_status":"MATCH_ALT_EFFECT",
             "p":"1e-9","variant_id":"2:100:A:G"},
            {"dataset":"S0","group_id":"G0","qc_status":"MATCH_ALT_EFFECT",
             "p":"1e-5","variant_id":"2:101:A:G"},
            {"dataset":"S0","group_id":"G1","qc_status":"MATCH_ALT_EFFECT",
             "p":"1e-6","variant_id":"2:200:C:G"},
            {"dataset":"S1","group_id":"G0","qc_status":"MATCH_ALT_EFFECT",
             "p":"1e-4","variant_id":"2:100:A:G"}]
   with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","wt") as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
    w.writeheader();w.writerows(records)
   out=prepare(root)
   self.assertEqual(len(out),1)
   self.assertEqual(out[0]["dataset"],"S0")
   self.assertEqual(out[0]["group_id"],"G0")
   self.assertEqual(out[0]["tested_variant_count"],2)
   self.assertEqual(out[0]["status"],"BLOCKED_NO_REFERENCE")
if __name__=="__main__":
 unittest.main()
