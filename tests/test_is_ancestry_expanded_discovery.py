"""Ancestry expansion keeps original candidate scope, separates EUR provisional loci."""
import csv,gzip,tempfile,unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from build_ancestry_expanded_discovery import build
class TestEURExpansion(unittest.TestCase):
 def test_eur_gws_suggestive_and_zero_p_underflow(self):
  with tempfile.TemporaryDirectory() as td:
   base=Path(td)/"base";base.mkdir();out=Path(td)/"out"
   e=Path(td)/"eur.gz"
   src=[]
   for i in range(30):
    src.append({"group_id":f"G{i}","chr":str(i%22+1),
      "region_start":"100","region_end":"200",
      "window_start":"1","window_end":"500","lead_variant":"1:100:A:T",
      "lead_p":"1e-7","phenotypes":"AIS","has_gws":"0"})
   with (base/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open("w") as h:
    w=csv.DictWriter(h,fieldnames=list(src[0]),delimiter="\t");w.writeheader();w.writerows(src)
   r=[{"chr":"1","pos":"150","p":"0","eaf":"0.3","variant_pair_id":"1:150:A:G"},
      {"chr":"1","pos":"300","p":"1e-7","eaf":"0.4","variant_pair_id":"1:300:C:T"},
      {"chr":"2","pos":"2000000","p":"8e-7","eaf":"0.2","variant_pair_id":"2:2000000:C:G"},
      {"chr":"3","pos":"400","p":"1e-9","eaf":"0.0001","variant_pair_id":"3:400:A:C"}]
   with gzip.open(e,"wt") as h:
    w=csv.DictWriter(h,fieldnames=list(r[0]),delimiter="\t");w.writeheader();w.writerows(r)
   summary=build(base,e,out,1000000)
   self.assertEqual(summary["original_eas_japanese_components"],30)
   self.assertEqual(summary["eur_distance_clusters"],2)
   self.assertEqual(summary["eur_gws_clusters"],1)
   self.assertEqual(summary["eur_suggestive_only_clusters"],1)
   self.assertEqual(summary["combined_rows"],32)
   self.assertEqual(summary["input_variant_qc"]["maf_lt_001"],1)
   with (out/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as h:
    rr=list(csv.DictReader(h,delimiter="\t"))
   self.assertEqual(len(rr),32)
   self.assertEqual(sum(r["source"]=="BBJ_GIGASTROKE_EAS" for r in rr),30)
   self.assertTrue(all(r["provisional_status"].endswith("NOT_LD_INDEPENDENT") for r in rr))
if __name__=="__main__": unittest.main()
