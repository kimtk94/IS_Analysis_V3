"""Coverage comparison must be exactly same GWAS windows and never lose matches."""
import csv
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from compare_expanded_gws_ld_coverage import compare

def save(path,rows):
    with path.open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def entry(d,g,n,missing):
    return {"dataset":d,"group_id":g,"window_variants":"200",
            "matched_oriented_nonpal":str(n),"no_reference_position":str(missing)}
class CoverageTest(unittest.TestCase):
 def test_new_coverage_increase_and_other_stable(self):
  with tempfile.TemporaryDirectory() as td:
   r=Path(td);v=r/"expanded_reference_v2";v.mkdir()
   old=[];new=[]
   for i in range(42):
    g="IS_XDATA_G0015" if i<6 else "IS_XDATA_G0022" if i<12 else "IS_XDATA_G0004"
    ds=f"D{i}"
    old.append(entry(ds,g,100,100))
    new.append(entry(ds,g,140,60) if i<12 else entry(ds,g,100,100))
   save(r/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv",old)
   save(v/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv",new)
   summary=compare(r)
   self.assertEqual(summary["newly_oriented_total"],480)
   self.assertEqual(summary["new_missing_total"],42*100-12*40)
   new[4]["matched_oriented_nonpal"]="99"
   save(v/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv",new)
   with self.assertRaisesRegex(ValueError,"lost"):
    compare(r)
if __name__=="__main__":unittest.main()
