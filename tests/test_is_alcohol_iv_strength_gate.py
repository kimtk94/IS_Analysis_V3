"""Non-ALDH2 positional sentinel F proxies are NOT conditional independent IVs."""
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_alcohol_iv_strength_gate import audit
class Strength(unittest.TestCase):
 def source(self,path):
  rows=[]
  for i in range(6):
   status="NOT_IN_EAS_AIS_SOURCE" if i==5 else "ALLELE_HARMONIZED"
   rows.append({"variant_grch37":f"{i+1}:100000:A:G",
      "alcohol_beta_ALT":".2","alcohol_se":".02","alcohol_p":"1e-10",
      "alcohol_ALT_eaf":".20","ais_ALT_eaf":".203" if i<5 else "",
      "ais_p":".21" if i<5 else "",
      "ais_qc":status,"independent_genetic_IV_certified":"False"})
  with path.open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
   w.writeheader();w.writerows(rows)
 def test_unconditional_F_proxy_and_overlap_gates(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);p=root/"s.tsv";self.source(p)
   x=audit(p,root/"out")
   self.assertEqual(x["candidate_positional_loci"],6)
   self.assertEqual(x["AIS_exact_allele_EAF_within_10pp"],5)
   self.assertEqual(x["AIS_SNP_missing"],1)
   self.assertAlmostEqual(x["F_proxy_min_across_all"],100.)
   self.assertFalse(x["LD_independence_validated"])
   self.assertFalse(x["MR_computed"])
 def test_false_source_IV_independence_is_refused(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);p=root/"s.tsv";self.source(p)
   s=p.read_text().replace("False","True",1);p.write_text(s)
   with self.assertRaisesRegex(ValueError,"falsely declared"):
    audit(p,root/"out")
if __name__=="__main__":unittest.main()
