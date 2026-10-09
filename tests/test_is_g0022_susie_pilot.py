"""G0022 pilot signed-LD dosage orientation and truncated-window safety tests."""
import csv,gzip,sys,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_g0022_susie_pilot import dosage_alt_orientation,source_variants,traw_rows

class G0022LDPilotTest(unittest.TestCase):
 def test_alt_dose_signed_flip(self):
  v=np.array([0.,1.,2.,np.nan])
  np.testing.assert_allclose(dosage_alt_orientation(v,"A","G","A","G"),np.array([2.,1.,0.,np.nan]))
  np.testing.assert_allclose(dosage_alt_orientation(v,"G","A","A","G"),v)
  with self.assertRaisesRegex(ValueError,"incompatible"):
   dosage_alt_orientation(v,"C","A","A","G")
 def test_skip_palindromic_and_keep_allele_oriented_gwas(self):
  with tempfile.TemporaryDirectory() as temp:
   p=Path(temp)/"gwas.gz"
   fields=["group_id","dataset","qc_status","ref","alt","variant_id",
           "alt_effect_beta","se","p","alt_effect_eaf"]
   rows=[dict(group_id="IS_XDATA_G0022",dataset="GCST90104545",
       qc_status="MATCH_ALT_EFFECT",ref="A",alt="G",variant_id="12:100:A:G",
       alt_effect_beta="0.2",se="0.1",p="0.04",alt_effect_eaf="0.35"),
    dict(group_id="IS_XDATA_G0022",dataset="GCST90104545",
       qc_status="MATCH_ALT_EFFECT",ref="A",alt="T",variant_id="12:110:A:T",
       alt_effect_beta="0.2",se="0.1",p="0.04",alt_effect_eaf="0.35"),
    dict(group_id="IS_XDATA_G0022",dataset="GCST90104545",
       qc_status="PALINDROMIC_REVIEW",ref="G",alt="C",variant_id="12:111:G:C",
       alt_effect_beta="0.2",se="0.1",p="0.04",alt_effect_eaf="0.45")]
   with gzip.open(p,"wt") as h:
    w=csv.DictWriter(h,delimiter="\t",fieldnames=fields);w.writeheader();w.writerows(rows)
   out=source_variants(p)
   self.assertEqual(set(out["GCST90104545"]),{"12:100:A:G"})
   self.assertAlmostEqual(out["GCST90104545"]["12:100:A:G"][0],2.)
 def test_plink_dosage_length_guardrail(self):
  with tempfile.TemporaryDirectory() as t:
   f=Path(t)/"a.traw"
   f.write_text("CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT\tS1\tS2\n12\t12:1:A:G\t0\t1\tA\tG\t0\t1\n")
   rows=list(traw_rows(f))
   self.assertEqual(len(rows),1)
   self.assertEqual(len(rows[0][1]),2)
   f.write_text("CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT\tS1\tS2\n12\t12:1:A:G\t0\t1\tA\tG\t0\n")
   with self.assertRaisesRegex(ValueError,"Incomplete"):
    list(traw_rows(f))
if __name__=="__main__":unittest.main()
