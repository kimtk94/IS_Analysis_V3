import csv
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from render_is_priority4_reference_ld_figure import render

class ReferenceLDGraphicsTest(unittest.TestCase):
 def test_descriptive_svg_tracks_two_loci_and_causality_disclaimer(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);sp=d/"s.tsv";gp=d/"g.tsv";out=d/"scatter.svg"
   with sp.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["locus","EAS_r2_with_GWAS_lead","EUR_r2_with_GWAS_lead"],delimiter="\t")
    w.writeheader()
    for locus in ("BBJ_IS_L001","BBJ_IS_L002"):
     w.writerow(dict(locus=locus,EAS_r2_with_GWAS_lead=".4",EUR_r2_with_GWAS_lead=".5"))
   with gp.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["locus","symbol","EAS_reference_r2","EUR_reference_r2"],delimiter="\t")
    w.writeheader()
    w.writerow(dict(locus="BBJ_IS_L001",symbol="FGF5",EAS_reference_r2=".4",EUR_reference_r2=".5"))
    w.writerow(dict(locus="BBJ_IS_L002",symbol="CALHM2",EAS_reference_r2=".8",EUR_reference_r2=".9"))
   self.assertEqual(render(sp,gp,out),2)
   root=ET.parse(out).getroot()
   texts=" ".join(e.text or "" for e in root.iter())
   self.assertIn("Not original BBJ/GTEx",texts)
   self.assertIn("CALHM2",texts)
   self.assertIn("FGF5",texts)
   self.assertEqual(sum(1 for e in root.iter() if e.tag.endswith("circle")),2)
if __name__=="__main__":unittest.main()
