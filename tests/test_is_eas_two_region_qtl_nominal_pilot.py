"""Regression: public nominal eQTL Catalogue region source and GRCh37/38 gates."""
import gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/pilot_is_eas_two_region_public_nominal_qtl.py"
spec=importlib.util.spec_from_file_location("easqtlpilot",p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)
class GoodLifter:
 def convert_coordinate(self,c,pos):
  return [(c,pos+100,"+",42)]
class BadLifter:
 def convert_coordinate(self,c,pos):
  return [(c,pos+100,"-",42)]
class TestEASQTL(unittest.TestCase):
 def test_build_and_chain_safety(self):
  r={"group_id":"IS_XDATA_G0025","chr":"17","start":"7095198","end":"7128727"}
  chrom,start,end,control=m.lift_interval(GoodLifter(),r)
  self.assertEqual((chrom,start,end),("17",6595298,7628827))
  with self.assertRaises(ValueError):m.lift_interval(BadLifter(),r)
 def test_nominal_source_19_column_genes_not_coloc(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"qtl.tsv.gz"
   row=["ENSG000001.11","17","7120000","A","G",
        "chr17_7120000_A_G","40",".15",".001","0.1",".05",
        "SNP","32","200",".95","ENSG000001.11","ENSG000001.11","24","rs10"]
   with gzip.open(p,"wt") as f:f.write("\t".join(row)+"\n")
   counts,genes=m.audit_slice(p,"17",6000000,7500000,{"ENSG000001"})
   self.assertEqual(counts["rows"],1)
   self.assertEqual(counts["in_positional_genes"],1)
   self.assertEqual(len(genes),1)
   self.assertAlmostEqual(genes["ENSG000001"]["minp"],.001)
 def test_qtl_wrong_variant_tuple_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"qtl.tsv.gz"
   row=["ENSG000001","17","7120000","A","G",
        "chr17_7120000_G_A","40",".15",".001","0.1",".05",
        "SNP","32","200",".95","ENSG000001","ENSG000001","24","rs10"]
   with gzip.open(p,"wt") as f:f.write("\t".join(row)+"\n")
   with self.assertRaises(ValueError):
    m.audit_slice(p,"17",6000000,7500000,{"ENSG000001"})
 def test_no_retroactive_raw_overwrite(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"exists.tsv.gz";p.write_text("not a qtl")
   with self.assertRaises(FileExistsError):m.fetch("QTD000131","17:10-20",p)
if __name__=="__main__":unittest.main()
