"""Positional nearest-gene mapping is explicitly non-causal for G0022 pilot."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from annotate_g0022_pilot_signal_genes import build
def save(path,data):
 with path.open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t");w.writeheader();w.writerows(data)

class AnnotationTest(unittest.TestCase):
 def test_gene_overlap_and_distance_not_causal(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t)/"root";root.mkdir()
   base=Path(t)/"base";base.mkdir()
   gates=[{"signal":"S1","scientific_gate":"EXPLORATORY_SIGNAL_ONLY_NEEDS_EFFECTIVE_N_AND_FULL_LOCI"}]
   fits=[{"signal":"S1","hypothetical_n":"256274","top_pip_variant":"12:150:A:G","top_pip":"0.4"}]
   genes=[{"group_id":"IS_XDATA_G0022","gene_id":"ENSG1",
          "gene_symbol":"GENE1","biotype":"protein_coding","gene_start":"100","gene_end":"200"},
          {"group_id":"IS_XDATA_G0022","gene_id":"ENSG2",
          "gene_symbol":"GENE2","biotype":"lincRNA","gene_start":"250","gene_end":"350"}]
   save(root/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv",gates)
   save(root/"G0022_SUSIE_RSS_PILOT_SUMMARY.tsv",fits)
   save(base/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv",genes)
   report=build(root,base)
   self.assertEqual(report["gene_link_rows"],2)
   with (root/"G0022_PILOT_TOP_SNP_NEAREST_GENES.tsv").open() as f:
    r=list(csv.DictReader(f,delimiter="\t"))
   self.assertEqual(r[0]["gene_body_distance_bp"],"0")
   self.assertEqual(r[1]["gene_body_distance_bp"],"100")
   self.assertEqual(r[0]["overlaps_gene_body"],"1")
   self.assertIn("NOT_CAUSAL",r[0]["annotation"])
if __name__=="__main__":unittest.main()
