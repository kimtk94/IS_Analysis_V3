"""AS-only marker audit distinguishes absent positions and mismatched alleles."""
import csv,gzip,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_as_only_source import run

def write(path,rows,gz=False):
    op=gzip.open if gz else open
    with op(path,"wt") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
class OnlyMarkerTest(unittest.TestCase):
 def test_alternative_explanations_are_kept_separate(self):
  with tempfile.TemporaryDirectory() as td:
   base=Path(td);root=base/"root";root.mkdir()
   region=root/"GCST90104544_AS_signal1_111629389";region.mkdir()
   as_rows=[]
   for pos in (100,200,300):
    as_rows.append({"variant_id":f"12:{pos}:A:G","pos":str(pos),
      "gwas_z":"1.5","gwas_p":".1","reference_alt_af":".3","gwas_alt_eaf":".31"})
   write(region/"variants.tsv",as_rows)
   ais=base/"AIS.gz"
   write(ais,[{"chr":"12","pos":"100","effect_allele":"A",
     "other_allele":"G","p":"0.7"},
     {"chr":"12","pos":"200","effect_allele":"A",
      "other_allele":"T","p":"0.1"}],True)
   prior=base/"matched.gz"
   write(prior,[{"group_id":"OTHER","dataset":"OTHER",
     "qc_status":"MATCH_ALT_EFFECT","variant_id":"12:1:A:G"}],True)
   s=run(root,ais,prior)
   self.assertEqual(s["as_only_marker_rows"],3)
   self.assertEqual(s["source_presence_status"],{
     "ALLELE_PAIR_PRESENT_BUT_FAILED_JOIN_OR_QC":1,
     "POSITION_PRESENT_DIFFERENT_ALLELES":1,
     "NO_AIS_CANONICAL_POSITION":1})
if __name__=="__main__":unittest.main()
