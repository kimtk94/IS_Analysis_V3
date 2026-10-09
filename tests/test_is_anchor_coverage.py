"""End-to-end synthetic anchor comparison contract."""
import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/compare_anchor_gene_coverage.py"
class TestAnchorCoverage(unittest.TestCase):
 def test_all_anchors_preserved_and_absence_explicit(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td)
   with (root/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=["group_id","gene_symbol","gene_id","lead_distance_bp"],delimiter="\t");w.writeheader()
    w.writerow(dict(group_id="G1",gene_symbol="FGF5",gene_id="ENSG1",lead_distance_bp=0))
    w.writerow(dict(group_id="G1",gene_symbol="OTHER",gene_id="ENSG2",lead_distance_bp=100))
   with (root/"LEGACY_ANCHOR_CANDIDATES.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=["gene","locus"],delimiter="\t");w.writeheader()
    w.writerow(dict(gene="FGF5",locus="L001"));w.writerow(dict(gene="ALDH2",locus="L003"))
   p=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root)],capture_output=True,text=True)
   self.assertEqual(p.returncode,0,p.stderr)
   d=json.loads((root/"IS_CANDIDATE_COVERAGE_SUMMARY.json").read_text())
   self.assertEqual(d["legacy_anchor_count"],2)
   self.assertEqual(d["anchors_overlapping_expanded_windows"],1)
   self.assertEqual(d["unique_positional_genes"],2)
if __name__=="__main__":unittest.main()
