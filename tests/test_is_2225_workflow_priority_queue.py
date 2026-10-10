import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_existing_2225_evidence_priority_queue import run

class TestQueue(unittest.TestCase):
 def test_bad_counts_fail_closed(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/"genes.tsv").write_text("gene_id_stable\nENSG000001\n")
   (p/"regions.tsv").write_text("gene_id_stable\nENSG000001\n")
   (p/"manifest.json").write_text(json.dumps({"stable_gene_ids":2225,"gene_region_pairs":2425}))
   with self.assertRaisesRegex(ValueError,"count"):
    run(p/"genes.tsv",p/"regions.tsv",p/"manifest.json",p/"out")
 def test_real_universe_if_source_mounted(self):
  p=Path("/srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1")
  if not p.is_dir(): self.skipTest("original IS server source not mounted")
  with tempfile.TemporaryDirectory() as d:
   x=run(p/"IS_2225_GENE_EVIDENCE_MATRIX.tsv",p/"IS_2425_GENE_REGION_EVIDENCE.tsv",p/"IS_GENE_EVIDENCE_MATRIX_MANIFEST.json",Path(d))
   self.assertEqual(x["total_source_unique_genes"],2225)
   self.assertEqual(x["causal_genes_established"],0)
   self.assertFalse(x["research_wide_hard_candidate_filter_applied"])
   self.assertEqual(sum(x["tier_counts"].values()),2225)
   with open(Path(d)/"IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv") as f:
    rows=list(csv.DictReader(f,delimiter="\t"))
   self.assertTrue(all(r["direct_causal_gene_status"]=="NOT_ESTABLISHED" for r in rows))
   self.assertTrue(all(r["not_assessed_is_not_negative"]=="True" for r in rows))
if __name__=="__main__":unittest.main()
