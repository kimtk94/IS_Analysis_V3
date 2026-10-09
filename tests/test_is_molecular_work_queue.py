"""Ensure no discovered gene is dropped from the molecular-screen work queue."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/build_molecular_work_queue.py"

def save(path, data):
    with path.open("w", newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t")
        w.writeheader();w.writerows(data)

class TestMolecularQueue(unittest.TestCase):
    def test_unassessed_genes_remain_in_work_queue(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            save(root/"IS_LOCUS_GENE_EVIDENCE_V2.tsv",[
                dict(group_id="G1",gene_id="ENSG01.1",gene_id_stable="ENSG01",gene_symbol="ALPHA",biotype="protein_coding",group_phenotypes="AIS",group_has_gws="1",group_best_p="1e-10",coloc_status="LOCUS_AND_GENE_ID_MATCHED"),
                dict(group_id="G1",gene_id="ENSG02.1",gene_id_stable="ENSG02",gene_symbol="BETA",biotype="lincRNA",group_phenotypes="AIS",group_has_gws="1",group_best_p="1e-10",coloc_status="NOT_TESTED_IN_THIS_GROUP"),
                dict(group_id="",gene_id="",gene_id_stable="",gene_symbol="ANCHOR_ONLY",biotype="",group_phenotypes="",group_has_gws="",group_best_p="",coloc_status="NOT_TESTED_IN_THIS_GROUP")])
            save(root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv",[
                dict(group_id="G1",chr="1",window_start="100",window_end="300",lead_variant="1:200:A:G",has_gws="1",phenotypes="AIS")])
            save(root/"LEAD_ALLELE_REFERENCE_AUDIT.tsv",[dict(group_id="G1",status="EXACT_REF_ALT_PRESENT")])
            p=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root)],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            report=json.loads((root/"IS_MOLECULAR_WORK_QUEUE_SUMMARY.json").read_text())
            self.assertEqual(report["gene_region_tasks"],2)
            self.assertEqual(report["new_eqtl_candidates"],1)
            self.assertEqual(report["legacy_coloc_reassess"],1)
            self.assertEqual(report["genotype_ready_regions"],0)
            with (root/"IS_ALL_CANDIDATE_MOLECULAR_WORK_QUEUE.tsv").open() as f:
                rows=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual({r["gene_symbol"] for r in rows},{"ALPHA","BETA"})
            self.assertTrue(all(r["ld_finemap_gate"].startswith("BLOCKED") for r in rows))
if __name__=="__main__":
    unittest.main()
