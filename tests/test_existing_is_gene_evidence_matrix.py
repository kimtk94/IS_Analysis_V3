"""Scientific guard tests for 2,225-gene IS evidence joins."""
from __future__ import annotations
import importlib.util
import sys
import unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/is/build_existing_is_gene_evidence_matrix.py"
spec=importlib.util.spec_from_file_location("is_2225_evidence",SOURCE)
module=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=module
spec.loader.exec_module(module)

class TestExistingGeneEvidence(unittest.TestCase):
    def test_stable_id_drops_only_ensembl_version(self):
        self.assertEqual(module.stable("ENSG00000107954.12"),"ENSG00000107954")
        self.assertEqual(module.stable("ENSG00000107954"),"ENSG00000107954")

    def test_qualified_key_uses_ensembl_id_not_gene_symbol(self):
        x={"group_id":"IS_XDATA_G0015","gene_symbol":"NEURL",
           "gene_id_stable":"ENSG00000107954"}
        y={"group_id":"IS_XDATA_G0015","gene_symbol":"NEURL1",
           "gene_id_stable":"ENSG00000107954"}
        self.assertEqual(module.exact_row_key(x),module.exact_row_key(y))

    def test_short_or_broken_input_fails_closed(self):
        with self.assertRaises((KeyError,ValueError)):
            module.build({
                "expanded":[],"windows":[],"coloc":[],"replay":[],
                "priority8":[],"literature":[],"cell_gene":[],
                "cell_feature":[],"cell_pairs":[],"universe":{"unique_gene_ids":2225}
            })

    def test_scalar_numeric_conversion(self):
        self.assertEqual(module.fnum("1e-9"),1e-9)
        self.assertIsNone(module.fnum(""))

if __name__=="__main__":
    unittest.main()
