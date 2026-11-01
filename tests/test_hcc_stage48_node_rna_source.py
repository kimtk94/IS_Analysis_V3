#!/usr/bin/env python3
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage48_node_rna_source import inspect_html,public_href,inspect_manifest_rows,review

class Stage48Tests(unittest.TestCase):
    def test_script_only_html_does_not_prove_matrix(self):
        x=inspect_html('<html>OEP000321<script src="/a.js"></script></html>')
        self.assertTrue(x["landing_embeds_project_accession"])
        self.assertFalse(x["matrix_file_or_download_verified"])
        self.assertEqual(x["literal_public_download_anchor_n"],0)
    def test_explicit_public_download_is_still_uncertified(self):
        x=inspect_html('<html>OEP000321<a href="/download/node/data/public/OED063528">file</a></html>')
        self.assertEqual(x["literal_public_download_anchor_n"],1)
        self.assertFalse(x["linked_file_ids_proven_to_belong_to_project"])
        self.assertFalse(x["matrix_file_or_download_verified"])
    def test_disallow_external_private_and_guessed_urls(self):
        self.assertFalse(public_href("https://other.example/download/node/data/public/OED063528"))
        self.assertFalse(public_href("https://www.biosino.org/download/node/data/private/OED063528"))
        self.assertFalse(public_href("https://www.biosino.org/download/node/data/public/OED063528?secret=a"))
        self.assertTrue(public_href("https://www.biosino.org/download/node/data/public/OED063528"))
    def test_optional_original_export_is_not_expression_verified(self):
        row={"Project ID":"OEP000321","Sample ID":"OES012345","Data ID":"OED063528",
             "Data Name":"RSEM_gene_expression.tsv",
             "Download Link":"https://www.biosino.org/download/node/data/public/OED063528"}
        x=inspect_manifest_rows([row])
        self.assertEqual(x["processed_candidate_rows_n"],1)
        self.assertEqual(x["node_valid_unique_OES_sample_ids_n"],1)
        self.assertFalse(x["processed_RNA_matrix_dimensions_verified"])
        self.assertFalse(x["PDC_same_biopsy_sample_mapping_verified"])
    def test_wrong_cohort_rejected(self):
        with self.assertRaises(ValueError):
            inspect_manifest_rows([{"Project ID":"OEP001000","Data ID":"OED012345"}])
    def test_missing_manifest_project_column_hold(self):
        x=inspect_manifest_rows([{"Data ID":"OED012345"}])
        self.assertEqual(x["status"],"PROJECT_COLUMN_MISSING_HOLD")
    def test_biological_gates_remain_closed(self):
        source={"project":{"http_status":200},
                "raw_file_or_expression_matrix_downloaded":False,
                "PDC_aliquot_to_NODE_RNA_sample_join_verified":False}
        r=review(source)
        self.assertTrue(r["G30_C"].startswith("HOLD"))
        self.assertTrue(r["G30_D"].startswith("STOP"))
        source["raw_file_or_expression_matrix_downloaded"]=True
        with self.assertRaises(ValueError):review(source)
if __name__=="__main__":unittest.main()
