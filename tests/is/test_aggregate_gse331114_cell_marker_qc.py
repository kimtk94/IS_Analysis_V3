#!/usr/bin/env python3
"""Synthetic eight-GSM summaries ensure missing cell labels cannot cause fake pseudobulk."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from aggregate_gse331114_cell_marker_qc import (  # noqa: E402
    EXPECTED,PANELS,load,descriptive,process,
)


def fixtures(root):
    for (group,idx),gsm in EXPECTED.items():
        dest=root/f"is-gse331114-marker-{group}-{idx}"
        dest.mkdir(parents=True)
        categories={x:0 for x in PANELS}
        categories["T_cell_marker_support"]=20+idx
        categories["UNASSIGNED"]=80-idx
        dat={
            "source":"GSE331114","GSM":gsm,"group":group,"index":idx,
            "status":"SAMPLE_CELL_QC_AND_MARKER_SUPPORT_PASS",
            "source_file_sha256_matches_prior":True,
            "not_a_celltype_annotation":True,
            "no_p_values_or_degs":True,
            "celltype_pseudobulk":"NOT_PERFORMED",
            "donor_independence":"NOT_VERIFIED_WITH_INDIVIDUAL_IDENTIFIERS",
            "cells_pass_qc":100, "cells_fail_qc":15,"total_barcodes":115,
            "pct_pass_qc":86.957,"median_mito_percent_all_barcodes":4.5,
            "median_n_genes_all_barcodes":350,"median_n_umis_all_barcodes":900,
            "marker_panel_support_unvalidated_counts_after_QC":categories,
            "metadata_characteristic_field_names":["group"],
        }
        (dest/"marker_qc.json").write_text(json.dumps(dat),encoding="utf-8")


class EightMarkerQCTests(unittest.TestCase):
    def test_per_sample_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fixtures(root)
            x=load(root)
            self.assertEqual(len(x),8)
            self.assertEqual(set(r["GSM"] for r in x),set(EXPECTED.values()))
            d=descriptive(x)
            self.assertEqual(len(d),len(PANELS))
            self.assertTrue(all(z["p_value"]=="NOT_TESTED" for z in d))

    def test_outputs_hold_inference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/"data"; out=Path(tmp)/"result"
            fixtures(root)
            process(root,out)
            evidence=json.loads((out/"evidence.json").read_text(encoding="utf-8"))
            self.assertFalse(evidence["cell_type_pseudobulk_ready"])
            self.assertEqual(evidence["cell_type_DEG"],"NOT_RUN")
            self.assertEqual(evidence["inferred_cell_type_labels"],"NONE")
            self.assertTrue(evidence["all_eight_actual_source_artifacts_validated"])
            self.assertTrue((out/"marker_group_descriptive.tsv").exists())

    def test_missing_one_artifact_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures(root)
            (root/"is-gse331114-marker-AIS-3"/"marker_qc.json").unlink()
            with self.assertRaisesRegex(ValueError,"Missing"):
                load(root)

    def test_wrong_gsm_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures(root)
            p=root/"is-gse331114-marker-Control-2"/"marker_qc.json"
            x=json.loads(p.read_text(encoding="utf-8"));x["GSM"]="GSMWRONG"
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"identity"):
                load(root)

    def test_marker_counts_must_partition_passed_cells(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures(root)
            p=root/"is-gse331114-marker-AIS-2"/"marker_qc.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["marker_panel_support_unvalidated_counts_after_QC"]["UNASSIGNED"]=1
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"category total"):
                load(root)

    def test_no_invented_cell_type_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures(root)
            p=root/"is-gse331114-marker-AIS-2"/"marker_qc.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["not_a_celltype_annotation"]=False
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"guardrail"):
                load(root)


if __name__=="__main__":
    unittest.main()
