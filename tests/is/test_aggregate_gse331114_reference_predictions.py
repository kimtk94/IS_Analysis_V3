#!/usr/bin/env python3
"""Synthetic tests. Test fixtures are not real IS expression measurements."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from aggregate_gse331114_reference_predictions import (
    EXPECTED,MODEL,validate_and_collect,summarize,run,
)


def synth(root):
    for (g,i),gsm in EXPECTED.items():
        folder=root/f"is-celltypist-reference-pilot-{g}-{i}"
        folder.mkdir(parents=True)
        data={
            "GSM":gsm,"group":g,"group_index_0based":i,"source":"NCBI_GEO_GSE331114",
            "status":"REFERENCE_PILOT_PROVISIONAL_COMPLETE",
            "model":MODEL,"original_source_SHA256_match":True,
            "reference_predictions_are_provisional":True,
            "source_cell_labels_confirmed":False,
            "brain_microglia_versus_peripheral_myeloid_validated":False,
            "independent_donor_identity_confirmed":False,
            "pseudobulk_DE":"NOT_RUN",
            "reference_genes_matched_fraction":0.91,
            "model_sha256":"a"*64,
            "reference_genes_matched":6000,
            "n_predictions":1000,
            "qc_barcodes_passed":1000,
            "original_barcodes":1020,
            "top_score_median":0.8,
            "fraction_score_ge_threshold":0.75,
            "provisional_labels_top_20":[
                {"model_label":"Monocytes","n_cells":600 if g=="AIS" else 300},
                {"model_label":"Macrophages","n_cells":200 if g=="AIS" else 500},
                {"model_label":"T cells","n_cells":150},
            ],
        }
        (folder/"reference_pilot.json").write_text(json.dumps(data),encoding="utf-8")


class ReferenceCohortSynthesisTest(unittest.TestCase):
    def test_expected_samples(self):
        self.assertEqual(len(EXPECTED),8)
        self.assertEqual(EXPECTED[("AIS",0)],"GSM9739251")
        self.assertEqual(EXPECTED[("Control",3)],"GSM9739258")

    def test_eight_sample_merge(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);synth(root)
            rows,sha,n=validate_and_collect(root)
            self.assertEqual((len(rows),sha,n),(8,"a"*64,6000))
            summary=summarize(rows)
            mono=next(x for x in summary if x["model_label"]=="monocytes")
            self.assertEqual(mono["difference_of_group_medians"],0.3)
            self.assertEqual(mono["p_value"],"NOT_TESTED")

    def test_outputs_retain_hold(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/"in";out=Path(tmp)/"out";synth(root)
            run(root,out)
            evidence=json.loads((out/"evidence.json").read_text(encoding="utf-8"))
            self.assertFalse(evidence["cell_types_validated"])
            self.assertFalse(evidence["independent_donors_confirmed"])
            self.assertFalse(evidence["celltype_pseudobulk_gene_counts_ready"])
            self.assertEqual(evidence["celltype_DE"],"NOT_RUN")
            self.assertTrue((out/"per_GSM_reference_predictions.tsv").exists())

    def test_missing_GSM_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);synth(root)
            (root/"is-celltypist-reference-pilot-AIS-2"/"reference_pilot.json").unlink()
            with self.assertRaisesRegex(ValueError,"Missing"):
                validate_and_collect(root)

    def test_mismatched_model_sha_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);synth(root)
            p=root/"is-celltypist-reference-pilot-Control-0"/"reference_pilot.json"
            data=json.loads(p.read_text(encoding="utf-8"))
            data["model_sha256"]="b"*64
            p.write_text(json.dumps(data),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"model hashes"):
                validate_and_collect(root)

    def test_invalid_biological_causality_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);synth(root)
            p=root/"is-celltypist-reference-pilot-AIS-1"/"reference_pilot.json"
            data=json.loads(p.read_text(encoding="utf-8"))
            data["brain_microglia_versus_peripheral_myeloid_validated"]=True
            p.write_text(json.dumps(data),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"inference"):
                validate_and_collect(root)


if __name__=="__main__":
    unittest.main()
