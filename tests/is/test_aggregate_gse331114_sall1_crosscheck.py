#!/usr/bin/env python3
"""Synthetic SALL1 crosscheck aggregation tests; no GEO or clinical data."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from aggregate_gse331114_sall1_crosscheck import (
    EXPECTED,MARKER_FIELDS,KEY_LABELS,collect,descriptive,run
)


def fixtures(root):
    for (g,i),gsm in EXPECTED.items():
        folder=root/f"is-celltypist-reference-pilot-{g}-{i}"
        folder.mkdir(parents=True)
        counts={
            "SALL1_any":20,
            "SALL1_plus_one_canonical":12,
            "SALL1_TMEM119_P2RY12_triple":8,
            "TMEM119_P2RY12_pair_regardless_SALL1":30,
            "at_least_two_FCN1_CD14_S100A8":90,
            "both_SALL1_canonical_and_peripheral":5,
        }
        def marker_row(label,n):
            if n<1:
                return {"model_provisional_label":label,"n_barcodes":0,
                        "counts":None,"rates":None}
            adjusted={key:min(value,n) for key,value in counts.items()}
            return {"model_provisional_label":label,"n_barcodes":n,
                    "counts":adjusted,
                    "rates":{k:round(v/n,6) for k,v in adjusted.items()}}
        items=[marker_row("Monocytes",600),
               marker_row("Macrophages",200),
               marker_row("Mono-mac",0),
               marker_row("T cells",150),
               marker_row("DC",50),
               marker_row("ALL_PREDICTED_MY_MARKERS",800)]
        x={
            "GSM":gsm,"group":g,"group_index_0based":i,
            "source":"NCBI_GEO_GSE331114",
            "status":"REFERENCE_PILOT_PROVISIONAL_COMPLETE",
            "original_source_SHA256_match":True,
            "source_cell_labels_confirmed":False,
            "independent_donor_identity_confirmed":False,
            "pseudobulk_DE":"NOT_RUN",
            "model_sha256":"a"*64,
            "qc_barcodes_passed":1000,
            "SALL1_per_GSM_orthogonal_RNA_markers":{
                "scope":"PUBLIC_GEO_PER_GSM_ORTHOGONAL_BINARY_RNA_DETECTION",
                "paper_SALL1_high_classifier_reproduced":False,
                "celltypes_validated":False,
                "myeloid_resident_versus_infiltrating_resolved":False,
                "n_total_barcodes":1000,
                "model_categories":items,
            },
        }
        (folder/"reference_pilot.json").write_text(json.dumps(x),encoding="utf-8")


class EightSall1CrosscheckTests(unittest.TestCase):
    def test_exact_8_samples(self):
        self.assertEqual(len(EXPECTED),8)
        self.assertEqual(set(x[0] for x in EXPECTED),{"AIS","Control"})
        self.assertEqual(len(MARKER_FIELDS),6)

    def test_reconcile_all_eight_stroke_samples(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fixtures(root)
            rows,sha=collect(root)
            self.assertEqual(len(rows),8*len(KEY_LABELS))
            self.assertEqual(sha,"a"*64)
            self.assertEqual(len(descriptive(rows)),len(MARKER_FIELDS)*len(KEY_LABELS))
            self.assertTrue(all(r["p_value"]=="NOT_TESTED" for r in descriptive(rows)))

    def test_group_and_hold_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/"input";out=Path(d)/"out"
            fixtures(root)
            report=run(root,out)
            self.assertEqual(report["status"],"SALL1_8SAMPLE_DESCRIPTIVE_CROSSCHECK_PASS")
            self.assertEqual(report["actual_cell_type_validation"],"HOLD")
            self.assertEqual(report["pseudobulk_DEG"],"NOT_RUN")
            self.assertTrue((out/"GSM_label_SALL1_marker_detection.tsv").exists())

    def test_missing_GSM_fails(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fixtures(root)
            (root/"is-celltypist-reference-pilot-AIS-2"/"reference_pilot.json").unlink()
            with self.assertRaisesRegex(ValueError,"Missing"):
                collect(root)

    def test_false_sall1_threshold_claim_fails(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fixtures(root)
            p=root/"is-celltypist-reference-pilot-AIS-0"/"reference_pilot.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["SALL1_per_GSM_orthogonal_RNA_markers"]["paper_SALL1_high_classifier_reproduced"]=True
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"orthogonal SALL1"):
                collect(root)

    def test_model_sha_change_fails(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fixtures(root)
            p=root/"is-celltypist-reference-pilot-Control-0"/"reference_pilot.json"
            x=json.loads(p.read_text(encoding="utf-8"));x["model_sha256"]="b"*64
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"identical"):
                collect(root)

    def test_denominator_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fixtures(root)
            p=root/"is-celltypist-reference-pilot-Control-1"/"reference_pilot.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["SALL1_per_GSM_orthogonal_RNA_markers"]["n_total_barcodes"]=995
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"denominators"):
                collect(root)

    def test_rate_count_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);fixtures(root)
            p=root/"is-celltypist-reference-pilot-Control-1"/"reference_pilot.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["SALL1_per_GSM_orthogonal_RNA_markers"]["model_categories"][0]["rates"]["SALL1_any"]=0.8
            p.write_text(json.dumps(x),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"rate/count"):
                collect(root)


if __name__=="__main__":
    unittest.main()
