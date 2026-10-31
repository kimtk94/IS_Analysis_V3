#!/usr/bin/env python3
"""Synthetic positive/negative controls for HCC Stage47 strict source-ID bridge."""
import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage47_official_aliquot_bridge import alias_evidence,numerical_bridge,review

class Stage47Tests(unittest.TestCase):
    def test_exact_330_official_aliquot_ids_pass(self):
        labels=[f"al-{i:03d}" for i in range(330)]
        records=[{"aliquot_id":x,"aliquot_submitter_id":f"S{i:03d}",
                  "sample_id":f"sample{i:03d}","sample_submitter_id":f"CS{i:03d}"}
                 for i,x in enumerate(labels)]
        records.append({"aliquot_id":"ref","aliquot_submitter_id":"POOL"})
        x=alias_evidence(labels,records,{"ref"},set(labels)|{"ref"})
        self.assertTrue(x["strict_330_official_aliquot_id_bridge_verified"])
        self.assertEqual(x["alias_namespace_results"]["aliquot_id"]["exact_unique_header_matches_n"],330)

    def test_duplicate_aliquot_id_in_source_fails(self):
        labels=[f"al-{i:03d}" for i in range(330)]
        records=[{"aliquot_id":x,"aliquot_submitter_id":f"S{i:03d}"} for i,x in enumerate(labels)]
        records[-1]["aliquot_id"]=records[0]["aliquot_id"]
        records.append({"aliquot_id":"ref"})
        x=alias_evidence(labels,records,{"ref"},set(labels)|{"ref"})
        self.assertFalse(x["strict_330_official_aliquot_id_bridge_verified"])

    def test_sample_level_alias_not_enough_for_aliquot_certification(self):
        labels=[f"sample{i:03d}" for i in range(330)]
        records=[{"aliquot_id":f"al-{i:03d}","sample_id":s} for i,s in enumerate(labels)]
        records.append({"aliquot_id":"ref"})
        design={f"al-{i:03d}" for i in range(330)}|{"ref"}
        x=alias_evidence(labels,records,{"ref"},design)
        self.assertEqual(x["alias_namespace_results"]["sample_id"]["exact_unique_header_matches_n"],330)
        self.assertFalse(x["strict_330_official_aliquot_id_bridge_verified"])

    def test_numeric_mapping_uses_independent_heldout_values(self):
        rng=np.random.default_rng(47)
        cdap=rng.normal(size=(24,330))
        pdc=np.round(cdap,4)
        headers=[f"sample{i:03d} Log Ratio" for i in range(330)]
        labels,proof=numerical_bridge(pdc,cdap,headers)
        self.assertEqual(labels,[f"sample{i:03d}" for i in range(330)])
        self.assertEqual(proof["independent_permutation_agreement_n"],330)
        self.assertEqual(proof["independent_four_decimal_equal_cells_n"],12*330)
        with self.assertRaises(ValueError):
            bad=pdc.copy()
            bad[1,0]+=0.1
            numerical_bridge(bad,cdap,headers)

    def test_review_remains_hold_until_official_match(self):
        data={"schema":"hcc-stage47-official-PDC-aliquot-bridge-v1",
              "source":{"status":"PDC_NUMERIC_TO_CDAP_LABEL_TO_OFFICIAL_ID_EXACT_AUDITED",
              "original_labels_or_values_exported":False,
              "NODE_RNA_same_biopsy_certified":False,
              "official_identity_bridge":{"strict_330_official_aliquot_id_bridge_verified":False}}}
        self.assertTrue(review(data)["G30_B"].startswith("HOLD"))
        self.assertTrue(review(data)["G30_D"].startswith("STOP"))
        data["source"]["original_labels_or_values_exported"]=True
        with self.assertRaises(ValueError):review(data)

if __name__=="__main__":
    unittest.main()
