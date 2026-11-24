#!/usr/bin/env python3
"""Negative controls for donor metadata: GSM IDs != participant IDs."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_gse331114_donor_metadata import public_metadata_summary

def data(with_donor=False,duplicate=False):
    samples={}
    for i in range(8):
        group="AIS" if i<4 else "Control"
        chars=["group: "+group]
        if with_donor:
            donor="P1" if duplicate and i in (0,1) else "P"+str(i)
            chars.append("donor id: "+donor)
        if i<2:
            chars.extend(["sex: female","age: 65"])
        samples["GSM"+str(100+i)]={
            "accession":"GSM"+str(100+i),
            "title":group+"_brain"+str(i),
            "characteristics":chars,
        }
    return {"accession":"GSE331114"},samples

class TestDonorGates(unittest.TestCase):
    def test_missing_donor_does_not_imply_independence(self):
        info,samples=data()
        result=public_metadata_summary(info,samples)
        self.assertEqual(result["independent_donor_evidence_status"],"NO_EXPLICIT_DONOR_IDS_IN_GEO")
        self.assertFalse(result["biological_donor_independence_verified"])
        self.assertFalse(result["pseudobulk_DEG_permitted_from_this_metadata"])
        self.assertEqual(result["available_clinical_field_GSM_counts"]["age"],2)

    def test_explicit_unique_donor_still_not_independent_verified(self):
        info,samples=data(with_donor=True)
        result=public_metadata_summary(info,samples)
        self.assertEqual(result["n_distinct_explicit_ids_without_values"],8)
        self.assertFalse(result["biological_donor_independence_verified"])
        self.assertEqual(result["independent_donor_evidence_status"],"EXPLICIT_IDS_UNIQUE_IN_GEO_NOT_BIOLOGICALLY_VERIFIED")

    def test_repeated_donor_flags_possible_pseudoreplication(self):
        info,samples=data(with_donor=True,duplicate=True)
        result=public_metadata_summary(info,samples)
        self.assertEqual(result["independent_donor_evidence_status"],"POSSIBLE_REPEATED_DONOR_IDS")
        self.assertEqual(result["n_distinct_explicit_ids_without_values"],7)

    def test_wrong_accession_fails(self):
        info,samples=data()
        with self.assertRaises(ValueError):
            public_metadata_summary({"accession":"GSE1"},samples)

    def test_missing_one_GSM_fails(self):
        info,samples=data()
        samples.pop("GSM100")
        with self.assertRaises(ValueError):
            public_metadata_summary(info,samples)

    def test_raw_donor_values_not_exported(self):
        info,samples=data(with_donor=True)
        result=public_metadata_summary(info,samples)
        serial=__import__("json").dumps(result)
        self.assertNotIn("P0",serial)
        self.assertNotIn("P1",serial)

if __name__=="__main__":
    unittest.main()
