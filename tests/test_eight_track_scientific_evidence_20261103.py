"""Offline contracts for Nov 3 CI source tasks; run by GitHub Actions only."""
import importlib.util
import json
import pathlib
import struct
import unittest
from unittest.mock import patch

ROOT=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("gate",ROOT/"scripts/research/validate_eight_tracks_20261103.py")
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)

class Contracts(unittest.TestCase):
    def test_exact_eight(self):
        self.assertEqual(len(gate.TASKS),8)
        self.assertEqual(set(gate.TASKS),set(gate.RUNNERS))
        self.assertEqual(set(gate.TASKS),set(gate.NEXT))
    def test_all_biological_gates_remain_hold(self):
        self.assertNotIn("PASS",gate.NEXT.values())
        self.assertIn("HOLD",gate.execute.__code__.co_consts)
    def test_bh_expected(self):
        vals=gate._bh([.01,.04,.03,.002])
        self.assertEqual([round(x,3) for x in vals],[.02,.04,.04,.008])
    def test_bh_ties(self):
        vals=gate._bh([.1,.1,.1])
        self.assertEqual([round(x,6) for x in vals],[.1,.1,.1])
    def test_bh_zero(self):
        self.assertEqual(gate._bh([0.,1.]),[0.,1.])
    def test_nifti_spacing(self):
        buf=bytearray(352);struct.pack_into("<I",buf,0,348)
        struct.pack_into("<8f",buf,76,1,0.8,0.8,3.0,0,0,0,0)
        self.assertEqual([round(x,3) for x in gate._nifti_header(bytes(buf),False)],[.8,.8,3.])
    def test_nifti_invalid_spacing(self):
        buf=bytearray(352);struct.pack_into("<I",buf,0,348)
        struct.pack_into("<8f",buf,76,1,0,1,3,0,0,0,0)
        with self.assertRaises(gate.IntegrityFail):gate._nifti_header(bytes(buf),False)
    def test_unknown_source_host_fails_closed(self):
        with self.assertRaises(gate.Hold):gate.fetch("https://untrusted.invalid/file")
    def test_20_loci_are_in_source_repo(self):
        import csv
        p=ROOT/"config/muscle/stage1_seed_loci.tsv"
        with p.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
        self.assertEqual(len([r for r in rows if r["cohort_type"]=="RT"]),20)
        self.assertEqual(len(set(r["rsid"] for r in rows)),20)
    def test_hcc_conflicting_parents_rejected(self):
        rows=[{"aliquot_id":"a","sample_id":"s1","case_id":"c1","sample_type":"Primary Tumor"},
              {"aliquot_id":"a","sample_id":"s2","case_id":"c2","sample_type":"Normal"}]
        with patch.object(gate,"fetch",return_value=(json.dumps({"data":{"biospecimenPerStudy":rows}}).encode(),{})):
            with self.assertRaises(gate.IntegrityFail):gate._pdc("PDC000198")
    def test_no_research_execution_off_ci(self):
        with patch.dict("os.environ",{"GITHUB_ACTIONS":"false"}):
            with self.assertRaises(SystemExit):gate.execute("is")
    def test_source_manifest_never_claims_full_raw_qc(self):
        self.assertIn("whole_file_checksum_verified",gate.is_task.__code__.co_consts)
        self.assertEqual(gate.NEXT["skin"].startswith("Whole-file checksum"),True)

if __name__=="__main__":unittest.main()
