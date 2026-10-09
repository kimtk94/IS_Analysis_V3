import csv
import tempfile
import unittest
from pathlib import Path
from scripts.audit_is_existing_outputs import audit, FILES

class TestExistingOutputsAudit(unittest.TestCase):
    def test_missing_sources_and_gates(self):
        with tempfile.TemporaryDirectory() as d:
            result=audit(Path(d))
            self.assertTrue(all(s["status"]=="MISSING" for s in result["sources"].values()))
            self.assertEqual(result["gates"]["variant_canonical_qc"], "BLOCK_MISSING_BUILD_AND_SE")
    def test_h4_is_not_promoted_to_causal(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/FILES["signal_coloc"]
            path.parent.mkdir(parents=True)
            with path.open("w") as f:
                f.write("locus\tPP.H4.abf\nL001\t0.98\n")
            result=audit(Path(d))
            self.assertEqual(result["loci"]["L001"]["signal_coloc_evidence_status"], "PROVISIONAL_LD_UNVERIFIED")
            self.assertEqual(result["loci"]["L001"]["signal_coloc_max_h4_observed"], 0.98)

if __name__=="__main__": unittest.main()
