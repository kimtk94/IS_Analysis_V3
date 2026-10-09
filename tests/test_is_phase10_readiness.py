import tempfile
import unittest
from pathlib import Path
from scripts.audit_is_phase10_readiness import audit, LD, CELL

class ReadinessTests(unittest.TestCase):
    def test_no_files_never_pass(self):
        with tempfile.TemporaryDirectory() as d:
            r=audit(Path(d))
            self.assertIn("LD_AUDIT_FILE_MISSING", r["warnings"])
            self.assertEqual(r["genes"]["FGF5"]["human_validation"], "NOT_ESTABLISHED")
    def test_species_and_ld_are_not_overclaimed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            l=root/LD;l.parent.mkdir(parents=True)
            l.write_text("locus\tstatus\tn_variants\nBBJ_IS_L001\tPASS\t100\n")
            c=root/CELL;c.parent.mkdir(parents=True)
            c.write_text("human_gene\tdataset\tspecies\ttop_celltype\nSH3PXD2A\tGSE225948\tMouse\tMg8\n")
            r=audit(root)
            self.assertFalse(r["ld_matrix_audit"]["BBJ_IS_L001"]["cohort_matched_ld_verified"])
            self.assertEqual(r["genes"]["SH3PXD2A"]["human_validation"], "NOT_ESTABLISHED")
            self.assertEqual(r["genes"]["SH3PXD2A"]["mouse_evidence"], "DESCRIPTIVE_ONLY")
if __name__ == "__main__": unittest.main()
