import tempfile
import unittest
from pathlib import Path
from scripts.build_is_phase10_evidence_matrix import build, HANDOFF, DRIVE
from scripts.audit_is_phase10_readiness import LD, CELL

class MatrixTests(unittest.TestCase):
    def test_unexecuted_handoff_does_not_count_as_human_validation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name, data in {
                LD:"locus\tstatus\tn_variants\nBBJ_IS_L003\tMISSING\t0\n",
                CELL:"human_gene\tdataset\tspecies\ttop_celltype\nALDH2\tGSE225948\tMouse\tpre3\n",
                HANDOFF:"component\tstatus\nSTATUS\tREADY_TO_RUN_NOT_EXECUTED\n",
                DRIVE:"component\tstatus\nLOCAL_RDS\tPRESENT\nREMOTE_RDS\tPRESENT_SIZE_MATCH\n",
            }.items():
                path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(data)
            rows=build(root)
            self.assertEqual(len(rows),4)
            aldh=next(r for r in rows if r["gene"]=="ALDH2")
            self.assertEqual(aldh["human_celltype_validation"],"NOT_ESTABLISHED")
            self.assertEqual(aldh["ld_matrix_status"],"MISSING")
            self.assertEqual(aldh["mouse_rows"],1)
            self.assertEqual(aldh["phase10_gate"],"AWAITING_COLAB_OUTPUT_VALIDATION")
if __name__=="__main__":unittest.main()
