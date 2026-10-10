"""PLINK EAS regional clump runner default never executes; contract QC."""
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from run_alcohol_eas_regional_clumping import plan,run
class ClumpingRunner(unittest.TestCase):
    def test_complete_sensitivity_plan_and_no_implicit_execution(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/"result";ref=Path(td)/"ref"
            tasks=plan(root,ref)
            self.assertEqual(len(tasks),20)
            self.assertEqual(sum("--clump" in t["command"] for t in tasks),16)
            self.assertEqual(sum("--freq" in t["command"] for t in tasks),4)
            self.assertEqual(sum("--maf" in t["command"] for t in tasks),12)
            self.assertEqual(run(root,ref,False)["task_count"],20)
            with self.assertRaises(FileNotFoundError):
                run(root,ref,True)
if __name__=="__main__":unittest.main()
