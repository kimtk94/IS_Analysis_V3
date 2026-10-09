"""LD extractor is planning-only without opt-in and rejects inappropriate sources."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/extract_broad_eas_ld.py"
class TestLDExtraction(unittest.TestCase):
    def test_plan_only_does_not_make_files(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            samples=root/"EAS.samples.txt"
            samples.write_text("".join(f"HG{i:05d}\n" for i in range(504)))
            inp=root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv"
            with inp.open("w",newline="") as h:
                writer=csv.DictWriter(h,fieldnames=["group_id","has_gws","lead_variant",
                    "chr","window_start","window_end"],delimiter="\t")
                writer.writeheader()
                writer.writerow(dict(group_id="G1",has_gws=1,lead_variant="2:200:A:G",chr=2,
                    window_start=100,window_end=300))
                writer.writerow(dict(group_id="G2",has_gws=0,lead_variant="2:250:A:G",chr=2,
                    window_start=100,window_end=300))
            out=root/"not_created"
            base=[sys.executable,str(SCRIPT),"--root",str(root),"--samples",str(samples),
                "--out",str(out)]
            proc=subprocess.run(base+["--group","G1"],capture_output=True,text=True)
            self.assertEqual(proc.returncode,0,proc.stderr)
            self.assertEqual(json.loads(proc.stdout)["status"],"PLAN_ONLY")
            self.assertFalse(out.exists())
            blocked=subprocess.run(base+["--group","G2"],capture_output=True,text=True)
            self.assertNotEqual(blocked.returncode,0)
            invalid=subprocess.run(base+["--group","OTHER"],capture_output=True,text=True)
            self.assertNotEqual(invalid.returncode,0)
            self.assertFalse(out.exists())
if __name__=="__main__":
    unittest.main()
