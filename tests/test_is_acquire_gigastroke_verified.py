"""Acquisition plans must not perform network writes without --execute."""
import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/acquire_gigastroke_verified.py"
class AcquisitionTest(unittest.TestCase):
 def test_plan_only_and_reject_unlisted(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t);output=r/"results";output.mkdir()
   rows=[dict(accession="GCST90104540",phenotype="AIS",ancestry_verified="EUR",
      source_md5="a"*32,source_url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/x/file.gz")]
   with (r/"manifest.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
   d=r/"not-created"
   cmd=[sys.executable,str(SCRIPT),"--manifest",str(r/"manifest.tsv"),
        "--dest",str(d),"--accession","GCST90104540"]
   p=subprocess.run(cmd,capture_output=True,text=True)
   self.assertEqual(p.returncode,0,p.stderr)
   self.assertEqual(json.loads(p.stdout)["status"],"PLAN_ONLY")
   self.assertFalse(d.exists())
   bad=subprocess.run(cmd[:-1]+["GCST90104541"],capture_output=True,text=True)
   self.assertNotEqual(bad.returncode,0)
   self.assertFalse(d.exists())
if __name__=="__main__":unittest.main()
