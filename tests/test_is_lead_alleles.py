"""Synthetic integration test for allele matching and orientation flags."""
import csv, subprocess, sys, tempfile, unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/audit_lead_alleles.py"
class TestLeadAlleles(unittest.TestCase):
 def test_exact_swap_mismatch_missing(self):
  with tempfile.TemporaryDirectory() as temp:
   root=Path(temp)/"out";root.mkdir();ref=Path(temp)/"refs";ref.mkdir()
   fields=["group_id","lead_variant"]
   with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open("w") as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader()
    for n,v in enumerate(["1:10:A:G","1:20:T:C","1:30:A:G","1:40:A:T"],1):w.writerow({"group_id":str(n),"lead_variant":v})
   prefix=ref/"demo"
   Path(str(prefix)+".pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n1\t10\tx\tA\tG\n1\t20\tx\tC\tT\n1\t30\tx\tA\tC\n")
   Path(str(prefix)+".pgen").write_bytes(b"stub")
   Path(str(prefix)+".psam").write_text("#IID\nx\n")
   p=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root),"--reference",str(ref)],capture_output=True,text=True)
   self.assertEqual(p.returncode,0,p.stderr)
   with (root/"LEAD_ALLELE_REFERENCE_AUDIT.tsv").open() as h: rows=list(csv.DictReader(h,delimiter="\t"))
   self.assertEqual([r["status"] for r in rows],["EXACT_REF_ALT_PRESENT","SWAPPED_REF_ALT_PRESENT_REQUIRES_ORIENTATION","POSITION_ONLY_ALLELES_DIFFER","NO_POSITION_MATCH"])
   self.assertTrue(all(r["fine_mapping_ready"]=="NO" for r in rows))
if __name__=="__main__":unittest.main()
