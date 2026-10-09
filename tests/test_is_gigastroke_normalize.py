"""GIGASTROKE harmonizer: checksum, effect allele orientation, QC, no overwrites."""
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/normalize_gigastroke_verified.py"
HEAD="chromosome\tbase_pair_location\teffect_allele_frequency\tbeta\tstandard_error\tp_value\todds_ratio\tci_lower\tci_upper\teffect_allele\tother_allele\n"
class NormalizeTest(unittest.TestCase):
 def test_grch37_unordered_allele_pair_and_checksum_gate(self):
  with tempfile.TemporaryDirectory() as t:
   base=Path(t);raw=base/"raw";raw.mkdir();out=base/"out"
   name="GCST12345";f=raw/f"{name}_buildGRCh37.tsv.gz"
   with gzip.open(f,"wt") as h:
    h.write(HEAD)
    h.write("2\t100\t0.4\t-0.2\t0.03\t5e-9\t0.82\t.\t.\tT\tC\n")
    h.write("2\t101\t0.5\t0.4\t0.05\t0\t1.49\t.\t.\tG\tA\n")
    h.write("2\t102\t1.2\t0.1\t0.03\t0.5\t1.1\t.\t.\tG\tA\n")
   checksum=hashlib.md5(f.read_bytes()).hexdigest()
   manifest=base/"manifest.tsv"
   rows=[dict(accession=name,phenotype="AIS",ancestry_verified="EUR",
              source_total_n=1000,source_md5=checksum,genome_build="GRCh37")]
   with manifest.open("w") as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
   cmd=[sys.executable,str(SCRIPT),"--manifest",str(manifest),"--raw",str(raw),"--out",str(out),"--accession",name]
   p=subprocess.run(cmd,capture_output=True,text=True)
   self.assertEqual(p.returncode,0,p.stderr)
   report=json.loads(p.stdout)
   self.assertEqual(report["counts"]["total"],3)
   self.assertEqual(report["counts"]["valid"],2)
   self.assertEqual(report["counts"]["invalid"],1)
   self.assertEqual(report["counts"]["p_underflow_zero"],1)
   self.assertEqual(report["counts"]["gws_variants"],2)
   with gzip.open(report["normalized_file"],"rt") as h:records=list(csv.DictReader(h,delimiter="\t"))
   self.assertEqual(records[0]["variant_pair_id"],"2:100:C:T")
   self.assertEqual(records[0]["effect_allele"],"T")
   self.assertEqual(records[0]["ref_alt_status"],"UNKNOWN_REQUIRES_GENOME_REFERENCE")
   again=subprocess.run(cmd,capture_output=True,text=True)
   self.assertNotEqual(again.returncode,0)
   self.assertIn("Refusing to overwrite",again.stderr)
if __name__=="__main__":unittest.main()
