"""Integration test expanded IS region manifest using synthetic PVAR fixtures."""
import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/build_expanded_region_manifest.py"
class TestExpandedManifest(unittest.TestCase):
 def test_reference_position_missingness_and_no_causal_claim(self):
  with tempfile.TemporaryDirectory() as tmp:
   base=Path(tmp);root=base/"results";root.mkdir();ref=base/"ref";ref.mkdir()
   with (root/"CROSS_DATASET_INTERVAL_GROUPS.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=["group_id","chr","start","end","lead_variant","lead_p","phenotypes","has_gws"],delimiter="\t");w.writeheader()
    w.writerow(dict(group_id="G1",chr="4",start=90,end=110,lead_variant="4:100:A:C",lead_p="1e-9",phenotypes="AIS",has_gws=1))
    w.writerow(dict(group_id="G2",chr="5",start=200,end=220,lead_variant="5:205:T:G",lead_p="8e-7",phenotypes="LAS",has_gws=0))
   stem=ref/"example"
   stem.with_suffix(".pvar").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n4\t100\trsX\tA\tG\n")
   stem.with_suffix(".pgen").write_bytes(b"stub")
   stem.with_suffix(".psam").write_text("#IID\nS1\n")
   p=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root),"--ref",str(ref)],capture_output=True,text=True)
   self.assertEqual(p.returncode,0,p.stderr)
   meta=json.loads((root/"IS_EXPANDED_EXECUTION_STATUS.json").read_text())
   self.assertEqual(meta["lead_position_present"],1)
   self.assertEqual(meta["lead_position_absent"],1)
   with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f: lines=list(csv.DictReader(f,delimiter="\t"))
   self.assertEqual(lines[0]["ld_status"],"LEAD_POSITION_FOUND_ALLELE_QC_REQUIRED")
   self.assertEqual(lines[1]["ld_status"],"NEEDS_NEW_ANCESTRY_MATCHED_REFERENCE")
   self.assertIn("BLOCKED",lines[0]["fine_mapping_status"])
if __name__=="__main__":unittest.main()
