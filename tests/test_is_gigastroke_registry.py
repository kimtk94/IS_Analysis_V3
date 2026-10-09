"""A catalog metadata UNKNOWN ancestry must never be silently relabeled."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/audit_gigastroke_full_registry.py"

class TestStudyInventory(unittest.TestCase):
    def test_unknown_kept_unknown_without_local_canonical(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)
            manifest=base/"src.tsv"
            columns=["accession","phenotype_class","ancestry_class","sample_size"]
            with manifest.open("w",newline="") as h:
                w=csv.DictWriter(h,delimiter="\t",fieldnames=columns)
                w.writeheader()
                w.writerow(dict(accession="GCST123",phenotype_class="AIS",ancestry_class="UNKNOWN",sample_size="100"))
                w.writerow(dict(accession="GCST124",phenotype_class="LAS",ancestry_class="HIS",sample_size="150"))
            root=base/"out"
            root.mkdir()
            p=subprocess.run([sys.executable,str(SCRIPT),"--registry",str(manifest),"--root",str(root),"--data",str(base/"data")],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            data=json.loads((root/"GIGASTROKE_STUDY_INVENTORY_SUMMARY.json").read_text())
            self.assertEqual(data["registry_accessions"],2)
            self.assertEqual(data["local_eas_canonical"],0)
            with (root/"GIGASTROKE_FULL_STUDY_INVENTORY.tsv").open() as h:
                rows=list(csv.DictReader(h,delimiter="\t"))
            self.assertEqual(rows[0]["population_resolution"],"UNRESOLVED_SOURCE_METADATA")
            self.assertEqual(rows[1]["population_resolution"],"HIS_REGISTRY_LABEL_ONLY")
if __name__=="__main__":
    unittest.main()
