"""EUR GIGASTROKE harmonization and allele-pair regression tests."""
import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/is/audit_alcohol_ais_eur_transfer.py"
spec=importlib.util.spec_from_file_location("eur_audit",SOURCE)
audit=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=audit
spec.loader.exec_module(audit)

class TestEUR(unittest.TestCase):
    def target(self):
        return {"chrom":"4","pos":"100","ref":"T","alt":"C"}

    def row(self):
        return {"accession":"GCST90104540","ancestry":"EUR","phenotype":"AIS",
                "build":"GRCh37","chr":"4","pos":"100",
                "variant_pair_id":"4:100:C:T","effect_allele":"T",
                "other_allele":"C","beta":"-0.2","se":"0.05",
                "p":"0.0001","eaf":"0.8","reported_n":"1296908"}

    def test_beta_flip(self):
        x=audit.harmonize_record(self.row(),self.target())
        self.assertAlmostEqual(x["eur_beta_ALT"],0.2)
        self.assertAlmostEqual(x["eur_eaf_ALT"],0.2)
        self.assertEqual(x["eur_source_pair_id_order"],"LEXICOGRAPHIC_ALLELE_PAIR")

    def test_lexicographic_pair_even_if_effect_first(self):
        r=self.row()
        r["effect_allele"]="C";r["other_allele"]="T"
        r["variant_pair_id"]="4:100:C:T"
        x=audit.harmonize_record(r,self.target())
        self.assertAlmostEqual(x["eur_beta_ALT"],-0.2)

    def test_ultrarare_frequency(self):
        r=self.row()
        r["eaf"]="0.9998"
        self.assertEqual(audit.harmonize_record(r,self.target())["eur_status"],
                         "EUR_ULTRA_RARE_MAF_LT_0_001")

    def test_conflicting_allele_pair_fails(self):
        r=self.row();r["variant_pair_id"]="4:100:A:T"
        with self.assertRaises(ValueError):
            audit.harmonize_record(r,self.target())

    def test_position_scan(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"test.tsv.gz"
            h=["accession","ancestry","phenotype","build","chr","pos",
                "variant_pair_id","effect_allele","other_allele","beta",
                "se","p","eaf","reported_n"]
            r=self.row()
            with gzip.open(p,"wt") as f:
                f.write("\t".join(h)+"\n")
                f.write("\t".join(r[k] for k in h)+"\n")
                r["pos"]="101";r["variant_pair_id"]="4:101:C:T"
                f.write("\t".join(r[k] for k in h)+"\n")
            found=audit.scan_eur(p,[self.target()])
            self.assertEqual(len(found[("4",100)]),1)
            self.assertEqual(len(found),1)

if __name__=="__main__":
    unittest.main()
