import csv
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

FILE=Path(__file__).resolve().parents[1]/"scripts/is/audit_legacy_gtexv8_original_source_identity.py"
SPEC=importlib.util.spec_from_file_location("qtl_legacy",FILE)
m=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(m)

HEADER=["molecular_trait_id","chromosome","position","ref","alt","pvalue","beta","se","maf","ac","an"]

class MockTabix:
    def __init__(self,archive):
        self.archive=archive
        self.contigs=["4","10"]
        self.header=[]
    def fetch(self,chrom,start,end):
        for key,row in self.archive.items():
            gene,c,pos,ref,alt=key
            if c!=chrom or not(start<pos<=end):continue
            yield "\t".join([gene,c,str(pos),ref,alt,str(row["p"]),
                 str(row["beta"]),"0.05",str(row["maf"]),str(row["ac"]),str(row["an"])])
    def close(self):pass

class MockModule:
    def __init__(self,archive):self.archive=archive
    def TabixFile(self,url):return MockTabix(self.archive)

class SourceIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.archive=m.load_legacy()
        cls.live={k:(v["beta"],0.05,v["p"]) for k,v in cls.archive.items()}
        cls.qc={k:{"maf":v["maf"],"ac":v["ac"],"an":v["an"]} for k,v in cls.archive.items()}
    def test_recovered_historic_files_verified_with_hash(self):
        self.assertEqual(len(self.archive),6242)
        self.assertEqual(sum(key[1]=="4" for key in self.archive),3136)
        self.assertEqual(sum(key[1]=="10" for key in self.archive),3106)
    def test_full_mock_source_with_exact_beta_and_p_is_pass(self):
        o=m.compare(self.archive,self.live,self.qc)
        self.assertEqual(o["archived_source_p_and_beta_match_official_original"],6242)
        self.assertEqual(o["legacy_QTL_statistical_source_identity_gate"],
          "6242_ARCHIVED_GTEX_V8_QTL_STATS_MATCH_EBI_IMPORTED_ORIGINAL")
        self.assertFalse(o["MR_or_coloc_performed"])
    def test_beta_sign_swap_does_not_pass(self):
        l=dict(self.live);k=next(k for k,v in l.items() if abs(v[0])>.01)
        b,se,p=l[k];l[k]=(-b,se,p)
        r=m.compare(self.archive,l,self.qc)
        self.assertEqual(r["difference_classes"]["beta_fail"],1)
        self.assertIn("HOLD",r["legacy_QTL_statistical_source_identity_gate"])
    def test_pvalue_change_does_not_pass(self):
        l=dict(self.live);k=next(k for k,v in l.items() if 0<v[2]<.6)
        b,se,p=l[k];l[k]=(b,se,p*1.15)
        z=m.compare(self.archive,l,self.qc)
        self.assertEqual(z["difference_classes"]["p_fail"],1)
    def test_deleted_variant_holds(self):
        l=dict(self.live);l.pop(next(iter(l)))
        z=m.compare(self.archive,l,self.qc)
        self.assertEqual(z["difference_classes"]["missing_official_original"],1)
        self.assertIn("HOLD",z["legacy_QTL_statistical_source_identity_gate"])
    def test_reversed_ref_alt_does_not_match(self):
        l=dict(self.live);k=next(iter(l));v=l.pop(k)
        l[(k[0],k[1],k[2],k[4],k[3])]=v
        z=m.compare(self.archive,l,self.qc)
        self.assertEqual(z["difference_classes"]["missing_official_original"],1)
    def test_ac_complement_is_not_beta_sign_flip(self):
        q={k:dict(v) for k,v in self.qc.items()}
        k=next(iter(q))
        q[k]["ac"]=q[k]["an"]-q[k]["ac"]
        z=m.compare(self.archive,self.live,q)
        self.assertEqual(z["archived_source_p_and_beta_match_official_original"],6242)
        self.assertTrue(z["no_allele_beta_flipping"])
        self.assertTrue(any(x["ac_complement"] for x in z["by_gene"].values()))
    def test_malformed_frozen_panel_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            for chrom,(path,_) in m.FROZEN.items():
                target=Path(d)/path
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(m.ROOT/path,target)
            with (Path(d)/m.FROZEN["10"][0]).open("a") as f:
                f.write("tamper\n")
            with self.assertRaises(ValueError):m.load_legacy(d)
    def test_full_workflow_mock_has_science_hold(self):
        z=m.run(backend=MockModule(self.archive),
                header_provider=lambda _:{"header":HEADER},wait=0)
        self.assertEqual(z["legacy_source_gate"],
         "6242_ARCHIVED_GTEX_V8_QTL_STATS_MATCH_EBI_IMPORTED_ORIGINAL")
        self.assertFalse(z["formal_MR_or_coloc_performed"])
        self.assertEqual(z["scientific_status"],"HOLD")
    def test_missing_live_source_reports_hold(self):
        class Unavailable:
            def TabixFile(self,url):raise OSError("temporary outage")
        o=m.run(backend=Unavailable(),header_provider=lambda _:{"header":HEADER},wait=0)
        self.assertEqual(o["legacy_source_gate"],"HOLD_SOURCE_NETWORK_OR_SCHEMA_FAILURE")
    def test_report_samples_limited(self):
        z=m.compare(self.archive,{}, {})
        self.assertEqual(len(z["bounded_discordant_examples"]),m.MAX_EXAMPLES)
        self.assertEqual(z["difference_classes"]["missing_official_original"],6242)
if __name__=="__main__":
    unittest.main()
