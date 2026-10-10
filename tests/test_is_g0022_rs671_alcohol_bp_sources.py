"""Public Zenodo alcohol exposure source download and source/mediation gates."""
import argparse,hashlib,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from download_g0022_japanese_alcohol_unstratified import planned,run as downloader,FILES
from audit_g0022_rs671_alcohol_bp_readiness import build as audit
class AlcoholSourceTests(unittest.TestCase):
    def fixture(self,root,with_complete=False):
        data1=b"unit-test-A";data2=b"unit-test-B"
        metas=[]
        for key,payload in zip(FILES.values(),(data1,data2)):
            metas.append({"key":key,"size":len(payload),
                "checksum":"md5:"+hashlib.md5(payload).hexdigest(),
                "links":{"self":"https://zenodo.org/api/records/10038152/files/"+key+"/content"}})
        root.mkdir(parents=True,exist_ok=True)
        (root/"ZENODO_10038152_METADATA.json").write_text(json.dumps({"id":10038152,"files":metas}))
        if with_complete:
            (root/"1_Alcohol_intake_Unstratified.tsv.gz").write_bytes(data1)
        else:
            (root/"1_Alcohol_intake_Unstratified.tsv.gz.part").write_bytes(data1[:3])
        return json.loads((root/"ZENODO_10038152_METADATA.json").read_text())
    def test_default_plan_only_does_not_download(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);m=self.fixture(root)
            result=downloader(argparse.Namespace(root=root,
              meta=root/"ZENODO_10038152_METADATA.json",
              file="both",execute=False))
            self.assertEqual(len(result),2)
            self.assertEqual(result[0]["status"],"PLAN_ONLY_NOT_VERIFIED")
            self.assertFalse((root/"5_Drinking_Unstratified.tsv.gz").exists())
    def test_allowlist_and_md5_integrity(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);m=self.fixture(root)
            m["files"][0]["links"]["self"]="https://example.org/fake"
            with self.assertRaisesRegex(ValueError,"Non-official"):
                planned(m,root,["alcohol"])
    def test_readiness_reports_partial_download_as_unverified(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)/"base";base.mkdir()
            exposure=Path(td)/"exposure";self.fixture(exposure)
            (base/"G0022_RS671_BBJ_BP_PHEWEB_SOURCE_AUDIT.json").write_text(json.dumps({
               "BP_traits":4,"variant":"12:112241766:G:A","mmHg_trait_units_assumed":False}))
            (base/"G0022_RS671_STROKE_SUBTYPE_SUMMARY.json").write_text(json.dumps({
               "rows":24,"rs671_phenotypes":6}))
            (base/"G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json").write_text(json.dumps({
               "molecular_colocalization_complete":False}))
            result=audit(base,exposure)
            self.assertEqual(result["alcohol_metaGWAS_full_source_md5_verified"],0)
            self.assertFalse(result["eligibility_for_ALDH2_to_alcohol_to_BP_to_stroke_MVMR"])
            self.assertFalse(result["eligible_one_SNP_MR_Egger"])
    def test_a_complete_file_is_still_unvalidated_allele_effect(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)/"base";base.mkdir()
            exposure=Path(td)/"exposure";self.fixture(exposure,True)
            for name,doc in (
              ("G0022_RS671_BBJ_BP_PHEWEB_SOURCE_AUDIT.json",
               {"BP_traits":4,"variant":"12:112241766:G:A","mmHg_trait_units_assumed":False}),
              ("G0022_RS671_STROKE_SUBTYPE_SUMMARY.json",{"rows":24,"rs671_phenotypes":6}),
              ("G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json",{"molecular_colocalization_complete":False})):
                (base/name).write_text(json.dumps(doc))
            result=audit(base,exposure)
            self.assertEqual(result["alcohol_metaGWAS_full_source_md5_verified"],1)
            self.assertEqual(result["alcohol_metaGWAS_rs671_direction_allele_validated"],0)
if __name__=="__main__":unittest.main()
