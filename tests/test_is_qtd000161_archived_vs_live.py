import importlib.util
from pathlib import Path
import tempfile
import unittest

FILE=Path(__file__).resolve().parents[1]/"scripts/is/audit_qtd000161_archived_20261011_vs_live.py"
SPEC=importlib.util.spec_from_file_location("qtd_historic_vs_live",FILE)
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

class MockTabix:
    def __init__(self, archive):
        self.data=archive
        self.contigs=["4","10"]
        self.header=[]
    def fetch(self,chrom,start,end):
        # Two official predeclared chromosome-specific windows; emulates
        # exact source rows but input derives from frozen archived source.
        for k,v in self.data.items():
            if k[1]==chrom and start<k[2]<=end:
                gene,c,pos,ref,alt=k
                beta,se,p=v
                yield "\t".join((gene,c,str(pos),ref,alt,
                                str(p),str(beta),str(se)))
    def close(self):pass

class FakePysam:
    def __init__(self,data):self.data=data
    def TabixFile(self,url):return MockTabix(self.data)

HEADER=["molecular_trait_id","chromosome","position","ref","alt","pvalue","beta","se"]
class ArchiveVsLiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.archive=m.load_archived()
    def test_actual_frozen_archive_count_and_sha(self):
        self.assertEqual(len(self.archive),20689)
        self.assertEqual(sum(k[1]=="4" for k in self.archive),11001)
        self.assertEqual(sum(k[1]=="10" for k in self.archive),9688)
    def test_all_identical_source_statistics_pass(self):
        report=m.compare_archived_with_live(self.archive,dict(self.archive))
        self.assertEqual(report["numerically_concordant_beta_se_p"],20689)
        self.assertEqual(report["archived_to_current_QTD_stats_gate"],
          "20689_HISTORICAL_RECONSTRUCTED_QTL_VARIANT_STATISTICS_MATCH_LIVE_QTD")
        self.assertEqual(report["scientific_status"],"HOLD")
        self.assertFalse(report["original_2023_EBI_release6_gzip_sha256_authenticated"])
    def test_wrong_beta_sign_fails_closed(self):
        live=dict(self.archive)
        k=next(k for k,v in live.items() if abs(v[0])>0.01)
        beta,se,p=live[k];live[k]=(-beta,se,p)
        result=m.compare_archived_with_live(self.archive,live)
        self.assertEqual(result["difference_classes"]["beta_discordance"],1)
        self.assertEqual(result["archived_to_current_QTD_stats_gate"],
         "HOLD_HISTORICAL_QTD_STATISTICS_CHANGED_OR_COVERAGE_INCOMPLETE")
    def test_changed_se_fails(self):
        live=dict(self.archive)
        k=next(iter(live))
        beta,se,p=live[k];live[k]=(beta,se*1.1,p)
        result=m.compare_archived_with_live(self.archive,live)
        self.assertEqual(result["difference_classes"]["se_discordance"],1)
    def test_changed_p_fails(self):
        live=dict(self.archive)
        k=next(k for k,v in live.items() if 0<v[2]<.8)
        b,se,p=live[k];live[k]=(b,se,p*1.1)
        result=m.compare_archived_with_live(self.archive,live)
        self.assertEqual(result["difference_classes"]["p_discordance"],1)
    def test_missing_archive_key_does_not_pass(self):
        live=dict(self.archive)
        k=next(iter(live));live.pop(k)
        report=m.compare_archived_with_live(self.archive,live)
        self.assertEqual(report["difference_classes"]["archive_only"],1)
        self.assertIn("HOLD",report["archived_to_current_QTD_stats_gate"])
    def test_modified_alleles_do_not_match(self):
        live=dict(self.archive)
        key=next(iter(live))
        val=live.pop(key)
        opposite=(key[0],key[1],key[2],key[4],key[3])
        live[opposite]=val
        rep=m.compare_archived_with_live(self.archive,live)
        self.assertEqual(rep["difference_classes"]["archive_only"],1)
        self.assertEqual(rep["difference_classes"]["live_only"],1)
    def test_source_mock_full_run_still_holds_coloc(self):
        report=m.live_audit(
            backend=FakePysam(self.archive),
            header_provider=lambda _:{"header":HEADER},
            wait=0)
        self.assertEqual(report["archive_to_live_gate"],
            "20689_HISTORICAL_RECONSTRUCTED_QTL_VARIANT_STATISTICS_MATCH_LIVE_QTD")
        self.assertFalse(report["historic_raw_EBI_QTD_gzip_hash_known"])
        self.assertFalse(report["coloc_performed"])
    def test_archived_panel_modified_bytes_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            import shutil
            for chrom,p in m.ARCHIVE.items():
                path=Path(tmp)/p
                path.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(m.ROOT/p,path)
            path=Path(tmp)/m.ARCHIVE["4"]
            with path.open("at") as f:f.write("corrupt\n")
            with self.assertRaises(ValueError):
                m.load_archived(tmp)
    def test_bounded_24_examples_only(self):
        report=m.compare_archived_with_live(self.archive,{})
        self.assertEqual(report["difference_classes"]["archive_only"],20689)
        self.assertEqual(len(report["bounded_difference_examples"]),24)
if __name__=="__main__":
    unittest.main()
