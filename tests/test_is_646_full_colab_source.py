"""Independent full-646 Colab source-provenance verifier: synthetic fail-closed fixtures."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
MODULE=ROOT/"scripts/is/audit_is_646_full_colab_source.py"
spec=importlib.util.spec_from_file_location("ISFullReplayVerifier",MODULE)
verifier=importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)

def fixture():
    rows,master,index,hashes=[],[],[],[]
    p=(0.01,0.09,0.10,0.60,0.20)
    for i in range(646):
        k=dict(locus="LOC1",dataset_key=f"GTEx_V8__Tissue_{i:04d}",
               gene_base=f"ENSG000000{i:05d}")
        filename="__".join(k.values())+".tsv"
        rows.append(dict(k,filename=filename,status="PASS",actual_nsnps="123",
            expected_nsnps="123",unique_snps="123",
            qtl_n_first_actual="200",qtl_n_first_expected="200",
            original_qtl_n_range="2",observed_qtl_n_range="2",
            max_abs_hypothesis_delta="0",
            **{f"PP_H{n}":str(p[n]) for n in range(5)}))
        master.append(dict(k,status="PASS",nsnps="123",qtl_n="200",
              **{f"PP.H{n}":str(p[n]) for n in range(5)}))
        index.append(dict(k,status="READY",file="/in/"+filename,
                          nsnps="123",qtl_n_first="200",qtl_n_min="198",qtl_n_max="200"))
        hashes.append(dict(k,filename=filename,status="PASS",source_sha256="a"*64))
    return rows,master,index,hashes

class TestISFull646SourceProof(unittest.TestCase):
    def test_synthetic_complete_646_success(self):
        x=fixture()
        result=verifier.validate_646_rows(*x)
        self.assertEqual(result["n_replayed"],646)
        self.assertEqual(result["n_PP_H4_ge_0_5"],0)
        self.assertEqual(result["QTL_n_varied_within_test_count"],646)

    def test_single_failed_status_invalidates_full_release(self):
        x=fixture()
        x[0][100]["status"]="FAIL"
        with self.assertRaisesRegex(ValueError,"NOT_646_REPLAYS_PASS"):
            verifier.validate_646_rows(*x)

    def test_mismatched_snp_cardinality_fails(self):
        x=fixture()
        x[0][200]["unique_snps"]="121"
        with self.assertRaisesRegex(ValueError,"SNP_COUNT_DISAGREE"):
            verifier.validate_646_rows(*x)

    def test_mismatched_posterior_fails(self):
        x=fixture()
        x[0][300]["PP_H4"]="0.3"
        with self.assertRaisesRegex(ValueError,"POSTERIOR_NOT_VALID"):
            verifier.validate_646_rows(*x)

    def test_missing_source_hash_fails(self):
        x=fixture()
        x[3][0]["source_sha256"]=""
        with self.assertRaisesRegex(ValueError,"SOURCE_SHA256_RECORD_INVALID"):
            verifier.validate_646_rows(*x)

    def test_no_output_overwrite(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            out=d/"existing"
            out.mkdir()
            sentinel=out/"old.txt"
            sentinel.write_text("PRESERVE")
            with self.assertRaises(FileExistsError):
                verifier.run(d/"none",d/"idx",d/"master",d/"script",d/"cache",out)
            self.assertEqual(sentinel.read_text(),"PRESERVE")

if __name__=="__main__":
    unittest.main()
