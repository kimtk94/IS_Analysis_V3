import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_v8_vs_qtd_direct_tabix.py"
spec=importlib.util.spec_from_file_location("gtex_direct_compare",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

HEAD=["molecular_trait_id","chromosome","position","ref","alt","pvalue","beta","se"]
def line(trait="ENSG00000138675.10",chrom="4",pos="79759933",ref="G",alt="C",
         p="1e-3",beta="0.3",se="0.06"):
    return "\t".join((trait,chrom,pos,ref,alt,p,beta,se))
def rec(*args,**kwargs):
    return m.select_qtl_rows(HEAD,[line(*args,**kwargs)])

class MockTabix:
    def __init__(self,url,broken=False):
        self.url=url
        self.header=["#"+"\t".join(HEAD)]
        self.contigs=["4","10"]
        self.broken=broken
        self.closed=False
    def fetch(self,contig,start,end):
        if self.broken:raise OSError("missing remote tabix blocks")
        b="0.4" if self.url==m.REPROCESSED else "0.3"
        yield line(chrom=contig,pos=str(start+1),beta=b)
    def close(self):
        self.closed=True

class FakePysam:
    def TabixFile(self,url):
        return MockTabix(url)

class FailPysam:
    def TabixFile(self,url):
        raise OSError("HTTPS range unavailable")

class DirectQTLTests(unittest.TestCase):
    def test_official_md5_of_tissue(self):
        text="a"*32+" Adipose.tsv.gz\n"+m.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n"
        self.assertTrue(m.verify_manifest(text))
    def test_md5_does_not_match_nearby_tissue(self):
        self.assertFalse(m.verify_manifest("c07a1731eb96aafbd7aed131cb97515a Brain_Cerebellum.tsv.gz"))
    def test_duplicate_md5_fails_closed(self):
        line=m.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n"
        self.assertFalse(m.verify_manifest(line+line))
    def test_header_and_parser(self):
        q=rec()
        self.assertEqual(q["status"],"SCREEN_VALID")
        self.assertEqual(q["observed"],1)
        self.assertEqual(len(q["rows"]),1)
        self.assertEqual(m.normalize_header(["#"+ "\t".join(HEAD)]),HEAD)
    def test_conflicting_identical_snp_fails(self):
        q=m.select_qtl_rows(HEAD,[line(beta="0.3"),line(beta="0.7")])
        self.assertEqual(q["status"],"HOLD_CONFLICTING_DUPLICATE_KEY")
    def test_missing_columns_fails_closed(self):
        q=m.select_qtl_rows(["chromosome"],["4"])
        self.assertEqual(q["status"],"HOLD_MISSING_QTL_FIELDS")
    def test_row_count_cap_fails(self):
        q=m.select_qtl_rows(HEAD,[line(),line()],line_limit=1)
        self.assertEqual(q["status"],"HOLD_REGION_ROW_CAP")
    def test_direct_same_variant_different_beta_quantified(self):
        a=rec(beta="0.3",p="1e-3")
        b=rec(beta="0.4",p="1e-4")
        x=m.compare_rows(a,b)
        self.assertEqual(x["matched_same_trait_and_exact_GRCh38_REF_ALT_pairs"],1)
        self.assertEqual(x["beta_sign_matches"],1)
        self.assertEqual(x["beta_numerically_near_equal"],0)
        self.assertAlmostEqual(x["median_abs_beta_difference"],0.1)
        self.assertAlmostEqual(x["median_abs_log10_p_difference"],1)
        self.assertFalse(x["coloc_performed"])
    def test_variant_different_ref_alt_no_match(self):
        a=rec(ref="A",alt="G")
        b=rec(ref="G",alt="A")
        x=m.compare_rows(a,b)
        self.assertEqual(x["matched_same_trait_and_exact_GRCh38_REF_ALT_pairs"],0)
        self.assertEqual(x["comparison_status"],"HOLD_NO_EXACT_GENE_VARIANT_OVERLAP")
    def test_two_remote_bounded_regions_fake(self):
        source=(m.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n").encode()
        result=m.live_audit(FakePysam(),lambda _:source,wait_seconds=0)
        self.assertEqual(result["total_exact_pairs_two_small_regions"],2)
        self.assertEqual(result["source_region_gate"],
                         "SMALL_REGION_GTEXv8_VS_QTD_QTL_NUMERIC_COMPARE_PASS_REVIEW_SOURCE_DIFFERENCES")
        self.assertFalse(result["GTEx_original_complete_cis_universe_authenticated"])
        self.assertEqual(result["scientific_status"],"HOLD")
    def test_unavailable_remote_blocks(self):
        source=(m.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n").encode()
        r=m.live_audit(FailPysam(),lambda _:source,wait_seconds=0)
        self.assertEqual(r["source_region_gate"],"HOLD_REMOTE_TABIX_SOURCE_UNAVAILABLE")
        self.assertFalse(r["coloc_performed"])
if __name__=="__main__":
    unittest.main()
