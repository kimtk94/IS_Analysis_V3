import unittest

from scripts.ckm_stroke import glgc_zero_p_probe as probe

HEADER=("rsID\tCHROM\tPOS_b37\tREF\tALT\tN\tN_studies\tPOOLED_ALT_AF\t"
        "EFFECT_SIZE\tSE\tpvalue_neg_log10\tpvalue\tpvalue_neg_log10_GC\tpvalue_GC")

def row(rsid="rs1",pos=1234,ref="A",alt="G",n="1000",af="0.2",
        beta="1.0",se="0.01",logp="350",p="0"):
    return "\t".join([rsid,"19",str(pos),ref,alt,n,"1",af,
                      beta,se,logp,p,logp,p])

class Reader:
    ROWS=[]
    def __init__(self,url):
        self.url=url
    def __enter__(self):
        return self
    def __exit__(self,*args):
        return False
    def fetch(self,chrom):
        assert chrom=="19"
        return iter(self.ROWS)

def header(url):
    return {"status":"PASS","header":HEADER,"origin":"official_gzip_prefix"}


class TestGLGCOriginalZeroPAudit(unittest.TestCase):
    def setUp(self):
        Reader.ROWS=[]

    def test_reads_original_logp_without_fabricating_p(self):
        Reader.ROWS=[row(),row(pos=1235,p="0.3",logp="0.5228"),
                     row(pos=1236,ref="R",alt="I"),
                     row(pos=1237,ref="A",alt="AT")]
        report,rows=probe.audit_pzero("LDL_C","19",
                                       factory=Reader,header_fetcher=header)
        self.assertEqual(report["status"],"PASS")
        self.assertEqual(report["scientific_status"],"HOLD")
        self.assertEqual(report["source_rows"],4)
        self.assertEqual(report["source_pzero_canonical_snp_count"],1)
        self.assertEqual(report["pzero_with_valid_original_logp"],1)
        self.assertEqual(float(rows[0]["pvalue_neg_log10"]),350)
        self.assertEqual(float(rows[0]["z_abs"]),100)
        self.assertFalse(report["source_full_genome_checksum_verified"])
        self.assertEqual(report["instrument_significance_or_harmonization"],"NOT_RUN")

    def test_invalid_original_neg_logp_forces_hold(self):
        Reader.ROWS=[row(logp="NA")]
        report,rows=probe.audit_pzero("LDL_C","19",
                                       factory=Reader,header_fetcher=header)
        self.assertEqual(report["status"],"HOLD")
        self.assertEqual(report["pzero_with_invalid_auxiliary_metadata"],1)
        self.assertEqual(len(rows),1)

    def test_duplicate_key_with_swapped_effect_alleles_forces_hold(self):
        Reader.ROWS=[row(ref="A",alt="G"),row(rsid="rsdup",ref="G",alt="A")]
        report,rows=probe.audit_pzero("LDL_C","19",
                                       factory=Reader,header_fetcher=header)
        self.assertEqual(report["status"],"HOLD")
        self.assertEqual(report["duplicate_zero_p_variant_keys"],1)

    def test_missing_source_header_is_not_success(self):
        report,rows=probe.audit_pzero("LDL_C","19",factory=Reader,
                        header_fetcher=lambda url:{"status":"HOLD"})
        self.assertEqual(report["status"],"HOLD")
        self.assertEqual(rows,[])

    def test_source_row_cap_is_fail_closed(self):
        Reader.ROWS=[row(pos=1000+k) for k in range(102)]
        report,rows=probe.audit_pzero("LDL_C","19",factory=Reader,
                                      header_fetcher=header,max_rows=100)
        self.assertEqual(report["status"],"HOLD")
        self.assertTrue(report["row_cap_hit"])
        self.assertEqual(report["source_rows"],101)

if __name__=="__main__":
    unittest.main()
