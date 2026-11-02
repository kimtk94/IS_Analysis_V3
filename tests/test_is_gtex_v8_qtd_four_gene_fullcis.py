import importlib.util
from pathlib import Path
import unittest

FILE=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_v8_qtd_four_gene_fullcis.py"
SPEC=importlib.util.spec_from_file_location("gtex_fullcis",FILE)
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

HEAD=["molecular_trait_id","chromosome","position","ref","alt","pvalue","beta","se"]
GENES=m.BASE.TARGET_GENES

def make_line(gene,chrom,pos,beta="0.3",se="0.05",p="0.0001",alt="C",ref="A"):
    return "\t".join([GENES[gene]+".8",chrom,str(pos),ref,alt,p,beta,se])

class MockTabix:
    def __init__(self,url,broken=False):
        self.url=url
        self.contigs=("4","10")
        self.header=["#"+"\t".join(HEAD)]
        self.broken=broken
    def fetch(self,c,start,end):
        if self.broken:raise OSError("remote interrupted")
        gene1,gene2=("FGF5","C4orf22") if c=="4" else ("CALHM2","NEURL")
        for i,g in enumerate((gene1,gene2)):
            pos=start+150+i*10
            yield make_line(g,c,pos,beta="0.23" if self.url==m.BASE.ORIGINAL else "0.38")
        # Untargeted neighbors should be excluded but count toward safety cap
        yield "ENSG9999999999.1\t"+c+"\t"+str(start+150)+"\tA\tC\t0.2\t0.1\t0.05"
    def close(self):pass

class Module:
    def TabixFile(self,url):return MockTabix(url)
class DownModule:
    def TabixFile(self,url):return MockTabix(url,True)

def source(rows):
    return {"status":"SCREEN_VALID","rows":rows,"scanned":len(rows),"target_rows":len(rows)}

class FullCisTests(unittest.TestCase):
    def test_only_four_approved_genes_retained(self):
        h=MockTabix(m.BASE.ORIGINAL)
        result=m.read_source_rows(h,HEAD,m.FULL_WINDOWS[0])
        self.assertEqual(result["status"],"SCREEN_VALID")
        self.assertEqual(result["scanned"],3)
        self.assertEqual(result["target_rows"],2)
        self.assertEqual(len(result["rows"]),2)

    def test_non_snp_explicit_indel_and_mnv_included_without_symbolic(self):
        class Mixed(MockTabix):
            def fetch(self,c,start,end):
                yield make_line("FGF5",c,start+50,ref="A",alt="ATG")
                yield make_line("FGF5",c,start+60,ref="GGA",alt="G")
                yield make_line("C4orf22",c,start+70,ref="AC",alt="GT")
                yield make_line("FGF5",c,start+80,ref="A",alt="<DEL>")
        out=m.read_source_rows(Mixed(m.BASE.ORIGINAL),HEAD,m.FULL_WINDOWS[0])
        self.assertEqual(out["status"],"SCREEN_VALID")
        self.assertEqual(len(out["rows"]),3)
        self.assertEqual(out["counts"]["complex_indel_mnv_rows"],3)
        self.assertEqual(out["counts"]["unsupported_allele_encoding_rows"],1)
        self.assertNotIn(("ENSG00000138675","4",79000079,"A","<DEL>"),out["rows"])

    def test_original_and_reprocessed_same_indel_exact_pair(self):
        key=(GENES["FGF5"],"4",80010000,"A","AT")
        a=source({key:(0.3,0.05,1e-5)})
        b=source({key:(-0.1,0.06,1e-3)})
        report=m.summarize_full_cis(a,b)
        self.assertEqual(report["FGF5"]["matched_same_trait_and_exact_GRCh38_REF_ALT_pairs"],1)
        self.assertEqual(report["FGF5"]["beta_sign_matches"],0)

    def test_source_scan_limit_is_fail_closed(self):
        h=MockTabix(m.BASE.ORIGINAL)
        result=m.read_source_rows(h,HEAD,m.FULL_WINDOWS[0],max_rows=1)
        self.assertEqual(result["status"],"HOLD_SCAN_CAP")
        self.assertFalse(result["rows"])

    def test_bad_header_stays_hold(self):
        x=m.read_source_rows(MockTabix(m.BASE.ORIGINAL),["chromosome"],m.FULL_WINDOWS[0])
        self.assertEqual(x["status"],"HOLD_REQUIRED_COLUMNS")

    def test_duplicate_conflicting_values_hold(self):
        class Dup(MockTabix):
            def fetch(self,c,start,end):
                yield make_line("FGF5",c,start+50,beta="0.1")
                yield make_line("FGF5",c,start+50,beta="-0.4")
        x=m.read_source_rows(Dup(m.BASE.ORIGINAL),HEAD,m.FULL_WINDOWS[0])
        self.assertEqual(x["status"],"HOLD_CONFLICTING_DUPLICATE_SOURCE_VARIANT")

    def test_invalid_target_statistics_are_recorded_not_silently_approved(self):
        class Bad(MockTabix):
            def fetch(self,c,start,end):
                yield make_line("FGF5",c,start+50,se="-1")
                yield make_line("FGF5",c,start+60,beta="nan")
        x=m.read_source_rows(Bad(m.BASE.ORIGINAL),HEAD,m.FULL_WINDOWS[0])
        self.assertEqual(x["counts"]["invalid_statistic_rows"],2)
        self.assertEqual(len(x["rows"]),0)

    def test_full_cis_gene_summary_uses_exact_variant_key(self):
        key=("ENSG00000138675","4",80000000,"A","C")
        opposite=("ENSG00000138675","4",80000000,"C","A")
        a=source({key:(0.3,0.05,1e-3)})
        b=source({opposite:(-0.3,0.05,1e-3)})
        out=m.summarize_full_cis(a,b)
        self.assertEqual(out["FGF5"]["matched_same_trait_and_exact_GRCh38_REF_ALT_pairs"],0)
        self.assertEqual(out["FGF5"]["original_only_variant"],1)
        self.assertEqual(out["FGF5"]["reprocessed_only_variant"],1)

    def test_mock_full_runner_produces_evidence_but_not_causal_claim(self):
        md5=(m.BASE.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n").encode()
        out=m.run(Module(),lambda _:{"header":HEAD},lambda _:md5,wait=0)
        self.assertEqual(out["total_exact_gene_variant_pairs"],4)
        self.assertEqual(set(out["four_gene_fullcis_screen"]),set(GENES))
        self.assertEqual(out["fullcis_source_gate"],"HOLD_FULL_CIS_DENOMINATOR_SOURCE_OR_ACCESS")
        self.assertFalse(out["expected_QTD_fullcis_reprocessed_gene_counts_equal"])
        self.assertEqual(out["scientific_status"],"HOLD")
        self.assertFalse(out["coloc_performed"])

    def test_unavailable_remote_ranges_stay_hold(self):
        result=m.run(DownModule(),lambda _:{"header":HEAD},
              lambda _:(m.BASE.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n").encode(),wait=0)
        self.assertEqual(result["fullcis_source_gate"],"HOLD_FULL_CIS_DENOMINATOR_SOURCE_OR_ACCESS")
        self.assertFalse(result["all_source_windows_complete"])

    def test_full_cis_boundaries_follow_previous_documented_source(self):
        self.assertEqual(m.FULL_WINDOWS,[("4",79000000,81500000),
                                         ("10",101500000,106500000)])
        self.assertEqual(m.EXPECTED_REPROCESSED_CIS,
           {"FGF5":5511,"C4orf22":5490,"CALHM2":4868,"NEURL":4820})

    def test_reported_full_screen_only_after_expected_denominators(self):
        save=dict(m.EXPECTED_REPROCESSED_CIS)
        try:
            for g in m.EXPECTED_REPROCESSED_CIS:
                m.EXPECTED_REPROCESSED_CIS[g]=1
            md5=(m.BASE.ORIGINAL_MD5+" Brain_Cerebellar_Hemisphere.tsv.gz\n").encode()
            out=m.run(Module(),lambda _:{"header":HEAD},lambda _:md5,wait=0)
            self.assertEqual(out["fullcis_source_gate"],
                 "FOUR_GENE_FULL_CIS_SOURCE_COMPARISON_SCREEN_PASS_NOT_COLOC")
            self.assertFalse(out["original_full_file_checksum_recomputed"])
            self.assertEqual(out["scientific_status"],"HOLD")
        finally:
            m.EXPECTED_REPROCESSED_CIS.clear()
            m.EXPECTED_REPROCESSED_CIS.update(save)

if __name__=="__main__":
    unittest.main()
