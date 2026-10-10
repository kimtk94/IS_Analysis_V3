"""Offline coordinate, allele, AF and absence regression tests for OASIS full-region join."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_oasis_regional_gwas_overlap import (
    best_primary_chr12_chain, make_mapper, map_ais_variants, compare_one_page
)


class TestOasisRegionalAudit(unittest.TestCase):
    def build_chain(self,tmp):
        # 0-based GRCh37 target 999..1009 -> GRCh38 query 1999..2009;
        # 10bp gap on both source and destination, next block 1019..1029.
        content=(
            "chain 100 chr12 20000 + 999 1029 chr12 20000 + 1999 2029 101\n"
            "10 10 10\n10\n\n"
            "chain 50 chr12 20000 + 999 1009 chr12 20000 + 4000 4010 102\n"
            "10\n\n"
        )
        p=Path(tmp)/"chain.gz"
        with gzip.open(p,"wt") as f: f.write(content)
        return p

    def test_forward_chain_off_by_one_gap_and_choice(self):
        with tempfile.TemporaryDirectory() as tmp:
            c=best_primary_chr12_chain(self.build_chain(tmp))
            self.assertEqual(c["id"],"101")
            lift=make_mapper(c)
            self.assertEqual(lift(1000),2000)
            self.assertEqual(lift(1009),2009)
            self.assertIsNone(lift(1010))
            self.assertIsNone(lift(1019))
            self.assertEqual(lift(1020),2020)
            self.assertEqual(lift(1029),2029)

    def test_source_liftover_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            lift=make_mapper(best_primary_chr12_chain(self.build_chain(tmp)))
            rows={"12:1000:A:G":{"pos_grch37":1000,"ref_grch37":"A","alt_grch37":"G"},
                  "12:1010:C:T":{"pos_grch37":1010,"ref_grch37":"C","alt_grch37":"T"}}
            mapped,counts=map_ais_variants(rows,lift)
            self.assertEqual(len(mapped),1)
            self.assertEqual(counts["NO_PRIMARY_CHAIN_MAPPING"],1)

    def q(self, af=0.8, allele="G", position=2000, rsid="rs1"):
        return dict(gene_name="ALDH2",variant_id_hg38=f"chr12:{position}:A:{allele}",
                    rsid=rsid,pval_nominal=0.30,effect_size=-0.4,
                    effect_size_SE=0.2,maf=af)

    def mapped(self):
        return {2000:dict(gwas_variant_grch37="12:1000:A:G",pos_grch37=1000,
                          ref_grch37="A",alt_grch37="G",gwas_ALT_beta=0.4,
                          gwas_se=0.2,gwas_p=0.02,gwas_ALT_EAF=0.13)}

    def test_bokeh_frequency_above_half_is_valid_af_not_maf(self):
        with patch("audit_g0022_oasis_regional_gwas_overlap.parse",return_value=[self.q()]):
            matched,overview=compare_one_page("ALDH2","Mono-L1",b"mock",self.mapped())
        self.assertEqual(len(matched),1)
        self.assertEqual(overview["bokeh_reported_maf_gt_0_5_count"],1)
        self.assertFalse(overview["bokeh_maf_column_means_true_minor_frequency"])
        self.assertEqual(matched[0]["oasis_browser_maf_column_AS_REPORTED_NOT_TRUE_MAF"],0.8)
        self.assertFalse(matched[0]["qtl_beta_alt_orientation_verified"])
        self.assertFalse(matched[0]["valid_for_full_locus_coloc"])

    def test_position_allele_mismatch_is_never_harmonized(self):
        with patch("audit_g0022_oasis_regional_gwas_overlap.parse",return_value=[self.q(allele="T")]):
            matched,overview=compare_one_page("ALDH2","Mono-L1",b"mock",self.mapped())
        self.assertEqual(len(matched),0)
        self.assertEqual(overview["plot_overlap_types"]["QTL_AND_GWAS_ALLELE_MISMATCH_REJECTED"],1)

    def test_non_source_page_is_not_biological_negative(self):
        with patch("audit_g0022_oasis_regional_gwas_overlap.parse",side_effect=ValueError("Missing embedded Bokeh data")):
            matched,overview=compare_one_page("RPH3A","B_Activated-L2",b"mock",self.mapped())
        self.assertEqual(matched,[])
        self.assertEqual(overview["source_state"],"NOT_ASSESSED_NO_EMBEDDED_BOKEH_ASSOCIATIONS")

    def test_duplicate_qtl_id_is_rejected(self):
        with patch("audit_g0022_oasis_regional_gwas_overlap.parse",
                   return_value=[self.q(rsid="rs1"),self.q(rsid="rs1")]):
            with self.assertRaisesRegex(ValueError,"Duplicated QTL variant"):
                compare_one_page("ALDH2","Mono-L1",b"mock",self.mapped())

    def test_invalid_af_rejected(self):
        with patch("audit_g0022_oasis_regional_gwas_overlap.parse",return_value=[self.q(af=1.1)]):
            with self.assertRaisesRegex(ValueError,"frequency"):
                compare_one_page("ALDH2","Mono-L1",b"mock",self.mapped())


if __name__=="__main__":
    unittest.main()
