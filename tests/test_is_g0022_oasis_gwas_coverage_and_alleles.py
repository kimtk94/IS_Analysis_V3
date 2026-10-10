import sys
import csv
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_oasis_gwas_coverage_and_alleles import CS, allele_set_result, audit_page, split_variant, load_ais_four

A = {x[1]: {"gwas_beta_alt":-.149, "gwas_se_alt":.0174, "gwas_alt_eaf":.2319} for x in CS}

class TestOasisAlleleCoverage(unittest.TestCase):
    def test_match(self):
        self.assertEqual(allele_set_result(CS[0][1],CS[0][2],CS[0][2]),"EXACT_REF_ALT_STRING_MATCH_AFTER_PRIOR_LIFTOVER")

    def test_swapped(self):
        self.assertEqual(allele_set_result(CS[0][1],"chr12:111730205:A:G",CS[0][2]),"SWAPPED_REF_ALT_REQUIRES_SIGN_FLIP_IF_ALT_EFFECT_VERIFIED")

    def test_mismatch(self):
        self.assertEqual(allele_set_result(CS[0][1],"chr12:111730205:G:C",CS[0][2]),"ALLELE_MISMATCH_DO_NOT_HARMONIZE")

    def test_invalid_build(self):
        self.assertEqual(allele_set_result(CS[0][1],"chr12:111730206:G:A",CS[0][2]),"MISMATCH_BUILD_OR_POSITION")

    def test_missing_browser_not_true_negative(self):
        info, rows=audit_page("ALDH2","Mono-L1",b"<html>No bokeh</html>",A)
        self.assertEqual(info["source_identified_cs_markers"],0)
        self.assertEqual(info["unobserved_bokeh_markers"],0)
        self.assertTrue(all(r["source_coverage_status"]=="NOT_ASSESSED_NO_SOURCE" for r in rows))

    def test_invalid_variant_ref(self):
        with self.assertRaises(ValueError):split_variant("12:112168009:A:A","GRCh37")

class TestPreHarmonizedAIS(unittest.TestCase):
    def test_gwas_original_ref_succeeds_only_when_oriented_alt(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"effects.tsv"
            fields=["study","trait","variant","ref","effect_allele","original_effect_allele",
                    "ALT_beta","se","ALT_EAF","p"]
            rows=[]
            for rsid,var,gr38,pip in CS:
                ref,alt=var.split(":")[2:]
                original=ref if rsid=="rs77768175" else alt
                rows.append(dict(study="GCST90104545",trait="EAS_AIS",variant=var,
                                 ref=ref,effect_allele="ALT",original_effect_allele=original,
                                 ALT_beta="-0.149",se="0.0174",ALT_EAF="0.23",p="9e-18"))
            with p.open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
                writer.writeheader();writer.writerows(rows)
            self.assertEqual(len(load_ais_four(p)),4)
            rows[-1]["original_effect_allele"]="C"
            with p.open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
                writer.writeheader();writer.writerows(rows)
            with self.assertRaisesRegex(ValueError,"ALT alignment"):
                load_ais_four(p)

if __name__=="__main__":unittest.main()
