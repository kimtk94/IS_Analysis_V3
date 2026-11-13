#!/usr/bin/env python3
"""Unit-test source eQTL nominal row structure using synthetic compressed data."""
import gzip
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_ebi_eqtl_nominal_real_rows import REQUIRED,parse_nominal_rows


def fixture(pvalues=("0.9","0.001"), alt="G", se="0.2"):
    header="\t".join(REQUIRED)
    rows=[]
    for index,p in enumerate(pvalues):
        fields={
            "variant":f"chr1_{100+index}_A_{alt}",
            "pvalue":p,"molecular_trait_id":"ENSG00000111111",
            "maf":"0.23","gene_id":"ENSG00000111111",
            "beta":"0.3","se":se,"chromosome":"1",
            "position":str(100+index),"ref":"A","alt":alt,"an":"200","ac":"46"
        }
        rows.append("\t".join(fields[k] for k in REQUIRED))
    return gzip.compress((header+"\n"+"\n".join(rows)+"\n").encode())


class EbiNominalSourceRowsTests(unittest.TestCase):
    def test_complete_nonsignificant_original_nominal_like(self):
        r=parse_nominal_rows(fixture())
        self.assertEqual(r["assessed_complete_SNP_rows"],2)
        self.assertTrue(r["has_nonsignificant_association_in_prefix"])
        self.assertTrue(r["all_checked_rows_have_beta_se_maf_and_ref_alt"])
        self.assertEqual(r["full_cis_distribution_tested"],"NOT_VERIFIED")

    def test_beta_se_missing_must_fail(self):
        with self.assertRaisesRegex(ValueError,"missing|Missing|numeric|statistical"):
            parse_nominal_rows(fixture(se="NA"))

    def test_invalid_variant_key_fails(self):
        with self.assertRaisesRegex(ValueError,"chromosome"):
            parse_nominal_rows(fixture(alt="T").replace(b"aa",b"bb")) if False else self._malformed_variant()

    def _malformed_variant(self):
        cols=list(REQUIRED)
        values={k:("chr2_100_A_G" if k=="variant" else
                   "1" if k=="chromosome" else
                   "100" if k=="position" else
                   "A" if k=="ref" else
                   "G" if k=="alt" else
                   "ENSG00000111111" if k in ("gene_id","molecular_trait_id") else
                   "0.3" if k=="beta" else
                   "0.1" if k=="se" else
                   "0.4" if k=="maf" else
                   "0.02" if k=="pvalue" else
                   "200" if k=="an" else
                   "80" if k=="ac" else "")
                for k in cols}
        text="\t".join(cols)+"\n"+"\t".join(values[k] for k in cols)+"\n"
        parse_nominal_rows(gzip.compress(text.encode()))

    def test_significant_only_rows_do_not_imply_all_cis(self):
        r=parse_nominal_rows(fixture(pvalues=("0.01","0.03")))
        self.assertFalse(r["has_nonsignificant_association_in_prefix"])
        self.assertEqual(r["full_cis_distribution_tested"],"NOT_VERIFIED")

    def test_corrupted_compressed_input_fails(self):
        with self.assertRaises(ValueError):
            parse_nominal_rows(b"x"*700)

    def test_bad_maf_fails(self):
        cols=list(REQUIRED)
        header="\t".join(cols)
        base={k:("0.7" if k=="maf" else
                   "chr1_100_A_G" if k=="variant" else
                   "1" if k=="chromosome" else
                   "100" if k=="position" else
                   "A" if k=="ref" else "G" if k=="alt" else
                   "ENSG00000111111" if k in ("gene_id","molecular_trait_id") else
                   "0.3" if k=="beta" else
                   "0.1" if k=="se" else
                   "0.01" if k=="pvalue" else
                   "200" if k=="an" else "80" if k=="ac" else "") for k in cols}
        with self.assertRaisesRegex(ValueError,"statistical"):
            parse_nominal_rows(gzip.compress((header+"\n"+"\t".join(base[k] for k in cols)+"\n").encode()))


if __name__=="__main__":
    unittest.main()
