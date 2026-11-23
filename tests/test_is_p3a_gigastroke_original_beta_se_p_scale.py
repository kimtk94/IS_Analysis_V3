"""Fail-closed source and statistical model tests for original GIGASTROKE EAS/EUR."""
import gzip
import importlib
import io
import json
import math
import struct
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_gigastroke_original_beta_se_p_scale")
COLS=(
"chromosome","base_pair_location","effect_allele","other_allele",
"effect_allele_frequency","beta","standard_error","p_value",
"rsid"
)


def manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def gw_source(rows=8,beta="0.05",se="0.025",p=None,eaf="0.35",
              ea="A",other="G",duplicates=False):
    if p is None:
        try:p=math.erfc(abs(float(beta)/float(se))/math.sqrt(2))
        except (ValueError, ZeroDivisionError, OverflowError):p=0.5
    header="\t".join(COLS)
    src=[header]
    for i in range(rows):
        r={"chromosome":"10","base_pair_location":str(103000000+i),
           "effect_allele":ea,"other_allele":other,
           "effect_allele_frequency":eaf,
           "beta":beta,"standard_error":se,"p_value":str(p),
           "rsid":f"rs{i}"}
        src.append("\t".join(r[k] for k in COLS))
    return gzip.compress(("\n".join(src)+"\n").encode())


def gwas_with_reported_or(rows=100,reported_OR="valid",include_ci=True):
    original=gzip.decompress(gw_source(rows=rows)).decode("utf-8").splitlines()
    effect=0.05
    se=0.025
    OR=(math.exp(effect) if reported_OR=="valid" else
        None if reported_OR=="missing" else
        1.5 if reported_OR=="wrong" else -0.3)
    source=[original[0]+"\todds_ratio\tci_lower\tci_upper"]
    for row in original[1:]:
        value="" if OR is None else str(OR)
        lo=(str(math.exp(effect-1.96*se)) if include_ci else "")
        hi=(str(math.exp(effect+1.96*se)) if include_ci else "")
        source.append(row+"\t"+value+"\t"+lo+"\t"+hi)
    return gzip.compress(("\n".join(source)+"\n").encode())


def original_paper():
    return (b'<article><front><article-meta><article-id pub-id-type="doi">'
            b'10.1038/s41586-022-05165-3</article-id></article-meta></front>'
            b'<body><p>Individual studies performed a GWAS using logistic regression '
            b'(or Cox regression in some longitudinal population-based cohorts).</p>'
            b'<p>We removed extreme effect sizes with log[OR] > 5, and with '
            b'log[OR] &lt; -5.</p></body></article>')


class FakeHTTP:
    def __init__(self,body,accession,status=206,range_header=None,actual_url=None):
        self.body=body
        self.status=status
        self.accession=accession
        self.headers={"Content-Range":range_header or
            f"bytes 0-{M.MAX_PREFIX_BYTES-1}/{M.PINNED.SOURCES[accession]['bytes']}"}
        self.url=actual_url or M.PINNED.build_urls(
            M.PINNED.SOURCES[accession])[0]
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def geturl(self):return self.url
    def read(self,n):return self.body[:n]


class OriginalPublisherGWASEffectScaleTests(unittest.TestCase):
    def test_real_source_header_numeric_and_z_p_coherence(self):
        r=M.parse_source_effect_stats(gw_source(),"EAS")
        self.assertEqual(r["original_source_preview_rows"],8)
        self.assertEqual(r["z_p_abs_log10_error_gt_1"],0)
        self.assertAlmostEqual(r["z_p_log10_error_median"],0,places=8)
        self.assertFalse(r["source_file_beta_exact_log_OR_scale_verified"])

    def test_original_odds_ratio_populated_validates_lnOR_in_preview(self):
        qc=M.parse_source_effect_stats(gwas_with_reported_or(rows=100),"EAS")
        self.assertTrue(qc["source_OR_and_CI_columns_present"])
        self.assertEqual(qc["source_OR_populated_rows"],100)
        self.assertEqual(qc["source_OR_95CI_populated_rows"],100)
        self.assertEqual(qc["source_beta_ln_OR_max_abs_error_le_0_015_rows"],100)
        self.assertEqual(qc["source_OR_95CI_exp_beta_pm_1_96SE_le_0_04_rows"],100)
        self.assertEqual(qc["source_direct_ln_OR_vs_beta_preview_status"],
            "PASS_ORIGINAL_POPULATED_OR_AGREES_WITH_BETA_IN_PREVIEW")
        self.assertFalse(qc["source_file_beta_exact_log_OR_scale_verified"])

    def test_original_OR_header_but_all_values_missing_keeps_hold(self):
        qc=M.parse_source_effect_stats(gwas_with_reported_or(
            rows=100,reported_OR="missing",include_ci=False),"EUR")
        self.assertEqual(qc["source_OR_populated_rows"],0)
        self.assertEqual(qc["source_OR_missing_or_unparseable_rows"],100)
        self.assertEqual(qc["source_direct_ln_OR_vs_beta_preview_status"],
            "HOLD_OR_SOURCE_NOT_POPULATED_OR_INCONSISTENT")

    def test_original_OR_disagreement_prevents_scale_approval(self):
        qc=M.parse_source_effect_stats(gwas_with_reported_or(
            rows=100,reported_OR="wrong"),"EAS")
        self.assertEqual(qc["source_OR_populated_rows"],100)
        self.assertEqual(qc["source_beta_ln_OR_max_abs_error_le_0_015_rows"],0)
        self.assertEqual(qc["source_direct_ln_OR_vs_beta_preview_status"],
            "HOLD_OR_SOURCE_NOT_POPULATED_OR_INCONSISTENT")

    def test_original_negative_OR_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_source_effect_stats(gwas_with_reported_or(
                rows=100,reported_OR="negative"),"EUR")

    def test_same_numeric_p_not_proof_of_effect_size_unit(self):
        r=M.parse_source_effect_stats(gw_source(beta="-0.03",se="0.02"),"EUR")
        self.assertEqual(r["n_real_effects_negative"],8)
        self.assertTrue(r["source_file_beta_numeric_consistency_is_not_log_OR_unit_proof"])

    def test_p_conflicts_detected_not_fabricated(self):
        r=M.parse_source_effect_stats(gw_source(p=0.5),"EAS")
        self.assertGreater(r["z_p_abs_log10_error_gt_1"],0)

    def test_negative_beta_is_valid_effect_not_or(self):
        r=M.parse_source_effect_stats(gw_source(beta="-0.08",se="0.02"),"EAS")
        self.assertEqual(r["n_real_effects_negative"],8)
        self.assertFalse(r["source_file_beta_exact_log_OR_scale_verified"])

    def test_invalid_effect_inputs_block(self):
        for kwargs in ({"beta":"NaN"},{"se":"0"},{"se":"NaN"},
                       {"eaf":"-0.2"},{"eaf":"NaN"},{"ea":"A","other":"A"},
                       {"p":"2"},{"p":"-0.1"}):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                M.parse_source_effect_stats(gw_source(**kwargs),"EAS")

    def test_missing_header_field_block(self):
        raw=gzip.decompress(gw_source()).decode().replace("standard_error","sd")
        with self.assertRaises(ValueError):
            M.parse_source_effect_stats(gzip.compress(raw.encode()),"EUR")

    def test_malformed_row_block(self):
        raw=gzip.decompress(gw_source()).decode().splitlines()
        raw[2]="\t".join(raw[2].split("\t")[:-2])
        with self.assertRaises(ValueError):
            M.parse_source_effect_stats(gzip.compress(("\n".join(raw)+"\n").encode()),"EAS")

    def test_wrong_accession_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_source_effect_stats(gw_source(),"AFR")

    def test_first_prefix_only_range_206(self):
        spec="EAS"
        body=b"\x1f\x8b"+b"0"*(M.MAX_PREFIX_BYTES-2)
        def opener(req,timeout):
            self.assertEqual(req.headers.get("Range"),"bytes=0-262143")
            return FakeHTTP(body,spec)
        raw,p=M.source_small_gwas_prefix(spec,opener)
        self.assertEqual(len(raw),M.MAX_PREFIX_BYTES)
        self.assertFalse(p["original_full_file_downloaded_in_this_run"])
        self.assertEqual(p["original_full_compressed_file_bytes_pinned"],
                         260_689_208)

    def test_source_prefix_wrong_status_rejected(self):
        body=b"\x1f\x8b"+b"0"*(M.MAX_PREFIX_BYTES-2)
        with self.assertRaises(ValueError):
            M.source_small_gwas_prefix("EUR",lambda req,timeout:
                FakeHTTP(body,"EUR",status=200))

    def test_source_prefix_invalid_release_length_rejected(self):
        body=b"\x1f\x8b"+b"0"*(M.MAX_PREFIX_BYTES-2)
        with self.assertRaises(ValueError):
            M.source_small_gwas_prefix("EAS",lambda req,timeout:
                FakeHTTP(body,"EAS",
                    range_header="bytes 0-262143/260689209"))

    def test_source_prefix_unapproved_redirect_rejected(self):
        body=b"\x1f\x8b"+b"0"*(M.MAX_PREFIX_BYTES-2)
        with self.assertRaises(ValueError):
            M.source_small_gwas_prefix("EAS",lambda req,timeout:
                FakeHTTP(body,"EAS",actual_url="https://evil.example/gwas.gz"))

    def test_primary_paper_log_or_and_mixed_study_design(self):
        d=M.primary_methods(original_paper())
        self.assertTrue(d["original_primary_paper_method_has_logistic_regression"])
        self.assertTrue(d["original_primary_paper_method_has_Cox_regression"])
        self.assertTrue(d["original_primary_paper_variant_quality_control_mentions_log_OR"])
        self.assertTrue(d["cannot_determine_exact_uploaded_EAS_EUR_file_beta_scale_from_paper_alone"])

    def test_wrong_original_study_doi_rejected(self):
        with self.assertRaises(ValueError):
            M.primary_methods(original_paper().replace(
                b"10.1038/s41586-022-05165-3",b"10.1038/s41586-022-00000-0"))

    def test_no_cox_studies_claim_is_not_accepted(self):
        with self.assertRaises(ValueError):
            M.primary_methods(original_paper().replace(b"Cox regression",b"Poisson regression"))

    def test_binary_case_cohort_counts_have_reference_not_APPROVED_s(self):
        r=M.PRIMARY_CASES
        self.assertEqual(r["EAS"]["cases"]+r["EAS"]["controls"],256274)
        self.assertEqual(r["EUR"]["cases"]+r["EUR"]["controls"],1296908)
        self.assertTrue(r["EAS"]["BBJ_is_discovery"])

    def test_yaml_no_numeric_phenotype_N_requires_HOLD(self):
        from unittest.mock import patch
        class FakeY:
            status=200
            def __init__(self,u):self.u=u
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def geturl(self):return self.u
            def read(self,n):return (
                "gwas_id: GCST90104545\n"
                "trait_description: Ischemic stroke\n"
                "genome_assembly: GRCh37\n"
                "data_file_name: GCST90104545_buildGRCh37.tsv.gz\n"
                "data_file_md5sum: 904e12e33ec3c2eeefe4dc93bec4995c\n").encode()
        import urllib.request
        from unittest.mock import patch
        with patch("urllib.request.urlopen",side_effect=lambda req,timeout:
                   FakeY(req.full_url)):
            source=M.source_yaml("EAS")
        self.assertFalse(source["publisher_exact_file_case_control_N_complete"])
        self.assertFalse(source["publisher_exact_file_case_fraction_s_directly_verified"])
        self.assertFalse(source["YAML_identity_verified_not_effect_scale_verified"] is False)

    def test_meta_manifest_holds_all_36_after_real_evidence(self):
        base=manifest()
        self.assertEqual(len(base["groups"]),36)
        self.assertEqual(len({x["gene_id"] for x in base["groups"]}),20)
        self.assertEqual(M.PRIMARY_CASES["EAS"]["BBJ_is_discovery"],True)

    def test_outputs_omit_snp_betas_and_donor_ids(self):
        r=M.parse_source_effect_stats(gw_source(),"EAS")
        s=json.dumps(r)
        self.assertNotIn("103000000",s)
        self.assertNotIn("rs0",s)
        self.assertNotIn('"beta":0.05',s)


if __name__=="__main__":
    unittest.main()
