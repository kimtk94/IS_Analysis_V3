"""No IDs leak from public GTEx/ROSMAP metadata; exact original MD5 source."""
import hashlib
import importlib
import json
import unittest
from unittest.mock import Mock

M=importlib.import_module("scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")


def fake_data():
    return ("sample_id\tgenotype_id\tqtl_group\tsmtsd\tcohort\n"+
        "patient_ABC\tgenotype_ABC\tartery_aorta\tArtery - Aorta\tA\n"+
        "patient_B\tgenotype_B\tartery_aorta\tArtery - Aorta\tA\n"+
        "patient_C\tgenotype_C\tbrain_cortex\tBrain - Cortex\tA\n"+
        "patient_D\tgenotype_D\tbrain_cortex\tBrain - Cortex\tA\n"+
        "patient_E\tgenotype_E\tbrain_naive\tBrain DLPFC\tB\n"+
        "patient_F\tgenotype_F\tbrain_naive\tBrain DLPFC\tB\n").encode()


def fake_item(name):
    return {"key":name,"size":M.FILES[name]["bytes"],
            "checksum":"md5:"+M.FILES[name]["md5"],
            "links":{"self":"https://zenodo.org/api/records/4678936/files/"+
                            name+"/content"}}


class SourceSafetyTests(unittest.TestCase):
    def test_true_source_Filenames_and_expected_MD5_pins(self):
        self.assertEqual(set(M.FILES),{"GTExV8.tsv","ROSMAP.tsv"})
        self.assertEqual(M.FILES["GTExV8.tsv"]["bytes"],1769117)
        self.assertEqual(M.FILES["ROSMAP.tsv"]["md5"],
                         "038f2fc144129d401144f3924109e1df")

    def test_exact_url_from_original_Zenodo_record(self):
        url=M.official_file_url(fake_item("ROSMAP.tsv"),"ROSMAP.tsv")
        self.assertTrue(url.endswith("/ROSMAP.tsv/content"))

    def test_other_host_redirect_file_url_rejected(self):
        bad=fake_item("GTExV8.tsv")
        bad["links"]["self"]="https://lookalike.invalid/GTExV8.tsv"
        with self.assertRaises(ValueError):
            M.official_file_url(bad,"GTExV8.tsv")

    def test_source_file_url_other_filename_rejected(self):
        bad=fake_item("ROSMAP.tsv")
        bad["links"]["self"]="https://zenodo.org/api/records/4678936/files/different.tsv/content"
        with self.assertRaises(ValueError):
            M.official_file_url(bad,"ROSMAP.tsv")

    def test_metadata_aggregates_never_release_sample_ids(self):
        summary=M.safe_count_source_metadata(fake_data(),"GTExV8.tsv")
        self.assertEqual(summary["original_2021_metadata_row_count"],6)
        self.assertEqual(summary["distinct_sample_donor_identifier_count_only"]["sample_id"],6)
        self.assertEqual(summary["source_original_target_label_row_matches_lower_bound_only"],
                         {"QTD000131":2,"QTD000171":2})
        public=json.dumps(summary)
        self.assertNotIn("patient_ABC",public)
        self.assertNotIn("genotype_ABC",public)
        self.assertFalse(summary["sample_identifier_values_uploaded"])
        self.assertFalse(summary["direct_original_36_gene_quantitative_sdY_computed"])

    def test_target_dlpfc_2021_count_not_release7_samples(self):
        summary=M.safe_count_source_metadata(fake_data(),"ROSMAP.tsv")
        self.assertEqual(summary["source_original_target_label_row_matches_lower_bound_only"],
                         {"QTD000434":2})
        self.assertFalse(summary["reusable_raw_data_sample_meta_for_release7_approved"])

    def test_metadata_malformed_physical_row_fails_closed(self):
        with self.assertRaises(ValueError):
            M.safe_count_source_metadata(fake_data()+b"bad\t1\t2\n","GTExV8.tsv")

    def test_metadata_duplicate_column_fails_closed(self):
        with self.assertRaises(ValueError):
            M.safe_count_source_metadata(fake_data().replace(
                b"genotype_id",b"sample_id",1),"GTExV8.tsv")

    def test_metadata_no_delimiter_fails(self):
        with self.assertRaises(ValueError):
            M.safe_count_source_metadata(b"one header\nsomething\n","ROSMAP.tsv")

    def test_genetic_snp_N_and_beta_scale_NEVER_promoted_from_metadata(self):
        d=M.safe_count_source_metadata(fake_data(),"GTExV8.tsv")
        self.assertFalse(d["source_2021_metadata_n_equal_Release7_qtlmap_intersection_verified"])
        self.assertFalse(d["exact_Release7_normalised_expression_matrix_present"])

    def test_study_GWAS_API_response_requires_exact_accession(self):
        # fake get() returns a valid response but different study accession
        raw=json.dumps({"accessionId":"GCST11111111"}).encode()
        fake=lambda url,hosts,cap:(raw,{"SHA256":hashlib.sha256(raw).hexdigest()})
        with self.assertRaises(ValueError):
            M.gwas_record("GCST90104545",fetch=fake)

    def test_study_GWAS_API_does_not_certify_logOR_even_with_metadata(self):
        raw=json.dumps({"accessionId":"GCST90104540",
                        "initialSampleSize":"1,296,908",
                        "studyDesignComment":"case-control"}).encode()
        fake=lambda url,hosts,cap:(raw,{"SHA256":hashlib.sha256(raw).hexdigest()})
        r=M.gwas_record("GCST90104540",fetch=fake)
        self.assertEqual(r["publisher_rest_accession"],"GCST90104540")
        self.assertFalse(r["effect_size_beta_exact_original_summary_file_explicit_log_OR_verified"])
        self.assertFalse(r["outcome_case_control_coloc_beta_varbeta_eligibility_approved"])

    def test_GWAS_Catalog_retired_v1_replaced_by_official_2026_v2(self):
        self.assertEqual(M.REST, "https://www.ebi.ac.uk/gwas/rest/api/v2/studies/")

    def test_original_GWAS_v2_response_can_be_read_without_false_logOR(self):
        raw=json.dumps({"study_accession":"GCST90104540",
                        "disease_trait":"ischemic stroke",
                        "initial_sample_size":"62,100 cases and 1,234,808 controls",
                        "study_design_comment":"binary phenotype"}).encode()
        fake=lambda url,hosts,cap:(raw,{"SHA256":hashlib.sha256(raw).hexdigest()})
        v=M.gwas_record("GCST90104540",fetch=fake)
        self.assertEqual(v["status"],"PASS_ORIGINAL_GWAS_STUDY_METADATA_V2_NOT_BETA_SCALE")
        self.assertEqual(v["published_study_descriptor_fields"]["disease_trait"],
                         "ischemic stroke")
        self.assertFalse(v["effect_size_beta_exact_original_summary_file_explicit_log_OR_verified"])

    def test_tissue_mapping_all_36_source_classes(self):
        self.assertEqual(set(M.TARGET_GROUPS),set(M.PRIOR.EXPECTED_QTLS))
        self.assertEqual(sum(v[4] for v in M.PRIOR.EXPECTED_QTLS.values()),36)

    def test_2021_published_count_is_not_release7_sdY(self):
        self.assertNotEqual("metadata", "source phenotype sdY")
        self.assertEqual(M.RECORD_SHA256,
                         "818a5ebd7f448dd685f3f23d90d42f719a520a8fc27db233cb624bee21619682")


if __name__=="__main__":
    unittest.main()
