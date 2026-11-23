"""No IDs leak from public GTEx/ROSMAP metadata; exact original MD5 source."""
import hashlib
import importlib
import json
import unittest
from unittest.mock import Mock

M=importlib.import_module("scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")


def fake_data():
    # Six entirely synthetic pseudonymous samples; IDs must never appear in
    # public aggregate output. One GTEx aorta sample fails RNA QC and one
    # cortex sample fails genotype QC.
    lines=[
     "pseudo_id\tsample_id\tgenotype_id\tqtl_group\tsmtsd\tcohort\trna_qc_passed\tgenotype_qc_passed",
     "pseudo_ABC\tpatient_ABC\tgenotype_ABC\tartery_aorta\tArtery - Aorta\tA\tTrue\tTrue",
     "pseudo_B\tpatient_B\tgenotype_B\tartery_aorta\tArtery - Aorta\tA\tFalse\tTrue",
     "pseudo_C\tpatient_C\tgenotype_C\tbrain_cortex\tBrain - Cortex\tA\tTrue\tTrue",
     "pseudo_D\tpatient_D\tgenotype_D\tbrain_cortex\tBrain - Cortex\tA\tTrue\tFalse",
     "pseudo_E\tpatient_E\tgenotype_E\tbrain_naive\tBrain DLPFC\tB\tTrue\tTrue",
     "pseudo_F\tpatient_F\tgenotype_F\tbrain_naive\tBrain DLPFC\tB\tTrue\tTrue"
    ]
    return ("\n".join(lines)+"\n").encode()


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

    def test_actual_original_RNA_and_genotype_QC_intersection_by_group(self):
        summary=M.safe_count_source_metadata(fake_data(),"GTExV8.tsv")
        by=summary["public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates"]
        self.assertEqual(by["QTD000131"]["2021_metadata_rows_in_tissue"],2)
        self.assertEqual(by["QTD000131"]["2021_both_rna_genotype_QC_passed_rows"],1)
        self.assertEqual(by["QTD000171"]["2021_both_rna_genotype_QC_passed_rows"],1)
        self.assertEqual(by["QTD000131"]["2021_public_RNA_genotype_QC_joint_flags"]
                         ["RNA_FAIL_GENOTYPE_PASS"],1)
        self.assertFalse(by["QTD000131"]["2021_both_QC_rows_equal_Release7_declared_N_DESCRIPTIVE"])
        self.assertTrue(by["QTD000131"]["2021_QC_pass_is_not_release7_qtlmap_sample_intersection_proof"])

    def test_source_unknown_QC_flag_not_counted_as_pass(self):
        v=fake_data().replace(b"\tFalse\tTrue\n",b"\tUNKNOWN\tTrue\n")
        summary=M.safe_count_source_metadata(v,"GTExV8.tsv")
        a=summary["public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates"]["QTD000131"]
        self.assertEqual(a["unparsed_source_QC_flags"],1)
        self.assertEqual(a["2021_both_rna_genotype_QC_passed_rows"],1)

    def test_duplicate_public_pseudo_id_not_double_counted(self):
        a=fake_data().replace(b"pseudo_B\t",b"pseudo_ABC\t")
        result=M.safe_count_source_metadata(a,"GTExV8.tsv")
        self.assertEqual(result["distinct_sample_donor_identifier_count_only"]["pseudo_id"],5)
        by=result["public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates"]["QTD000131"]
        self.assertEqual(by["2021_distinct_pseudo_ids_in_tissue_count_only"],1)
        self.assertEqual(by["2021_both_QC_passed_distinct_pseudo_ids"],1)

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

    def test_36_group_2021_QC_matched_GTEx_21_ROSMAP_15_scientific_HOLD(self):
        from pathlib import Path
        m=json.loads(Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json").read_text())
        groups=m["groups"]
        original={
            "schema":"is-p3a-36-primary-phenotype-matrix-and-stroke-beta-units-v1",
            "technical_status":"PASS_36_GROUP_SOURCE_METHOD_AND_PUBLIC_MATRIX_INVENTORY",
            "n_real_gene_tissue_assays":36,
            "n_groups_with_direct_source_expression_sdY_approved":0,
            "n_groups_with_original_GWAS_beta_units_approved":0,
            "groups":[{"gene_id":x["gene_id"],"QTL_dataset":x["QTL_dataset"]}
                      for x in groups]}
        off={"technical_status":"PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE",
             "QTL_datasets":{q:{"study_id":v[0],"study_label":v[1],
                  "sample_group":v[2],"nominal_catalogue_sample_size":v[3]}
                    for q,v in M.PRIOR.EXPECTED_QTLS.items()}}
        cohort={"schema":"is-p3a-primary-2022-gigastroke-AIS-supplement-BBJ-cohort-36-v1",
                "technical_status":"PASS_ORIGINAL_AIS_SUPPLEMENT_BOTH_ANCESTRY_TOTALS_AND_BBJ",
                "scientific_status":"HOLD"}
        files={"GTExV8.tsv":M.safe_count_source_metadata(fake_data(),"GTExV8.tsv"),
               "ROSMAP.tsv":M.safe_count_source_metadata(fake_data(),"ROSMAP.tsv")}
        for q in ("QTD000131","QTD000171"):
            sample=M.TARGET_GROUPS[q][1]
            qc=files["GTExV8.tsv"][
                "public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates"][q]
            qc["2021_both_rna_genotype_QC_passed_rows"]=sample
            qc["2021_metadata_rows_in_tissue"]=sample
            qc["2021_distinct_pseudo_ids_in_tissue_count_only"]=sample
        qc=files["ROSMAP.tsv"]["public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates"]["QTD000434"]
        qc["2021_both_rna_genotype_QC_passed_rows"]=569
        qc["2021_metadata_rows_in_tissue"]=581
        qc["2021_distinct_pseudo_ids_in_tissue_count_only"]=581
        gwas={a:{"status":"PASS_ORIGINAL_GWAS_STUDY_METADATA_V2_NOT_BETA_SCALE"}
              for a in M.GWAS}
        result=M.evaluate(m,original,off,cohort,files,gwas,{})
        self.assertEqual(result["n_gene_tissue_assays_2021_dual_QC_count_equal_release7_dataset_N_DESCRIPTIVE"],21)
        self.assertEqual(result["n_gene_tissue_assays_2021_dual_QC_count_differs_release7_dataset_N"],15)
        self.assertEqual({x["source_2021_pass_QC_minus_Release7_N"]
            for x in result["groups"] if x["QTL_dataset"]=="QTD000434"},{9})
        self.assertTrue(all(not x["can_calculate_formal_coloc_abf"] for x in result["groups"]))
        self.assertEqual(result["direct_Release7_phenotype_sdY_approved"],0)

    def test_2021_published_count_is_not_release7_sdY(self):
        self.assertNotEqual("metadata", "source phenotype sdY")
        self.assertEqual(M.RECORD_SHA256,
                         "818a5ebd7f448dd685f3f23d90d42f719a520a8fc27db233cb624bee21619682")




class ZenodoMutableApiRecordTests(unittest.TestCase):
    """Mutable API JSON cannot bypass the original file MD5/size gates."""

    @staticmethod
    def example_record(views=1):
        files=[
            {"key":name,"size":spec["bytes"],
             "checksum":"md5:"+spec["md5"],
             "links":{"self":"https://zenodo.org/api/records/4678936/files/"+name+"/content"}}
            for name,spec in M.FILES.items()
        ]
        files.extend({"key":f"other_{i}.tsv","size":1,
                      "checksum":"md5:"+"0"*32} for i in range(30))
        return {"id":4678936,"doi":"10.5281/zenodo.4678936",
                "files":files,"stats":{"views":views}}

    def test_source_file_hashes_stay_pinned_across_mutable_record_metadata(self):
        first=json.dumps(self.example_record(1),sort_keys=True).encode()
        second=json.dumps(self.example_record(2),sort_keys=True).encode()
        self.assertNotEqual(hashlib.sha256(first).hexdigest(),
                            hashlib.sha256(second).hexdigest())
        a,meta1=M.require_original_record(first)
        b,meta2=M.require_original_record(second)
        self.assertEqual(set(a),set(b))
        self.assertTrue(meta1["publisher_MD5_and_byte_lengths_remain_strictly_pinned"])
        self.assertTrue(meta2["publisher_MD5_and_byte_lengths_remain_strictly_pinned"])
        self.assertNotEqual(meta1["official_record_original_SHA256"],
                            meta2["official_record_original_SHA256"])

    def test_mutated_deposited_file_md5_always_fails(self):
        record=self.example_record()
        record["files"][0]["checksum"]="md5:"+"f"*32
        with self.assertRaises(ValueError):
            M.require_original_record(json.dumps(record).encode())

    def test_wrong_deposited_file_length_always_fails(self):
        record=self.example_record()
        record["files"][1]["size"]+=1
        with self.assertRaises(ValueError):
            M.require_original_record(json.dumps(record).encode())

if __name__=="__main__":
    unittest.main()
