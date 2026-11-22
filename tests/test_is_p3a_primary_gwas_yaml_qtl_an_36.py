"""Tests for exact original GWAS-SSF release metadata and QTL an/2 science gate."""
import importlib
import io
import unittest

M=importlib.import_module("scripts.is.audit_p3a_primary_gwas_yaml_qtl_an_36")


def row(gene="ENSG00000065613",pos=103000000,an="774",ac="154",maf="0.198966"):
    d={x:"." for x in M.QTL.OFFICIAL_COLUMN_ORDER}
    d.update({"molecular_trait_id":gene,"gene_id":gene,
      "chromosome":"10","position":str(pos),"ref":"A","alt":"C",
      "variant":f"chr10_{pos}_A_C",
      "beta":"0.12","se":"0.06","pvalue":"0.00001","maf":maf,
      "type":"SNP","an":an,"ac":ac,"ma_samples":"101",
      "rsid":"rs123"})
    return "\t".join(d[k] for k in M.QTL.OFFICIAL_COLUMN_ORDER)


class FakeTabix:
    contigs=("10",)
    def __init__(self,rows):self.rows=rows
    def __enter__(self):return self
    def __exit__(self,*a):return None
    def fetch(self,chrom,start,end):
        assert chrom=="10"
        assert start==M.CIS.EXPECTED_WINDOW["start"]-1
        assert end==M.CIS.EXPECTED_WINDOW["end"]
        return iter(self.rows)


class FakeHTTP:
    status=200
    def __init__(self,data):self.data=data
    def __enter__(self):return self
    def __exit__(self,*a):return None
    def read(self,n):return self.data[:n]


class MetadataTests(unittest.TestCase):
    def test_valid_autosomal_even_an_to_variant_N(self):
        group={"ENSG00000065613"}
        result=M.aggregate_an("QTD000131",group,[row()])
        x=result["ENSG00000065613"]
        self.assertEqual(x["count_valid_nonpal_original_QTL"],1)
        self.assertEqual(x["count_an_valid_even_autosomal"],1)
        self.assertEqual(x["min_variant_an_over_2"],387)
        self.assertEqual(x["max_variant_an_over_2"],387)
        self.assertEqual(x["count_an_above_dataset_N"],0)
        self.assertTrue(x["has_single_constant_n_across_SNPs"])
        self.assertFalse(x["QTL_model_regression_effective_N_independently_approved"])

    def test_variable_variant_an_values_stay_variable_not_imputed(self):
        records=[row(pos=103000000,an="774"),
                 row(pos=103000001,an="770",ac="144",maf="0.187013")]
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},records)
        x=r["ENSG00000065613"]
        self.assertEqual(x["n_distinct_variant_an_over_2"],2)
        self.assertEqual(x["min_variant_an_over_2"],385)
        self.assertEqual(x["max_variant_an_over_2"],387)
        self.assertFalse(x["has_single_constant_n_across_SNPs"])
        self.assertFalse(x["QTL_model_regression_effective_N_independently_approved"])

    def test_invalid_odd_an_never_used(self):
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},
                         [row(an="773")])
        x=r["ENSG00000065613"]
        self.assertEqual(x["count_an_missing_or_invalid"],1)
        self.assertEqual(x["count_an_valid_even_autosomal"],0)
        self.assertFalse(x["all_valid_SNP_have_an_over_2"])

    def test_missing_an_not_equated_to_catalogue_387(self):
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},[row(an=".")])
        x=r["ENSG00000065613"]
        self.assertIsNone(x["min_variant_an_over_2"])
        self.assertEqual(x["count_an_missing_or_invalid"],1)

    def test_over_study_dataset_n_preserved_for_review(self):
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},[row(an="900")])
        x=r["ENSG00000065613"]
        self.assertEqual(x["count_an_above_dataset_N"],1)
        self.assertEqual(x["max_variant_an_over_2"],450)

    def test_ac_over_an_is_invalid_not_frequency_evidence(self):
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},
                         [row(ac="900")])
        self.assertEqual(r["ENSG00000065613"]["count_ac_an_invalid"],1)

    def test_ma_samples_larger_than_N_reported(self):
        s=row().replace("\t101\t","\t800\t")
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},[s])
        self.assertEqual(r["ENSG00000065613"]["count_ma_samples_invalid"],1)

    def test_source_duplicate_rsid_record_no_doubled_N_observation(self):
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},[row(),row()])
        self.assertEqual(r["ENSG00000065613"]["count_valid_nonpal_original_QTL"],1)

    def test_record_other_gene_not_misattributed(self):
        r=M.aggregate_an("QTD000131",{"ENSG00000065613"},
                         [row("ENSG00000099999")])
        self.assertEqual(r["ENSG00000065613"]["count_valid_nonpal_original_QTL"],0)

    def test_invalid_original_19_column_source_fail_closed(self):
        with self.assertRaises(ValueError):
            M.aggregate_an("QTD000131",{"ENSG00000065613"},["chr10\t1"])

    def test_synthetic_primary_GWAS_YAML_identity_and_counts(self):
        from unittest.mock import patch
        import yaml
        original=M.PINNED.SOURCES["EAS"]
        sample=(
          "gwas_id: GCST90104545\n"
          "genome_assembly: GRCh37\n"
          "data_file_name: GCST90104545_buildGRCh37.tsv.gz\n"
          "data_file_md5sum: "+original["publisher_md5_pinned"]+"\n"
          "trait_description: Ischemic stroke\n"
          "sample_size: 256274\n"
          "case_control_study: true\n"
          "case_count: 19032\n"
          "control_count: 237242\n"
        )
        result=M.fetch_primary_gwas_metadata(
            "EAS",open_url=lambda request,timeout:FakeHTTP(sample.encode()))
        self.assertTrue(result["case_control_sum_reconciles_primary_file"])
        self.assertTrue(result["followup_publication_case_control_numbers_match"])
        self.assertEqual(result["declared_sample_counts"]["case_count"],19032)
        self.assertFalse(result["scientific_case_fraction_and_N_approved_for_coloc"])

    def test_wrong_publisher_original_md5_rejected(self):
        sample=("gwas_id: GCST90104545\n"
          "genome_assembly: GRCh37\n"
          "data_file_name: GCST90104545_buildGRCh37.tsv.gz\n"
          "data_file_md5sum: 00000000000000000000000000000000\n")
        with self.assertRaises(ValueError):
            M.fetch_primary_gwas_metadata(
                "EAS",open_url=lambda request,timeout:FakeHTTP(sample.encode()))

    def test_publisher_missing_case_counts_stays_scientifically_held(self):
        spec=M.PINNED.SOURCES["EAS"]
        sample=("gwas_id: GCST90104545\n"
          "genome_assembly: GRCh37\n"
          "data_file_name: GCST90104545_buildGRCh37.tsv.gz\n"
          "data_file_md5sum: "+spec["publisher_md5_pinned"]+"\n"
          "sample_size: 256274\n")
        r=M.fetch_primary_gwas_metadata(
           "EAS",open_url=lambda request,timeout:FakeHTTP(sample.encode()))
        self.assertEqual(r["source_status"],"HOLD_PUBLISHER_CASE_CONTROL_MISSING")
        self.assertFalse(r["primary_GWAS_N_metadata_technically_complete"])
        self.assertFalse(r["scientific_case_fraction_and_N_approved_for_coloc"])

    def test_metadata_counts_but_wrong_total_not_approved(self):
        spec=M.PINNED.SOURCES["EAS"]
        sample=("gwas_id: GCST90104545\n"
          "genome_assembly: GRCh37\n"
          "data_file_name: GCST90104545_buildGRCh37.tsv.gz\n"
          "data_file_md5sum: "+spec["publisher_md5_pinned"]+"\n"
          "sample_size: 260000\n"
          "case_count: 19032\n"
          "control_count: 237242\n")
        r=M.fetch_primary_gwas_metadata(
           "EAS",open_url=lambda request,timeout:FakeHTTP(sample.encode()))
        self.assertFalse(r["case_control_sum_reconciles_primary_file"])

    def test_qtl_and_gwas_n_not_interchangeable(self):
        self.assertEqual(M.EXPECTED_SAMPLES,
           {"QTD000131":387,"QTD000171":205,"QTD000434":560})
        self.assertEqual(M.EXPECTED_LITERATURE["EUR"]["case_count"],62100)


if __name__=="__main__":
    unittest.main()
