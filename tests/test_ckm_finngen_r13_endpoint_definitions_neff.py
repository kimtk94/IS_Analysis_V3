"""Offline mutation and provenance tests for R13 phenotype vs GWAS availability."""
import hashlib
import io
import json
import unittest
from unittest.mock import patch
import zipfile
from scripts.ckm_stroke import audit_finngen_r13_endpoint_definitions_neff as r
from scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff import EXPECTED_SUBTYPES

def manifest():
    header="phenocode\tn_cases\tn_controls\n"
    src=["I9_STR\t36461\t446552\n",
         "I9_STR_EMBOLIC\t1150\t440000\n",
         "I9_STR_EMBOLIC_WITH_FLIMMER\t900\t440000\n"]
    src+=["DUMMY_TEST_%04d\t20\t1000\n"%i for i in range(2751)]
    return (header+"".join(src)).encode("utf-8")

def workbook(rows,header=("NAME","HD_ICD_10","INCLUDE","CONTROL_EXCLUDE")):
    def cell(i,j,value):
        col=chr(65+j)
        return '<c r="%s%d" t="inlineStr"><is><t>%s</t></is></c>'%(col,i,value)
    a=['<?xml version="1.0" encoding="UTF-8"?>',
       '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>']
    a.append('<row r="1">%s</row>'%"".join(cell(1,j,x) for j,x in enumerate(header)))
    for i,row in enumerate(rows,2):
        a.append('<row r="%d">%s</row>'%(i,"".join(cell(i,j,x) for j,x in enumerate(row))))
    a.append('</sheetData></worksheet>')
    buff=io.BytesIO()
    with zipfile.ZipFile(buff,"w",zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/worksheets/sheet1.xml","".join(a))
    return buff.getvalue()

def definition():
    return workbook([
        ("I9_STR_EXH","I63|I64","", "I9_STR_SAH"),
        ("I9_STR","I61|I63|I64","","I9_STR_EXH"),
        ("I9_STR_EMBOLIC","I63.1","","I9_CEREBVASC"),
        ("I9_STR_EMBOLIC_WITH_FLIMMER","","I9_STR_EMBOLIC",""),
        ("I9_ICH","I61","",""),
        ("I9_SAH","I60","",""),
    ])

def controls():
    return workbook([
        ("I9_STR_EXH","","","I9_ICH"),
        ("I9_STR","","","I9_STR_EXH"),
        ("I9_STR_EMBOLIC","","","I9_CEREBVASC"),
        ("I9_STR_EMBOLIC_WITH_FLIMMER","","","I9_AF"),
    ])

def cohort():
    items=[]
    for ph,ref in EXPECTED_SUBTYPES.items():
        n=ref["cases"];controls=1234808
        pooled=round(4*n*controls/(n+controls),3)
        known=round(.83*pooled,3)
        items.append({"phenotype":ph,"accession":ref["accession"],
                      "EUR_case_count_from_original_table1":n,
                      "EUR_control_count_from_original_table1":controls,
                      "illustrative_harmonic_neff":pooled,
                      "illustrative_cohortwise_neff":{
                          "cohortwise_sum_harmonic_neff_illustrative":known,
                          "cases_missing_cohort_control_count":0,
                          "not_per_variant_N_and_not_approved_for_LDSC":True}})
    return {"technical_status":"PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY",
            "publisher_source_sha256":r.COHORT_SOURCE_SHA,
            "publisher_original_table1_eur_cohort_count":44,
            "ready_for_ldsc":False,"LDSC_h2_rg_executed":False,"outcomes":items}

class TestOfficialEndpointSources(unittest.TestCase):
    def test_names_not_always_gwas_eligible(self):
        artifact=r.extract_endpoint_rows(definition(),"definitions")
        self.assertEqual(len(artifact["endpoint_definition_rows"]["I9_STR_EXH"]),1)
        self.assertIn("I61",artifact["endpoint_definition_rows"]["I9_STR"][0]["source_fields"]["HD_ICD_10"])

    def test_case_controls_not_equal_per_variant_neff(self):
        out=r.validated_cohort_sensitivity(cohort())
        self.assertEqual(len(out["scenarios"]),4)
        self.assertFalse(any(x["primary_LDSC_allowed"] for x in out["scenarios"]))
        self.assertTrue(out["cohort_and_pooled_scenarios_are_not_true_per_variant_Neff_bounds"])

    def test_cohort_sha_tamper_fails_closed(self):
        c=cohort();c["publisher_source_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"publisher-supplement cohort evidence mismatch"):
            r.validated_cohort_sensitivity(c)

    def test_invalid_cases_never_approved(self):
        c=cohort();c["outcomes"][0]["EUR_case_count_from_original_table1"]+=1
        with self.assertRaisesRegex(ValueError,"case/control"):
            r.validated_cohort_sensitivity(c)

    def test_manifests_source_sha_required(self):
        data=manifest()
        with self.assertRaisesRegex(ValueError,"SHA256"):
            r.manifest_review(data)

    def test_r13_absent_ischemic_file_remains_unavailable_even_with_definition(self):
        m=manifest()
        with patch.object(r,"R13_MANIFEST_SHA",hashlib.sha256(m).hexdigest()):
            out=r.review(m,definition(),controls(),cohort())
        self.assertEqual(out["technical_status"],
            "PASS_FG_R13_PUBLISHED_CASE_CONTROL_DEFINITIONS_AND_GWAS_AVAILABILITY_AUDIT")
        self.assertFalse(out["manifest"]["I9_STR_EXH_in_R13_GWAS_manifest"])
        self.assertTrue(out["I9_STR_EXH_GWAS_data_unavailable_if_not_listed_in_R13_manifest"])
        self.assertEqual(out["endpoints"][0]["case_definition_rows_found"],1)
        self.assertFalse(out["endpoints"][0]["present_in_public_R13_GWAS_manifest"])
        self.assertEqual(len(out["GIGASTROKE_meta_N_eff_input_sensitivity"]["scenarios"]),4)
        self.assertFalse(out["ready_for_ldsc"])
        self.assertFalse(out["LDSC_h2_intercept_rg_MR_coloc_performed"])

    def test_invalid_oo_xml_duplicate_zip_member_rejected(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("xl/worksheets/sheet1.xml","<a/>")
            z.writestr("xl/worksheets/sheet1.xml","<a/>")
        with self.assertRaisesRegex(ValueError,"duplicated"):
            r.extract_endpoint_rows(b.getvalue(),"bad")

    def test_not_actual_xlsx_zip_rejected(self):
        with self.assertRaisesRegex(ValueError,"ZIP"):
            r.extract_endpoint_rows(b"NO XLSX FILE","bad")

    def test_manifest_wrong_row_count_rejected(self):
        x=b"phenocode\nI9_STR\n"
        with patch.object(r,"R13_MANIFEST_SHA",hashlib.sha256(x).hexdigest()):
            with self.assertRaisesRegex(ValueError,"2754"):
                r.manifest_review(x)

    def test_manifest_explicit_counts_are_not_snp_counts(self):
        m=manifest()
        with patch.object(r,"R13_MANIFEST_SHA",hashlib.sha256(m).hexdigest()):
            out=r.manifest_review(m)
        self.assertEqual(out["requested_endpoint_rows"]["I9_STR"]["n_cases"],"36461")
        self.assertIn("n_cases",out["manifest_N_case_column_candidates"])
        self.assertTrue(out["manifest_is_not_per_variant_Neff"])

if __name__=="__main__":
    unittest.main()
