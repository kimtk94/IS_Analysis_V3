"""Synthetic original-paper table/manifest contracts. No network or GWAS data."""
from __future__ import annotations
import hashlib
import io
from unittest.mock import patch
from zipfile import ZipFile
import unittest

from scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff import (
    source_evidence, workbook_cells, EXPECTED_SUBTYPES,
)


def fixture():
    t1={i:{"B":f"EUR_COHORT_{i}"} for i in range(5,49)}
    t1[3]={"J":"Incident & prevalent combined"}
    t1[4]={"J":"AS","K":"AIS","L":"CES","M":"SVS","N":"LAS","O":"ICH","P":"Controls"}
    t1[15]={"B":"UKBiobank","K":"2211","P":"369419",
            "L":"NA","M":"NA","N":"NA"}
    t1[5].update({"K":"59889","L":"10804","M":"6811","N":"6399","P":"865389"})
    t1[49]={"A":"TOTAL EUR","J":"73652","K":"62100","L":"10804",
            "M":"6811","N":"6399","P":"1234808"}
    t3={
        70:{"C":"AIS"},
        71:{"A":"EUR","C":"1.1","D":"1.08","E":"1.06","F":"1.09"},
        79:{"C":"AIS"},
        80:{"A":"EUR","C":"4.9E-3 (5E-4)","D":"3.3E-3 (1.1E-3)",
            "E":"7.4E-3 (1.2E-3)","F":"2E-3 (1.2E-3)"}
    }
    return t1,t3


def manifest():
    text="phenotype\taccession\tbuild\tsample_total\n"
    for phenotype,ref in EXPECTED_SUBTYPES.items():
        text+=f'{phenotype}\t{ref["accession"]}\tGRCh37\t{ref["total"]}\n'
    return text


class CohortAuditTests(unittest.TestCase):
    def mocked(self,t1=None,t3=None,manifest_text=None):
        a,b=fixture()
        with patch("scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff.workbook_cells",
                   return_value=(t1 or a,t3 or b)):
            return source_evidence(b"",manifest_text or manifest())

    def test_pass_primary_cohorts_and_published_ldsc_benchmarks(self):
        r=self.mocked()
        self.assertEqual(r["technical_status"],"PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY")
        self.assertEqual(r["publisher_original_table1_eur_cohort_count"],44)
        self.assertEqual(r["uk_biobank"]["AIS_case_count"],2211)
        self.assertEqual(r["uk_biobank"]["UKB_reported_control_count"],369419)
        self.assertFalse(r["uk_biobank"]["overlap_in_BMI_actual_count_known"])
        self.assertEqual(r["outcomes"][1]["accession"],"GCST90104542")
        self.assertEqual(r["outcomes"][1]["EUR_case_count_from_original_table1"],6399)
        self.assertEqual(r["outcomes"][2]["accession"],"GCST90104541")
        self.assertEqual(r["outcomes"][2]["original_CovLDSC_intercept_from_table3B"],1.06)
        self.assertEqual(r["outcomes"][1]["original_CovLDSC_observed_h2_SE"],0.0012)
        self.assertTrue(all(x["illustrative_neff_not_per_snp_or_approved"] for x in r["outcomes"]))
        self.assertFalse(r["ready_for_ldsc"])
        self.assertEqual(r["outcomes"][0]["illustrative_cohortwise_neff"]["reported_case_total_from_rows"],62100)
        self.assertTrue(r["outcomes"][1]["illustrative_cohortwise_neff"]["not_per_variant_N_and_not_approved_for_LDSC"])

    def test_manifest_swapped_LAA_CE_accessions_fails(self):
        m=manifest().replace("GCST90104542","WRONG_TMP").replace(
            "GCST90104541","GCST90104542").replace("WRONG_TMP","GCST90104541")
        with self.assertRaisesRegex(ValueError,"accession"):
            self.mocked(manifest_text=m)

    def test_bad_manifest_total_fails(self):
        m=manifest().replace("1241207","1241208")
        with self.assertRaisesRegex(ValueError,"GWAS YAML"):
            self.mocked(manifest_text=m)

    def test_ukb_cases_changed_fails(self):
        a,b=fixture()
        a[15]["K"]="2212"
        with self.assertRaisesRegex(ValueError,"UK Biobank"):
            self.mocked(t1=a,t3=b)

    def test_missing_original_EUR_controls_fails(self):
        a,b=fixture()
        a[49]["P"]="1234807"
        with self.assertRaisesRegex(ValueError,"control total"):
            self.mocked(t1=a,t3=b)

    def test_wrong_original_ldsc_table_schema_fails(self):
        a,b=fixture()
        b[70]["C"]="CE"
        with self.assertRaisesRegex(ValueError,"LDSC EUR schema"):
            self.mocked(t1=a,t3=b)

    def test_full_publisher_sha_pin_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"pin mismatch"):
            workbook_cells(b"unverified xlsx")

    def test_minimal_xlsx_parser_inlines_and_sheet_identification(self):
        ns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        wb=f'<workbook xmlns="{ns}"><sheets>'+''.join(
            f'<sheet name="{name}" sheetId="{i+1}"/>'
            for i,name in enumerate(["Table.of.content","Suppl.Table1","Suppl.Table2","Suppl.Table3"])
        )+'</sheets></workbook>'
        shared=f'<sst xmlns="{ns}"></sst>'
        s2=f'<worksheet xmlns="{ns}"><sheetData><row r="15"><c r="B15" t="inlineStr"><is><t>UKBiobank</t></is></c><c r="K15"><v>2211</v></c></row></sheetData></worksheet>'
        s4=f'<worksheet xmlns="{ns}"><sheetData><row r="71"><c r="C71"><v>1.10</v></c></row></sheetData></worksheet>'
        buf=io.BytesIO()
        with ZipFile(buf,"w") as z:
            for name,data in [
                ("xl/workbook.xml",wb),("xl/sharedStrings.xml",shared),
                ("xl/worksheets/sheet2.xml",s2),
                ("xl/worksheets/sheet4.xml",s4),
            ]:z.writestr(name,data.encode())
        payload=buf.getvalue()
        with patch("scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff.EXPECTED_SIZE",len(payload)),patch(
            "scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff.EXPECTED_SHA256",
            hashlib.sha256(payload).hexdigest()
        ):
            table1,table3=workbook_cells(payload)
        self.assertEqual(table1[15]["B"],"UKBiobank")
        self.assertEqual(table1[15]["K"],"2211")
        self.assertEqual(table3[71]["C"],"1.10")


if __name__=="__main__":
    unittest.main()
