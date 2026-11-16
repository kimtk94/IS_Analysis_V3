#!/usr/bin/env python3
"""Original expression ZIP source tests, with synthetic private-data-free fixtures."""
import io
import sys
import tarfile
import unittest
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_EGEAD397_expression_source import (
    EXPECTED_BYTES,EXPECTED_MD5,filelist,parse_zip_index,parse_zip_members,
    tar_first_member,summarize_clinical,count_tar_sample_header_probe,
    confidential_matrix_clinical_overlap,
)


class ExpressionSourceTests(unittest.TestCase):
    def test_publisher_exact_source(self):
        source=(
            "#Archive/File\tName\tTime\tSize\tMD5\n"
            "Archive\tE-GEAD-397.processed.zip\t2021-05-13\t"
            f"{EXPECTED_BYTES}\t{EXPECTED_MD5}\n"
            "File\tclinical_diagnosis_age_sex_v2.txt\t2021-05-13\t8256\t"
            "0fd7179b660713a32a4f7657ee675162\n"
            "File\tcount.tar\t2021-05-13\t1310883840\t"
            "d8f1d3af683d4d362671be16ca3ceab9\n"
        ).encode()
        parsed=filelist(source)
        self.assertEqual(parsed["count.tar"]["source_length"],1310883840)

    def test_fake_source_hash_rejected(self):
        source=(
            "Archive\tE-GEAD-397.processed.zip\t2021\t"
            f"{EXPECTED_BYTES}\t{'0'*32}\n"
            "File\tclinical_diagnosis_age_sex_v2.txt\t2021\t8256\t"
            "0fd7179b660713a32a4f7657ee675162\n"
            "File\tcount.tar\t2021\t1310883840\t"
            "d8f1d3af683d4d362671be16ca3ceab9\n"
        ).encode()
        with self.assertRaisesRegex(ValueError,"drift"):
            filelist(source)

    def test_regular_original_zip_index(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w",compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr("clinical_diagnosis_age_sex_v2.txt","group\tcell\nHC\tCD16\n")
            z.writestr("count.tar","dummy-count-tar")
        raw=b.getvalue()
        f=parse_zip_index(raw,{"start":0,"end":len(raw)-1,"total":len(raw)})
        idx=raw[f["start"]:f["end"]+1]
        members=parse_zip_members(idx,f["entries"])
        self.assertEqual(len(members),2)
        self.assertEqual(members[0]["name"],"clinical_diagnosis_age_sex_v2.txt")
        self.assertEqual(members[1]["name"],"count.tar")

    def test_wrong_zip_signature_fails(self):
        with self.assertRaisesRegex(ValueError,"EOCD"):
            parse_zip_index(b"not a ZIP file"*80,{"start":0,"end":1119,"total":1120})

    def test_real_tar_member_name_no_content_export(self):
        tar=io.BytesIO()
        with tarfile.open(fileobj=tar,mode="w") as z:
            data=b"fake data"
            i=tarfile.TarInfo("counts/CD16p_Mono.txt")
            i.size=len(data)
            z.addfile(i,io.BytesIO(data))
        out=tar_first_member(tar.getvalue()[:512])
        self.assertEqual(out["first_member_name"],"counts/CD16p_Mono.txt")
        self.assertFalse(out["full_expression_count_matrix_rows_read"])

    def test_true_count_matrix_header_width_only_no_identifiers(self):
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode="w") as tar:
            data=b"GeneID\tDonor1\tDonor2\tDonor3\nENSG0001\t0\t9\t12\n"
            member=tarfile.TarInfo("CD16p_Mono_count.txt")
            member.size=len(data)
            tar.addfile(member,io.BytesIO(data))
        report=count_tar_sample_header_probe(stream.getvalue()[:8192])
        self.assertEqual(report["source_first_matrix_header_fields"],4)
        self.assertEqual(report["numeric_nonnegative_values_in_first_data_row_excluding_first_field"],3)
        self.assertTrue(report["count_matrix_column_count_not_equivalent_to_qtl_regression_N"])
        self.assertNotIn("Donor1",str(report))

    def test_anonymized_original_416_member_count_and_clinical_intersection(self):
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode="w") as tar:
            data=b"gene\tDonorA\tDonorB\tDonorC\nENSG0001\t0\t4\t7\n"
            m=tarfile.TarInfo("CD16p_Mono_count.txt")
            m.size=len(data)
            tar.addfile(m,io.BytesIO(data))
        clinical=b"id\tdisease\tPhase\tage\tsex\nDonorA\tHC\t1\t50\tF\nDonorB\tRA\t2\t40\tM\nDonorD\tHC\t1\t51\tF\n"
        r=confidential_matrix_clinical_overlap(stream.getvalue()[:8192],clinical)
        self.assertEqual(r["original_count_matrix_unique_sample_columns"],3)
        self.assertEqual(r["exact_same_public_name_match_columns_to_clinical"],2)
        self.assertEqual(r["unmatched_count_matrix_columns"],1)
        self.assertNotIn("DonorA",str(r))
        self.assertEqual(r["per_celltype_eQTL_regression_N"],"NOT_VERIFIED")

    def test_duplicate_matrix_source_names_must_fail(self):
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode="w") as tar:
            data=b"gene\tDonorA\tDonorA\nENSG0001\t0\t7\n"
            m=tarfile.TarInfo("CD16p_Mono_count.txt")
            m.size=len(data)
            tar.addfile(m,io.BytesIO(data))
        clinical=b"id\tdisease\nDonorA\tHC\n"
        with self.assertRaisesRegex(ValueError,"duplicated"):
            confidential_matrix_clinical_overlap(stream.getvalue()[:8192],clinical)

    def test_missing_original_count_matrix_data_line_fails(self):
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode="w") as tar:
            data=b"GeneID\tDonor1\tDonor2\n"
            member=tarfile.TarInfo("CD16p_Mono_count.txt")
            member.size=len(data)
            tar.addfile(member,io.BytesIO(data))
        with self.assertRaisesRegex(ValueError,"two complete count matrix lines"):
            count_tar_sample_header_probe(stream.getvalue()[:8192])

    def test_wrong_first_count_file_is_rejected(self):
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode="w") as tar:
            data=b"GeneID\tDonor1\tDonor2\nENSG\t1\t2\n"
            member=tarfile.TarInfo("Other_immune_count.txt")
            member.size=len(data)
            tar.addfile(member,io.BytesIO(data))
        with self.assertRaisesRegex(ValueError,"Unexpected first expression"):
            count_tar_sample_header_probe(stream.getvalue()[:8192])

    def test_clinical_summary_only_column_names(self):
        sample=("Donor ID\tDisease\tCell type\n"
                "anonymous1\tHC\tCD16p_Mono\n"
                "anonymous2\tRA\tCD16p_Mono\n").encode()
        out=summarize_clinical(sample)
        self.assertEqual(out["public_metadata_rows"],2)
        self.assertIn("Cell type",out["celltype_columns_detected"])
        self.assertNotIn("anonymous1",str(out))
        self.assertFalse(out["public_clinical_member_contains_direct_eQTL_per_celltype_N"])

if __name__=="__main__":
    unittest.main()
