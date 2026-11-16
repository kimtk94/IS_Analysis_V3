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
    tar_first_member,summarize_clinical,
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
