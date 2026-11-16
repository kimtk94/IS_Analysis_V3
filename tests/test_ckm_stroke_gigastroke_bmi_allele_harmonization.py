"""Synthetic GIGASTROKE 11-field allele/MD5/regression fixtures. No external network."""
from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import (
    COLUMNS,orientation,aligned_effect,examine_outcome_rows,scan_full,audit
)
from scripts.ckm_stroke.verify_gigastroke_raw import source_urls

ACC="GCST90104542"
PHEN="LAA"
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"
URL=BASE+ACC+"/"+ACC+"_buildGRCh37.tsv.gz"


def panel():
    return {
        (1,100):("rs1","A","C"),
        (1,200):("rs2","A","C"),
        (1,300):("rs3","A","C"),
        (1,400):("rs4","A","C"),
        (1,500):("rs5","A","C"),
        (1,600):("rs6","A","C"),
        (1,700):("rs7","A","C"),
    }


def row(pos,ea="A",oa="C",beta="0.1",se="0.02"):
    return "\t".join(["1",str(pos),"0.2",beta,se,"0.05","1.1","1.02","1.2",ea,oa])


class Response(io.BytesIO):
    def __init__(self,data,url=URL):
        super().__init__(data)
        self.status=200
        self.headers={"Content-Length":str(len(data))}
        self._url=url

    def geturl(self):
        return self._url


def publisher_metadata(md5):
    return f"""gwas_id: {ACC}
genome_assembly: GRCh37
data_file_name: {ACC}_buildGRCh37.tsv.gz
data_file_md5sum: {md5}
trait_description:
  - Large artery stroke
sample_size: 1241207
"""


def fake_data():
    rows=[
        row(100,"A","C"),row(200,"C","A"),row(300,"T","G"),
        row(400,"G","T"),row(500,"A","G"),
        row(600,"A","C"),row(600,"C","A"),row(999,"G","A")
    ]
    raw=("\t".join(COLUMNS)+"\n"+"\n".join(rows)+"\n").encode("ascii")
    gz=gzip.compress(raw)
    md5=hashlib.md5(gz).hexdigest()
    return gz,md5


def source_record(md5):
    return {
        "phenotype":PHEN,"accession":ACC,"reported_trait":"Large artery stroke",
        "build":"GRCh37","ancestry":"EUR",
        "sample_total":"1241207",
        "published_raw_md5":md5,
        "metadata_url":URL+"-meta.yaml",
        "verification":"METADATA_VERIFIED_RAW_HASH_UNCHECKED",
    }


def opener_for(gz,md5):
    def op(request,timeout):
        url=request.full_url
        if url.endswith("-meta.yaml"):return Response(publisher_metadata(md5).encode(),url)
        if url.endswith("/md5sum.txt"):
            return Response((md5+"  "+ACC+"_buildGRCh37.tsv.gz\n").encode(),url)
        if url==URL:return Response(gz,url)
        raise ValueError("unexpected URL "+url)
    return op


class HarmonizationTests(unittest.TestCase):
    def test_orientation_all_four_directions(self):
        self.assertEqual(orientation("A","C","A","C"),"same")
        self.assertEqual(orientation("A","C","C","A"),"swapped")
        self.assertEqual(orientation("A","C","T","G"),"complement_same")
        self.assertEqual(orientation("A","C","G","T"),"complement_swapped")

    def test_oriented_beta_se_sign_transform(self):
        self.assertEqual(aligned_effect(0.4,0.1,"same"),(0.4,4.0))
        self.assertEqual(aligned_effect(0.4,0.1,"swapped"),(-0.4,-4.0))
        self.assertEqual(aligned_effect(-0.4,0.1,"complement_swapped"),(0.4,4.0))
        with self.assertRaises(ValueError):
            aligned_effect(0.1,0.1,"allele_incompatible")

    def test_orientation_mismatch_and_ambiguous(self):
        self.assertEqual(orientation("A","C","A","G"),"allele_incompatible")
        self.assertEqual(orientation("A","T","A","T"),"palindromic_ambiguous")
        self.assertEqual(orientation("A","C","N","C"),"non_ACGT_alleles")

    def test_count_orientation_missing_and_duplicate_coordinate(self):
        gz,_=fake_data()
        with gzip.GzipFile(fileobj=io.BytesIO(gz),mode="rb") as f:
            def lines():
                for x in f:yield x.decode()
            r=examine_outcome_rows(lines(),panel())
        c=r["counts"]
        self.assertEqual(c["original_outcome_rows"],8)
        self.assertEqual(c["unique_same"],1)
        self.assertEqual(c["unique_swapped"],1)
        self.assertEqual(c["unique_complement_same"],1)
        self.assertEqual(c["unique_complement_swapped"],1)
        self.assertEqual(c["unique_allele_incompatible"],1)
        self.assertEqual(c["outcome_multiple_rows_same_BMI_coordinate"],1)
        self.assertEqual(c["outcome_collision_coordinates_excluded"],1)
        self.assertEqual(c["candidate_coordinates_not_in_outcome"],1)
        self.assertEqual(c["unique_harmonizable_BMI_positions"],4)
        self.assertEqual(c["unique_orientation_sign_positive"],2)
        self.assertEqual(c["unique_orientation_sign_negative"],2)

    def test_real_gzip_publisher_md5_pass_but_h2_hold(self):
        gz,md5=fake_data()
        with patch("scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization.EXPECTED_GZIP_SIZES",{PHEN:len(gz)}):
            r=scan_full(source_record(md5),panel(),fetcher=opener_for(gz,md5))
        self.assertEqual(r["technical_status"],"PASS_RAW_MD5_AND_ALLELE_CLASSIFICATION_ONLY")
        self.assertEqual(r["original_compressed_bytes"],len(gz))
        self.assertTrue(r["raw_md5_verified"])
        self.assertEqual(r["counts"]["unique_harmonizable_BMI_positions"],4)
        self.assertFalse(r["ready_for_ldsc"])
        self.assertFalse(r["outcome_per_variant_N_available"])

    def test_raw_source_md5_mismatch_is_blocking(self):
        gz,md5=fake_data()
        published="0"*32
        with patch("scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization.EXPECTED_GZIP_SIZES",{PHEN:len(gz)}):
            with self.assertRaisesRegex(ValueError,"MD5"):
                scan_full(source_record(published),panel(),fetcher=opener_for(gz,published))

    def test_broken_outcome_header_blocks(self):
        gz=gzip.compress(b"SNP\tP\nrs1\t0.1\n")
        md5=hashlib.md5(gz).hexdigest()
        with patch("scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization.EXPECTED_GZIP_SIZES",{PHEN:len(gz)}):
            with self.assertRaisesRegex(ValueError,"header"):
                scan_full(source_record(md5),panel(),fetcher=opener_for(gz,md5))

    def test_invalid_effect_numeric_excluded(self):
        raw="\t".join(COLUMNS)+"\n"+row(100,beta="NA")+"\n"+row(200,se="0")+"\n"
        r=examine_outcome_rows(iter(raw.splitlines(keepends=True)),panel())
        self.assertEqual(r["counts"]["unique_matched_numeric_invalid"],2)
        self.assertEqual(r["counts"]["unique_harmonizable_BMI_positions"],0)

    def test_source_error_does_not_fake_harmonization(self):
        from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import scan_full
        gz,md5=fake_data()
        with patch("scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization.EXPECTED_GZIP_SIZES",{PHEN:len(gz)+1}):
            with self.assertRaisesRegex(ValueError,"Content-Length"):
                scan_full(source_record(md5),panel(),fetcher=opener_for(gz,md5))


if __name__=="__main__":
    unittest.main()
