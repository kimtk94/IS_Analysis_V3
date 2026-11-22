"""Synthetic fail-closed gates for UCSC official liftover/FASTA anchors."""
import copy
import hashlib
import importlib
import json
import tempfile
import unittest
from pathlib import Path

M=importlib.import_module("scripts.is.audit_p3a_chr10_ucsc_liftover_fasta")
MANIFEST=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")


class DNA:
    def __init__(self,shift=0,flip_lead=False):
        self.shift=shift
        self.flip_lead=flip_lead
    def __len__(self):return 150_000_000
    def __getitem__(self,key):
        if key==M.REGION["historical_lead_GRCh37_1based"]-1+self.shift:
            return ord("C" if self.flip_lead else "G")
        return ord("A")


class FakeLift:
    def __init__(self,shift=1000,duplicate=False,strand="+",bad=False):
        self.shift=shift
        self.duplicate=duplicate
        self.strand=strand
        self.bad=bad
    def convert_coordinate(self,chrom,pos):
        out=[("chr10",pos+self.shift,self.strand,1200)]
        if self.duplicate and pos==M.REGION["historical_lead_GRCh37_1based"]-1:
            out.append(("chr10",pos+self.shift+1,self.strand,4))
        if self.bad and pos==M.ANCHORS["cis_region_midpoint"]-1:
            return []
        return out


def manifest():
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


class UCSCTest(unittest.TestCase):
    def test_six_anchor_ref_pass_still_zero_coloc(self):
        a=M.anchors_audit(FakeLift(),DNA(),DNA(1000))
        self.assertEqual(len(a["mapped_anchors"]),6)
        self.assertTrue(a["all_unique_hg38_chr10_plus_monotonic_REF_concordant"])
        self.assertEqual(a["mapped_anchors"]["palindromic_lead"]["hg19_REF"],"G")
        self.assertTrue(a["mapped_anchors"]["palindromic_lead"]["archived_allele_pair_ref_member"])
        r=M.audit_results(manifest(),{"chain":{"verified":True}},{"GRCh37_chr10":"synthetic"},a)
        self.assertEqual(r["n_real_group_keys"],36)
        self.assertEqual(r["n_distinct_genes"],20)
        self.assertEqual(r["formal_coloc_abf_eligible"],0)
        self.assertEqual(r["formal_coloc_susie_executed"],0)
        self.assertEqual(len(r["groups"]),36)
        self.assertTrue(all(x["coloc_abf_eligible"] is False for x in r["groups"]))

    def test_six_anchor_map_is_not_full_cis_proof(self):
        x=M.anchors_audit(FakeLift(),DNA(),DNA(1000))
        self.assertFalse(x["full_cis_region_equivalence_verified"])
        self.assertFalse(x["all_variant_liftover_verified"])
        self.assertFalse(x["candidate_query_region_GRCh38_1based"][
            "is_lower_bound_of_full_cis_with_guaranteed_coverage"])

    def test_duplicated_lead_position_denied(self):
        x=M.anchors_audit(FakeLift(duplicate=True),DNA(),DNA(1000))
        self.assertEqual(x["anchor_status"],"HOLD_COORDINATE_OR_REFERENCE_DISCORDANCE")
        self.assertIsNone(x["candidate_query_region_GRCh38_1based"])

    def test_unmapped_middle_denied(self):
        x=M.anchors_audit(FakeLift(bad=True),DNA(),DNA(1000))
        self.assertFalse(x["all_unique_hg38_chr10_plus_monotonic_REF_concordant"])

    def test_lead_REF_discordance_denied(self):
        x=M.anchors_audit(FakeLift(),DNA(),DNA(1000,flip_lead=True))
        self.assertEqual(x["mapped_anchors"]["palindromic_lead"]["status"],
                         "REFERENCE_MISMATCH")
        self.assertFalse(x["all_unique_hg38_chr10_plus_monotonic_REF_concordant"])

    def test_wrong_chr_mapping_denied(self):
        class Wrong(FakeLift):
            def convert_coordinate(self,chrom,pos):
                return [("chr7",pos+1000,"+",1200)]
        x=M.anchors_audit(Wrong(),DNA(),DNA(1000))
        self.assertFalse(x["all_unique_hg38_chr10_plus_monotonic_REF_concordant"])

    def test_reverse_strand_requires_adjudication(self):
        x=M.anchors_audit(FakeLift(strand="-"),DNA(),DNA(1000))
        self.assertFalse(x["all_unique_hg38_chr10_plus_monotonic_REF_concordant"])

    def test_nonpinned_download_rejected(self):
        class Stream:
            status=200
            def __init__(self):self.remain=140000
            def __enter__(self):return self
            def __exit__(self,*args):return None
            def read(self,n):
                if not self.remain:return b""
                take=min(self.remain,n)
                self.remain-=take
                return b"x"*take
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                M.fetch_checked("chain",Path(d)/"x.gz",
                                lambda req,timeout:Stream())

    def test_unexpected_source_group_flag_denied(self):
        x=manifest()
        x["groups"][0]["approved_for_new_coloc"]=True
        with self.assertRaises(ValueError):
            M.audit_results(x,{}, {},{"all_unique_hg38_chr10_plus_monotonic_REF_concordant":False})

    def test_extra_QTL_group_denied(self):
        x=manifest()
        x["groups"][0]["QTL_dataset"]="QTD000434"
        with self.assertRaises(ValueError):
            M.audit_results(x,{}, {},{"all_unique_hg38_chr10_plus_monotonic_REF_concordant":False})

    def test_coordinate_conversion_zero_based_explicit(self):
        a=M.single_coordinate(FakeLift(shift=17),"chr10",105596641)
        self.assertEqual(a["source_position_1based"],105596641)
        self.assertEqual(a["target_position_1based"],105596658)

    def test_official_known_digests_are_full_md5(self):
        for x in M.FILES.values():
            self.assertEqual(len(x["md5"]),32)
            self.assertGreater(x["limit"],100_000)

if __name__=="__main__":
    unittest.main()
