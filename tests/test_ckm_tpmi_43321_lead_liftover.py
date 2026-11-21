"""Official assembly reverse-identity and reference-allele source unit tests."""
import copy
import unittest
from scripts.ckm_stroke import map_tpmi_43321_lead_ensembl as m

def source():
    return {
      "schema":"ckm-TPMI43321-pinned-source-genomewide-significant-lead-v1",
      "technical_status":"PASS_PINNED_1739048_SOURCE_ROWS_ONE_SIGNIFICANT_VARIANT",
      "pinned_source_sha256_verified":m.SHA,
      "all_original_rows_verified":1739048,
      "original_source_genome_build":"GRCh38",
      "source_genomewide_significant_count":1,
      "original_significant_associations":[{
          "chrom_GRCh38":"11","position_GRCh38":6700000,
          "REF_GRCh38":"A","ALT_GRCh38":"G","A1_effect_allele":"G",
          "reported_rsid":"rs123456","P_original":1e-8,
      }],
      "scientific_status":"HOLD","ready_for_ldsc":False
    }
def mapping(chrom1, pos1, asm1, chrom2, pos2, asm2,strand1=1,strand2=1):
    return {"mappings":[{
      "original":{"assembly":asm1,"seq_region_name":chrom1,"start":pos1,"end":pos1,"strand":strand1},
      "mapped":{"assembly":asm2,"seq_region_name":chrom2,"start":pos2,"end":pos2,"strand":strand2},
    }]}
def http_mock(url,host):
    if "/map/human/GRCh38" in url:
        return mapping("11",6700000,"GRCh38","11",6590000,"GRCh37")
    if "/map/human/GRCh37" in url:
        return mapping("11",6590000,"GRCh37","11",6700000,"GRCh38")
    if "/sequence/" in url:
        return {"seq":"A"}
    raise AssertionError("unexpected API URL"+url)

class MappingTests(unittest.TestCase):
    def test_forward_and_reverse_roundtrip_reference_and_a1(self):
        r=m.source_to_reference(source(),reader=http_mock)
        self.assertEqual(r["technical_status"],"PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE")
        self.assertEqual(r["mapped_GRCh37"]["pos_GRCh37"],6590000)
        self.assertEqual(r["mapped_GRCh37"]["A1_effect_allele_GRCh37"],"G")
        self.assertFalse(r["scientific_replication_GO"])
        self.assertFalse(r["case_control_Neff_available"])
    def test_source_GWAS_SHA_tamper_fails(self):
        x=source();x["pinned_source_sha256_verified"]="0"*64
        with self.assertRaisesRegex(ValueError,"unverified"):
            m.source_gate(x)
    def test_wrong_original_build_rejected(self):
        x=source();x["original_source_genome_build"]="GRCh37"
        with self.assertRaisesRegex(ValueError,"unverified"):
            m.source_gate(x)
    def test_zero_or_negative_mapped_position_rejected(self):
        for pos in (0,-1):
            with self.subTest(pos=pos):
                with self.assertRaisesRegex(ValueError,"precision"):
                    m.parse_map(mapping("11",6700000,"GRCh38","11",pos,"GRCh37"),
                                "11",6700000,"GRCh38","GRCh37")
    def test_missing_or_multi_Ensembl_map_refused(self):
        with self.assertRaisesRegex(ValueError,"no unique"):
            m.parse_map({"mappings":[]}, "11",6700000,"GRCh38","GRCh37")
        with self.assertRaisesRegex(ValueError,"no unique"):
            m.parse_map({"mappings":[{},{}]}, "11",6700000,"GRCh38","GRCh37")
    def test_mismatched_reference_allele_fails(self):
        hit=source()["original_significant_associations"][0]
        fwd=m.parse_map(mapping("11",6700000,"GRCh38","11",6590000,"GRCh37"),"11",6700000,"GRCh38","GRCh37")
        with self.assertRaisesRegex(ValueError,"mapped reference"):
            m.check_alleles(hit,fwd,{"seq":"A"},{"seq":"T"})
    def test_reverse_strand_complement_correct(self):
        hit=source()["original_significant_associations"][0]
        fwd=m.parse_map(mapping("11",6700000,"GRCh38","11",6590000,"GRCh37",strand2=-1),
                         "11",6700000,"GRCh38","GRCh37")
        result=m.check_alleles(hit,fwd,{"seq":"A"},{"seq":"T"})
        self.assertEqual(result["REF_GRCh37"],"T")
        self.assertEqual(result["A1_effect_allele_GRCh37"],"C")
    def test_reverse_map_wrong_coordinates_hard_fails(self):
        fwd=m.parse_map(mapping("11",6700000,"GRCh38","11",6590000,"GRCh37"),
                         "11",6700000,"GRCh38","GRCh37")
        wrong=m.parse_map(mapping("11",6590000,"GRCh37","11",6700001,"GRCh38"),
                           "11",6590000,"GRCh37","GRCh38")
        with self.assertRaisesRegex(ValueError,"roundtrip"):
            m.validate_roundtrip(fwd,wrong)
    def test_snv_only_refuses_unnormalized_indel(self):
        x=source()["original_significant_associations"][0]
        x["REF_GRCh38"]="AA"
        x["A1_effect_allele"]="AA"
        with self.assertRaisesRegex(ValueError,"SNV"):
            m.check_alleles(x,{"mapping_orientation":1,"target_chrom":"11","target_pos":12},
                            {"seq":"AA"},{"seq":"AA"})
    def test_offhost_ensembl_map_rejected_without_request(self):
        with self.assertRaisesRegex(ValueError,"host"):
            m.official_ensembl("https://evil.example/map/human/GRCh38/11:1..1:1/GRCh37",
                               "rest.ensembl.org")

if __name__=="__main__":
    unittest.main()
