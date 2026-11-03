import gzip
import hashlib
import importlib.util
import io
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/"scripts/is/audit_gigastroke_ais_original_ftp.py"
spec=importlib.util.spec_from_file_location("giga_ftp",S)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

FILE="GCST90104545.tsv.gz"
YAML=FILE+"-meta.yaml"
HEADER="chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\tp_value\teffect_allele_frequency"
ORIG=(HEADER+"\n1\t100\tG\tA\t0.2\t0.03\t0.001\t0.25\n1\t101\tC\tT\t-0.4\t0.05\t1e-5\t0.35\n").encode()

class Response:
    def __init__(self,url,data,head=False,redirect=False):
        self.url="https://evil.example/data" if redirect else url
        self.status=206 if "?" in url else 200
        self.headers={"Content-Length":str(len(data))}
        self.stream=io.BytesIO(b"" if head else data)
    def read(self,n=-1):
        return self.stream.read(n)
    def __enter__(self):return self
    def __exit__(self,*exc):self.close()
    def close(self):self.stream.close()

def fake_provider(original=ORIG, yaml_overrides="", header_only=False,
                  missing_yaml=False, invalid_host=False, extra_yaml=False):
    data=gzip.compress(original,mtime=0)
    md5=hashlib.md5(data).hexdigest()
    yaml=f"gwas_id: GCST90104545\ndata_file_name: {FILE}\ndata_file_md5sum: {md5}\nfile_type: GWAS-SSF v1.0\ncase_count: 19032\ncontrol_count: 237242\n"+yaml_overrides
    names=[FILE]
    if not missing_yaml:names.append(YAML)
    if extra_yaml:names.append("GCST90104545.other-meta.yaml")
    listing="<html>"+"".join(f'<a href="{name}">{name}</a>' for name in names)+"</html>"
    def opener(req,timeout=25):
        url=req.full_url
        if url==m.BASE:
            body=listing.encode()
        elif url.endswith(YAML):
            body=yaml.encode()
        elif url.endswith(FILE):
            body=data
        else:
            raise ValueError("unexpected mock URL")
        return Response(url,body,head=req.get_method()=="HEAD",redirect=invalid_host)
    return opener

class OriginalGIGASTROKEFTPTests(unittest.TestCase):
    def test_bounded_original_source_probe(self):
        out=m.probe(fake_provider())
        self.assertEqual(out["metadata_gate"],"METADATA_HEADER_SCREEN_PASS")
        self.assertEqual(out["yaml_selected_fields"]["case_count"],"19032")
        self.assertTrue(out["column_presence"]["standard_error"])
        self.assertFalse(out["original_file_byte_hash_verified"])
        self.assertEqual(out["scientific_status"],"HOLD")

    def test_full_hash_matches_official_yaml_but_not_scientific_claim(self):
        out=m.verify_full(m.probe(fake_provider()),fake_provider())
        self.assertEqual(out["full_file_gate"],"ORIGINAL_FILE_MD5_AND_GZIP_ROWS_PASS_REVIEW_REQUIRED")
        self.assertTrue(out["original_file_byte_hash_verified"])
        self.assertEqual(out["full_gwas_data_rows"],2)
        self.assertEqual(out["gwas_numeric_qc"]["finite_beta_se_p_positive_se_rows"],2)
        self.assertEqual(out["numeric_qc_status"],"BASIC_NUMERIC_FIELDS_SCREEN_PASS")
        self.assertFalse(out["original_gwas_N_SE_effect_allele_approved"])
        self.assertFalse(out["source_to_local_EAS5_canonical_allele_alignment_verified"])
        self.assertFalse(out["coloc_performed"])

    def test_negative_standard_error_is_not_numerically_approved(self):
        data=ORIG.replace(b"0.03",b"-0.03")
        opener=fake_provider(original=data)
        out=m.verify_full(m.probe(opener),opener)
        self.assertTrue(out["original_file_byte_hash_verified"])
        self.assertEqual(out["gwas_numeric_qc"]["finite_beta_se_p_positive_se_rows"],1)
        self.assertEqual(out["numeric_qc_status"],"HOLD_NUMERIC_OR_ALLELE_QC_REVIEW")
        self.assertFalse(out["original_gwas_N_SE_effect_allele_approved"])

    def test_malformed_eaf_is_not_approved(self):
        data=ORIG.replace(b"0.25",b"missing")
        opener=fake_provider(original=data)
        out=m.verify_full(m.probe(opener),opener)
        self.assertEqual(out["numeric_qc_status"],"HOLD_NUMERIC_OR_ALLELE_QC_REVIEW")
        self.assertFalse(out["original_gwas_N_SE_effect_allele_approved"])

    def test_missing_eaf_column_blocks_numeric_qc(self):
        data=ORIG.replace(b"\teffect_allele_frequency",b"\tmissing_eaf")
        opener=fake_provider(original=data)
        out=m.verify_full(m.probe(opener),opener)
        self.assertEqual(out["full_file_gate"],"HOLD_NUMERIC_QC_COLUMNS_MISSING")
        self.assertFalse(out["original_file_byte_hash_verified"])

    def test_gzip_row_scan_limit_is_fail_closed(self):
        opener=fake_provider()
        out=m.verify_full(m.probe(opener),opener,max_rows=1)
        self.assertEqual(out["full_file_gate"],"HOLD_EXCEEDS_ROW_SCAN_BOUND")
        self.assertFalse(out["original_file_byte_hash_verified"])

    def test_three_drive_reference_snps_independently_match(self):
        rows=(HEADER+"\n"
              "4\t81182554\tC\tT\t0.3\t0.05\t2.888e-8\t0.15\n"
              "10\t105447782\tA\tG\t0.2\t0.03\t1.332e-6\t0.25\n"
              "10\t105459834\tC\tT\t0.2\t0.03\t7.792e-7\t0.35\n").encode()
        provider=fake_provider(original=rows)
        out=m.verify_full(m.probe(provider),provider)
        self.assertEqual(out["original_to_drive_panel_gate"],
          "THREE_POSITION_ALLELE_AND_P_SCREEN_PASS_NOT_FULL_CANONICAL_VERIFICATION")
        self.assertTrue(all(x["unordered_alleles_and_p_match"] for x in
                            out["original_to_drive_three_anchor_panel"]))
        self.assertFalse(out["original_gwas_N_SE_effect_allele_approved"])
        self.assertFalse(out["source_to_local_EAS5_canonical_allele_alignment_verified"])

    def test_three_snp_panel_wrong_p_does_not_falsely_promote(self):
        rows=(HEADER+"\n"
              "4\t81182554\tC\tT\t0.3\t0.05\t2.888e-3\t0.15\n").encode()
        provider=fake_provider(original=rows)
        out=m.verify_full(m.probe(provider),provider)
        self.assertEqual(out["original_to_drive_panel_gate"],
          "HOLD_ORIGINAL_TO_DRIVE_3_POSITION_ANCHOR_MISMATCH")
        self.assertFalse(out["coloc_performed"])

    def test_32_drive_snp_manifest_is_frozen_and_valid(self):
        panel=m.load_extended_panel()
        self.assertEqual(len(panel),32)
        self.assertEqual(len(set(panel)),32)

    def test_32_snp_signed_beta_and_effect_allele_flip(self):
        panel=m.load_extended_panel()
        mock={}
        for key,row in panel.items():
            _,_,ref,alt=row["snp_GRCh37"].split(":")
            inverted=row["archive_allele_order_provenance"]=="OTHER_EFFECT_ORDER"
            mock[key]=[{
                "ea":ref if inverted else alt,
                "oa":alt if inverted else ref,
                "beta":-row["beta_ALT"] if inverted else row["beta_ALT"],
                "p":row["p"],
                "eaf":1-row["AF_ALT"] if inverted else row["AF_ALT"],
            }]
        result=m.compare_extended_panel(panel,mock)
        self.assertEqual(result["n_pass"],32)
        self.assertEqual(result["n_fail"],0)
        self.assertEqual(result["gate"],"32_SNPS_EXACT_REF_ALT_P_BETA_AF_SCREEN_PASS_NOT_FULL_ALIGNMENT")
        self.assertFalse(result["scientific_approved"])

    def test_32_snp_wrong_effect_direction_blocks_promotion(self):
        panel=m.load_extended_panel()
        row=next(iter(panel.values()))
        chrom,pos,ref,alt=row["snp_GRCh37"].split(":")
        data={(chrom,int(pos)):[{
            "ea":ref,"oa":alt,"beta":row["beta_ALT"],
            "p":row["p"],"eaf":row["AF_ALT"]
        }]}
        result=m.compare_extended_panel({(chrom,int(pos)):row},data)
        self.assertEqual(result["gate"],"HOLD_32_SNP_ORIENTATION_OR_SOURCE_MISMATCH")
        self.assertFalse(result["rows"][0]["checks"]["alt_oriented_beta_close"])

    def test_32_snp_duplicate_coordinate_blocks(self):
        panel=m.load_extended_panel()
        row=next(iter(panel.values()))
        chrom,pos,ref,alt=row["snp_GRCh37"].split(":")
        hit={"ea":alt,"oa":ref,"beta":row["beta_ALT"],
             "p":row["p"],"eaf":row["AF_ALT"]}
        result=m.compare_extended_panel({(chrom,int(pos)):row},
                      {(chrom,int(pos)):[hit,hit]})
        self.assertEqual(result["n_fail"],1)
        self.assertFalse(result["rows"][0]["checks"]["single_source_row"])

    def test_fullcis_6273_source_manifest_complete(self):
        sample=m.load_fullcis_source_panel()
        self.assertEqual(len(sample),6273)
        self.assertEqual(sum(not x["palindromic"] for x in sample.values()),5326)
        self.assertEqual(sum(x["palindromic"] for x in sample.values()),947)
        self.assertEqual(sum(key[0]=="4" for key in sample),3227)

    def test_fullcis_signed_alt_beta_and_af_both_allele_orders(self):
        sample=m.load_fullcis_source_panel()
        observed={}
        for key, ref in sample.items():
            alternate=ref["order"]=="OTHER_EFFECT_ORDER"
            observed[key]=[{
                "ea":ref["ref"] if alternate else ref["alt"],
                "oa":ref["alt"] if alternate else ref["ref"],
                "beta":-ref["beta_ALT"] if alternate else ref["beta_ALT"],
                "p":ref["p"],
                "eaf":1-ref["AF_ALT"] if alternate else ref["AF_ALT"],
            }]
        report=m.compare_fullcis_panel(sample,observed)
        self.assertEqual(report["screen_pass"],6273)
        self.assertEqual(report["screen_fail"],0)
        self.assertEqual(report["nonpalindromic_numerical_gate"],"PASS")
        self.assertEqual(report["palindromic_numeric_screen"],"PASS_ORIENTATION_STILL_REVIEW")
        self.assertEqual(report["scientific_status"],"HOLD")
        self.assertFalse(report["palindromic_strand_orientation_independently_resolved"])

    def test_fullcis_wrong_beta_is_hold(self):
        panel=m.load_fullcis_source_panel()
        key=next(k for k,v in panel.items() if not v["palindromic"] and abs(v["beta_ALT"])>0.03)
        ref=panel[key]
        observed={key:[{"ea":ref["alt"],"oa":ref["ref"],
                        "beta":-ref["beta_ALT"],"p":ref["p"],"eaf":ref["AF_ALT"]}]}
        report=m.compare_fullcis_panel({key:ref},observed)
        self.assertEqual(report["screen_fail"],1)
        self.assertEqual(report["full_reprocessed_ebi_cis_gwas_lineage_gate"],
          "HOLD_SOURCE_TO_DRIVE_DISCORDANCE")
        self.assertGreater(report["by_stratum"]["failed_beta_ALT_agrees"],0)

    def test_fullcis_duplicate_source_snp_is_hold(self):
        panel=m.load_fullcis_source_panel()
        key=next(iter(panel))
        ref=panel[key]
        item={"ea":ref["alt"],"oa":ref["ref"],"beta":ref["beta_ALT"],
              "p":ref["p"],"eaf":ref["AF_ALT"]}
        report=m.compare_fullcis_panel({key:ref},{key:[item,item]})
        self.assertEqual(report["screen_fail"],1)
        self.assertEqual(report["by_stratum"]["failed_unique_at_source_position"],1)

    def test_expected_md5_wrong_fails_source_verification(self):
        opener=fake_provider(yaml_overrides="data_file_md5sum: "+"0"*32+"\n")
        out=m.verify_full(m.probe(opener),opener)
        self.assertEqual(out["full_file_gate"],"HOLD_YAML_MD5_MISMATCH")
        self.assertFalse(out["original_file_byte_hash_verified"])

    def test_missing_se_stays_hold(self):
        orig=ORIG.replace(b"\tstandard_error",b"\tsome_other_info")
        out=m.probe(fake_provider(original=orig))
        self.assertEqual(out["metadata_gate"],"HOLD_REQUIRED_OR_SE_COLUMN_NOT_FOUND")

    def test_missing_metadata_yaml_stays_hold(self):
        out=m.probe(fake_provider(missing_yaml=True))
        self.assertEqual(out["metadata_gate"],"HOLD_ZERO_OR_MULTIPLE_YAML_FILES")

    def test_multiple_yaml_stays_hold(self):
        out=m.probe(fake_provider(extra_yaml=True))
        self.assertEqual(out["metadata_gate"],"HOLD_ZERO_OR_MULTIPLE_YAML_FILES")

    def test_wrong_origin_redirect_fails_closed(self):
        out=m.probe(fake_provider(invalid_host=True))
        self.assertEqual(out["metadata_gate"],"HOLD_NETWORK_OR_METADATA_PROBE")
        self.assertEqual(out["failure_type"],"ValueError")

    def test_reject_path_traversal(self):
        with self.assertRaises(ValueError):
            m.safe_asset_url("GCST90104545../secret.tsv.gz")

    def test_no_full_download_when_no_se(self):
        orig=ORIG.replace(b"\tstandard_error",b"\tx_other")
        opener=fake_provider(original=orig)
        out=m.verify_full(m.probe(opener),opener)
        self.assertEqual(out["full_file_gate"],"HOLD_HEADER_GATE_NOT_PASSED")

if __name__=="__main__":unittest.main()
