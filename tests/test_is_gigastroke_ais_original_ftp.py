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
HEADER="chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\tp_value"
ORIG=(HEADER+"\n1\t100\tG\tA\t0.2\t0.03\t0.001\n1\t101\tC\tT\t-0.4\t0.05\t1e-5\n").encode()

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
    data=gzip.compress(original)
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
        self.assertFalse(out["original_gwas_N_SE_effect_allele_approved"])
        self.assertFalse(out["source_to_local_EAS5_canonical_allele_alignment_verified"])
        self.assertFalse(out["coloc_performed"])

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
