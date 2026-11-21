"""Network-free pinpoint original 433.21 significant variant provenance tests."""
import gzip
import hashlib
import io
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
from scripts.ckm_stroke import identify_tpmi_43321_lead_variant as m

HEADER="\t".join(m.EXPECTED_HEADER)
def line(n=1,p="0.5",chrom="1",ref="A",alt="G",a1="G",obs="204099",orv="1.02",se="0.05",z="0.396",err="."):
    # P=1e-8 is generally incompatible with OR 1.02 Z=.396; set significant values consistently.
    if p=="1e-8":
        orv="1.4";se="0.059013";z="5.700"
    return "\t".join([chrom,str(100000+n),f"{chrom}:{100000+n}:{ref}:{alt}",
                      ref,alt,"N",a1,"NA","0.3","ADD",obs,orv,se,z,p,err,"rs123456"])
def fake(n=5,lead=True):
    body=[line(i,p="1e-8" if lead and i==3 else "0.7") for i in range(1,n+1)]
    # P of non-significant lines need not match z in scan; significant parse checks only hit.
    return gzip.compress((HEADER+"\n"+"\n".join(body)+"\n").encode())
class Response:
    status=200
    def __init__(self,data,url=m.URL):
        self.data=io.BytesIO(data)
        self.headers={"Content-Length":str(len(data))}
        self.url=url
    def __enter__(self):return self
    def __exit__(self,*args):return None
    def geturl(self):return self.url
    def read(self,n=-1):return self.data.read(n)
class TestLead(unittest.TestCase):
    def test_extract_one_true_significant_SNP(self):
        b=fake()
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"source.gz";p.write_bytes(b)
            with patch.object(m,"EXPECTED_ROWS",5):
                r=m.scan_original(p)
        self.assertEqual(r["original_significant_associations"][0]["reported_rsid"],"rs123456")
        self.assertEqual(r["original_significant_associations"][0]["position_GRCh38"],100003)
        self.assertEqual(r["source_genomewide_significant_count"],1)
        self.assertFalse(r["GIGASTROKE_EUR_GRCh37_same_variant_harmonization_verified"])
        self.assertFalse(r["ready_for_ldsc"])
    def test_source_SHA_must_match_prior_original_run(self):
        b=fake()
        with patch.object(m,"EXPECTED_ROWS",5),patch.object(m,"EXPECTED_SHA256",hashlib.sha256(b).hexdigest()):
            r=m.download_once(opener=lambda req,timeout:Response(b),expected_bytes=len(b))
        self.assertEqual(r["source_download_count"],1)
        self.assertTrue(r["temporary_full_file_deleted"])
        self.assertFalse(r["entire_GWAS_summary_statistics_NOT_redistributed"] is False)
    def test_original_source_checksum_mutation_rejected(self):
        b=fake()
        with self.assertRaisesRegex(ValueError,"SHA256"):
            m.download_once(opener=lambda req,timeout:Response(b),expected_bytes=len(b))
    def test_significance_count_changed_fails_closed(self):
        b=fake(lead=False)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"src.gz";p.write_bytes(b)
            with patch.object(m,"EXPECTED_ROWS",5):
                with self.assertRaisesRegex(ValueError,"single genomewide"):
                    m.scan_original(p)
    def test_significant_non_source_test_allele_rejected(self):
        cols=line(3,p="1e-8",a1="C").split("\t")
        with self.assertRaisesRegex(ValueError,"genotype"):
            m.parse_significant(cols)
    def test_significant_identical_REF_ALT_refused(self):
        cols=line(3,p="1e-8",ref="A",alt="A",a1="A").split("\t")
        with self.assertRaisesRegex(ValueError,"genotype"):
            m.parse_significant(cols)
    def test_significant_non_ACGT_allele_refused(self):
        cols=line(3,p="1e-8",ref="N",alt="G",a1="G").split("\t")
        with self.assertRaisesRegex(ValueError,"genotype"):
            m.parse_significant(cols)
    def test_significant_missing_or_invalid_N_refused(self):
        cols=line(3,p="1e-8",obs="99999").split("\t")
        with self.assertRaisesRegex(ValueError,"N mismatch"):
            m.parse_significant(cols)
    def test_pinned_original_header_change_rejected(self):
        b=fake()
        src=gzip.decompress(b).decode().replace("LOG(OR)_SE","SE")
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"src.gz";p.write_bytes(gzip.compress(src.encode()))
            with self.assertRaisesRegex(ValueError,"header changed"):
                m.scan_original(p)
    def test_no_multi_hit_or_p_zero_can_be_promoted(self):
        b=fake(5)
        txt=gzip.decompress(b).decode()
        txt=txt.replace("\t0.7\t.\trs123456","\t1e-8\t.\trs123456",1)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"src.gz";p.write_bytes(gzip.compress(txt.encode()))
            with patch.object(m,"EXPECTED_ROWS",5):
                with self.assertRaises(ValueError):
                    m.scan_original(p)

if __name__=="__main__":
    unittest.main()
