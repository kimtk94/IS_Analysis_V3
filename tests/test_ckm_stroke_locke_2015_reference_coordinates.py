"""Independent Locke 2015 EUR LDscore coordinate + 44-sentinel reference fixtures."""
from __future__ import annotations
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_locke_2015_reference_coordinates import (
    CANDIDATE_HEADER,read_candidate,parse_ldscore_member,
    scan_reference, audit,
)
from scripts.ckm_stroke.audit_eur_ldscore_zenodo import FILE_URL


def tab_header():
    return b"CHR SNP BP L2\n"


def mini_archive():
    out=io.BytesIO()
    panel={}
    with tarfile.open(fileobj=out,mode="w:gz") as tf:
        for ch in range(1,23):
            rows=[tab_header()]
            for i in range(1000):
                rows.append(f"{ch} other{ch}_{i} {5000+i} 1.0\n".encode("ascii"))
            for i in range(2):
                snp=f"rsTest{ch}_{i}"
                bp=20000+i
                rows.append(f"{ch} {snp} {bp} 1.0\n".encode("ascii"))
                panel[snp]=("A","C")
            compressed=gzip.compress(b"".join(rows),mtime=0)
            name=f"eur_w_ld_chr/{ch}.l2.ldscore.gz"
            meta=tarfile.TarInfo(name);meta.size=len(compressed)
            tf.addfile(meta,io.BytesIO(compressed))
            m_name=f"eur_w_ld_chr/{ch}.l2.M_5_50"
            scalar=b"1234\n"
            meta=tarfile.TarInfo(m_name);meta.size=len(scalar)
            tf.addfile(meta,io.BytesIO(scalar))
    return out.getvalue(),panel


class Response(io.BytesIO):
    def __init__(self,blob):
        super().__init__(blob)
        self.status=200
        self.headers={"Content-Length":str(len(blob))}
    def geturl(self):
        return FILE_URL


class IndependentCoordinatesTests(unittest.TestCase):
    def test_candidate_fractional_N_read_without_2018(self):
        with TemporaryDirectory() as td:
            f=Path(td)/"test.candidate.tsv.gz"; prov=Path(td)/"provenance.json"
            text=CANDIDATE_HEADER+"\n"+"".join(
                f"rsLocal{i}\t{100000+i/10}\t1.234\tA\tC\n" for i in range(120))
            raw=text.encode("ascii")
            f.write_bytes(gzip.compress(raw))
            prov.write_text(json.dumps({
                "technical_status":"PASS_2015_ORIGINAL_INDEPENDENT_SNP_PANEL_ONLY",
                "2018_exposure_candidates_used":False,
                "ready_for_ldsc":False,
                "counts":{"selected_rows":120},
                "2015_candidate_plaintext_sha256":hashlib.sha256(raw).hexdigest()
            }))
            ids,meta=read_candidate(f,prov,pin=False)
            self.assertEqual(meta["count"],120)
            self.assertIn("rsLocal21",ids)
            self.assertEqual(ids["rsLocal21"],("A","C"))
            self.assertEqual(meta["plaintext_sha256"],hashlib.sha256(raw).hexdigest())

    def test_unapproved_2018_candidate_is_rejected(self):
        with TemporaryDirectory() as td:
            f=Path(td)/"fake.candidate.tsv.gz";prov=Path(td)/"source.json"
            f.write_bytes(gzip.compress((CANDIDATE_HEADER+"\n").encode()))
            prov.write_text(json.dumps({"technical_status":"PASS_2015_ORIGINAL_INDEPENDENT_SNP_PANEL_ONLY",
                                         "2018_exposure_candidates_used":True,
                                         "ready_for_ldsc":False,"counts":{"selected_rows":120}}))
            with self.assertRaisesRegex(ValueError,"not authenticated"):
                read_candidate(f,prov,pin=False)

    def test_ldscore_snp_and_chromosome_position(self):
        panel={"rs1":("A","C"),"rs2":("A","C")}
        by_chr={k:[] for k in range(1,23)}
        rows=tab_header()+b"1 rs1 123 1.5\n1 rs2 125 2.0\n"
        seen=set();coords=set();counts=Counter()
        parse_ldscore_member(io.BytesIO(gzip.compress(rows)),1,panel,seen,coords,by_chr,counts,allow_small=True)
        self.assertEqual(seen,{"rs1","rs2"})
        self.assertEqual(coords,{(1,123),(1,125)})
        self.assertEqual(len(by_chr[1]),2)

    def test_wrong_ldscore_chr_rejected(self):
        panel={"rs1":("A","C")}
        by_chr={k:[] for k in range(1,23)}
        rows=tab_header()+b"2 rs1 123 1.5\n"
        with self.assertRaisesRegex(ValueError,"chromosome"):
            parse_ldscore_member(io.BytesIO(gzip.compress(rows)),1,panel,set(),set(),by_chr,Counter(),allow_small=True)

    def test_duplicate_ldscore_rsid_reported(self):
        panel={"rs1":("A","C")}
        by_chr={k:[] for k in range(1,23)}
        rows=tab_header()+b"1 rs1 123 1.5\n1 rs1 123 1.7\n"
        c=Counter()
        parse_ldscore_member(io.BytesIO(gzip.compress(rows)),1,panel,set(),set(),by_chr,c,allow_small=True)
        self.assertEqual(c["candidate_duplicate_ldscore_rsid"],1)

    def test_all_22_canonical_coords_and_44_sentinels(self):
        raw,panel=mini_archive()
        md5=hashlib.md5(raw).hexdigest()
        def opener(req,timeout):
            return Response(raw)
        def ensembl_fake(sample):
            self.assertEqual(len(sample),44)
            self.assertEqual({x["chr"] for x in sample},set(range(1,23)))
            return ([{"chr":v["chr"],"rsid":v["rsid"],
                      "orientation_classification":"reference_base_in_forward_allele_pair"}
                     for v in sample],
                    {"reference_base_in_forward_allele_pair":44})
        module="scripts.ckm_stroke.audit_locke_2015_reference_coordinates"
        with patch(module+".SIZE",len(raw)),patch(module+".EXPECTED_MD5",md5):
            result=scan_reference(panel,opener=opener,sequence_probe=ensembl_fake)
        self.assertEqual(result["technical_status"],
                         "PASS_2015_INDEPENDENT_CANDIDATE_22_CHR_REFERENCE_AND_44_SENTINELS_ONLY")
        self.assertEqual(result["counts"]["candidate_unique_ref_coordinates"],44)
        self.assertEqual(result["counts"]["candidate_missing_from_eur_ldscore"],0)
        self.assertEqual(result["sentinel_count"],44)
        self.assertFalse(result["ready_for_ldsc"])

    def test_bad_candidate_path_fails_closed(self):
        with TemporaryDirectory() as td:
            report=audit(Path(td)/"notthere.candidate.tsv.gz",Path(td)/"missing.json")
            self.assertEqual(report["technical_status"],"HOLD_REFERENCE_OR_POSITION")
            self.assertFalse(report["ready_for_ldsc"])
            self.assertFalse(report["2018_BMI_input_used"])


if __name__=="__main__":
    unittest.main()
