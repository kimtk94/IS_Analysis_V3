#!/usr/bin/env python3
import importlib.util
import unittest
from pathlib import Path

P = Path(__file__).resolve().parents[1] / "scripts" / "prepare_ckd_stage2c_ld.py"
spec = importlib.util.spec_from_file_location("m", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class Stage2CLDTests(unittest.TestCase):
    def test_choose_region_block(self):
        blocks={"1":[(100,200),(201,400)]}
        lo,hi,src,rawlo,rawhi,trunc=m.choose_region("1",150,[80,90,100,150,200,250],blocks,500)
        self.assertEqual((lo,hi),(100,200))
        self.assertEqual(src,"Berisa-Pickrell_EUR_hg19")
        self.assertEqual(trunc,0)

    def test_choose_region_truncated(self):
        blocks={"1":[(50,300)]}
        lo,hi,src,rawlo,rawhi,trunc=m.choose_region("1",150,[100,150,200],blocks,500)
        self.assertEqual((lo,hi),(100,200))
        self.assertEqual(trunc,1)

    def test_write_sample_list_uses_real_newlines(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "samples.txt"
            m.write_sample_list(p, ["HG00096", "NA12878"])
            self.assertEqual(p.read_bytes(), b"HG00096\nNA12878\n")
            self.assertEqual(p.read_text().splitlines(), ["HG00096", "NA12878"])

    def test_phase3_vcf_url(self):
        self.assertEqual(
            m.phase3_vcf_url("6"),
            "https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        )

    def test_parse_afreq(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.afreq"
            p.write_text("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\tv1\tG\tA\t0.8\t1000\n1\tv2\tC\tT\t0.2\t1000\n")
            got = m.parse_afreq(p)
            self.assertAlmostEqual(got["v1"]["alt_freq"], 0.8)
            self.assertEqual(got["v1"]["ref"], "G")
            self.assertEqual(got["v1"]["alt"], "A")

    def test_match_reference_orientation(self):
        summary=[
          {"pos37":"100","allele0_pqtl":"G","allele1_pqtl":"A","snp":"rs1"},
          {"pos37":"200","allele0_pqtl":"T","allele1_pqtl":"C","snp":"rs2"},
        ]
        pvar=[
          {"pos":100,"id":"1:100:G:A","ref":"G","alt":"A","chrom":"1"},
          {"pos":200,"id":"1:200:C:T","ref":"C","alt":"T","chrom":"1"},
        ]
        got,missing,amb=m.match_reference(summary,pvar)
        self.assertEqual((missing,amb,len(got)),(0,0,2))
        self.assertEqual(got[0]["_effect_vs_ref_sign"],-1)
        self.assertEqual(got[1]["_effect_vs_ref_sign"],1)


if __name__ == "__main__":
    unittest.main()
