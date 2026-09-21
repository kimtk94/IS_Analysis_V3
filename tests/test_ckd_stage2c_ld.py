#!/usr/bin/env python3
import importlib.util
from pathlib import Path

P = Path(__file__).resolve().parents[1] / "scripts" / "prepare_ckd_stage2c_ld.py"
spec = importlib.util.spec_from_file_location("m", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_choose_region_block():
    blocks={"1":[(100,200),(201,400)]}
    lo,hi,src,rawlo,rawhi,trunc=m.choose_region("1",150,[80,90,100,150,200,250],blocks,500)
    assert (lo,hi)==(100,200)
    assert src=="Berisa-Pickrell_EUR_hg19"
    assert trunc==0

def test_choose_region_truncated():
    blocks={"1":[(50,300)]}
    lo,hi,src,rawlo,rawhi,trunc=m.choose_region("1",150,[100,150,200],blocks,500)
    assert (lo,hi)==(100,200)
    assert trunc==1

def test_match_reference_orientation():
    summary=[
      {"pos37":"100","allele0_pqtl":"G","allele1_pqtl":"A","snp":"rs1"},
      {"pos37":"200","allele0_pqtl":"T","allele1_pqtl":"C","snp":"rs2"},
    ]
    pvar=[
      {"pos":100,"id":"1:100:G:A","ref":"G","alt":"A","chrom":"1"},
      {"pos":200,"id":"1:200:C:T","ref":"C","alt":"T","chrom":"1"},
    ]
    got,missing,amb=m.match_reference(summary,pvar)
    assert missing==0 and amb==0 and len(got)==2
    # rs1 effect allele=A=ALT -> inverse of REF-based LD; rs2 effect=C=REF.
    assert got[0]["_ld_sign"]==-1
    assert got[1]["_ld_sign"]==1
