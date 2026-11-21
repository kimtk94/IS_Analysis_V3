"""Fail-closed SAIGE allele orientation/QC and streaming-join fixtures."""
import sqlite3
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from scripts.ckm_stroke.tpmi_saige_eas_ais_full_allele_overlap import (
 SAIGE_HEADER, saige_candidate,make_sqlite,compare,
)
def row():
 out=["1","10353","1:10353:A:G","A","G","100",
      "0.40","0.98","0.15","0.02","1","2","0.001",
      "0.001","true","0.39","0.41","9249","304660","PASS"]
 assert len(out)==len(SAIGE_HEADER)
 return out

class TestSaigeEasJoin(unittest.TestCase):
 def test_saige_ref_alt_effect_allele2(self):
  status,out=saige_candidate(row())
  self.assertEqual(status,"CANDIDATE")
  self.assertEqual(out[:5],("1",10353,"A","G","G"))
 def test_effect_allele_can_be_marker_ref(self):
  r=row();r[3]="G";r[4]="A"
  status,out=saige_candidate(r)
  self.assertEqual(status,"CANDIDATE")
  self.assertEqual(out[4],"A")
 def test_exclude_source_qc_fail_even_with_valid_stats(self):
  r=row();r[19]="FAIL"
  self.assertEqual(saige_candidate(r)[0],"SOURCE_QC_FAIL")
 def test_exclude_unknown_qc(self):
  r=row();r[19]="."
  self.assertEqual(saige_candidate(r)[0],"SOURCE_QC_.")
 def test_exclude_palindrome(self):
  r=row();r[2]="1:10353:A:T";r[4]="T"
  self.assertEqual(saige_candidate(r)[0],"PALINDROMIC_UNRESOLVED")
 def test_malformed_markers(self):
  r=row();r[2]="1:10354:A:G"
  self.assertEqual(saige_candidate(r)[0],"MALFORMED_MARKER_ID")
 def test_non_snv(self):
  r=row();r[2]="1:10353:AT:G";r[3]="AT"
  self.assertEqual(saige_candidate(r)[0],"NOT_BIALLELIC_SNV")
 def test_bad_n(self):
  r=row();r[17]="9000"
  self.assertEqual(saige_candidate(r)[0],"BAD_NUMERIC_OR_N")
 def test_negative_p(self):
  r=row();r[12]="-0.01"
  self.assertEqual(saige_candidate(r)[0],"BAD_NUMERIC_OR_N")
 def test_coordinate_allele_join_and_freq_orientation(self):
  with tempfile.TemporaryDirectory() as d:
   conn=make_sqlite(Path(d)/"db.sqlite")
   conn.execute("INSERT INTO saige(chrom,pos,r,a,eff,freq,info,beta,se,pval) VALUES (?,?,?,?,?,?,?,?,?,?)",
      ("1",999,"A","G","G",.4,.99,.15,.02,.001))
   conn.execute("INSERT INTO outcome(chrom,pos,ea,oa,eaf,beta,se,pval) VALUES (?,?,?,?,?,?,?,?)",
      ("1",999,"A","G",.6,.17,.03,.002))
   stats,af,bychr=compare(conn)
   self.assertEqual(stats["MATCHED"],1)
   self.assertEqual(stats["INFO_GE_0.8_AND_SAIGE_MAF_GE_0.01"],1)
   self.assertEqual(stats["INFO_GE_0.8_MAF_GE_0.01_AF_DELTA_LE_0.05"],1)
   self.assertEqual(stats["INFO_GE_0.8_MAF_GE_0.01_AF_DELTA_LE_0.10"],1)
   self.assertEqual(af["DELTA_LE_0.05"],1)
   self.assertEqual(bychr["1"],1)
   conn.close()
 def test_joint_filter_excludes_large_af_discrepancy(self):
  with tempfile.TemporaryDirectory() as d:
   conn=make_sqlite(Path(d)/"db.sqlite")
   conn.execute("INSERT INTO saige(chrom,pos,r,a,eff,freq,info,beta,se,pval) VALUES (?,?,?,?,?,?,?,?,?,?)",
      ("1",999,"A","G","G",.4,.99,.15,.02,.001))
   conn.execute("INSERT INTO outcome(chrom,pos,ea,oa,eaf,beta,se,pval) VALUES (?,?,?,?,?,?,?,?)",
      ("1",999,"G","A",.9,.17,.03,.002))
   stats,af,_=compare(conn)
   self.assertEqual(stats["MATCHED"],1)
   self.assertEqual(stats.get("INFO_GE_0.8_MAF_GE_0.01_AF_DELTA_LE_0.05",0),0)
   self.assertEqual(stats["AF_DELTA_GT_0.1"],1)
   conn.close()
 def test_ambiguous_same_coordinate_is_not_counted(self):
  with tempfile.TemporaryDirectory() as d:
   conn=make_sqlite(Path(d)/"db.sqlite")
   conn.execute("INSERT INTO saige(chrom,pos,r,a,eff,freq,info,beta,se,pval,ambiguous) VALUES (?,?,?,?,?,?,?,?,?,?,1)",
      ("1",999,"A","G","G",.4,.99,.15,.02,.001))
   conn.execute("INSERT INTO outcome(chrom,pos,ea,oa,eaf,beta,se,pval) VALUES (?,?,?,?,?,?,?,?)",
      ("1",999,"G","A",.4,.17,.03,.002))
   self.assertEqual(compare(conn)[0].get("MATCHED",0),0)
   conn.close()
if __name__=="__main__":unittest.main()
