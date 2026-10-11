"""Strict local, no-network SNP coordinate/reference/allele overlap guards."""
import importlib.util,sys,tempfile,unittest
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas5_public_qtl_gwas_snp_overlap.py"
spec=importlib.util.spec_from_file_location("is_eas5_snp_audit",f)
a=importlib.util.module_from_spec(spec);sys.modules[spec.name]=a;spec.loader.exec_module(a)
class StubLifter:
 def __init__(self,chrom="chr1",position=19,strand="+",duplicate=False):
  self.chrom=chrom;self.position=position;self.strand=strand;self.duplicate=duplicate
 def convert_coordinate(self,c,p):
  result=[(self.chrom,self.position,self.strand,0)]
  return result+result if self.duplicate else result
class StubFasta:
 def __init__(self,allele):self.allele=allele
 def get(self,chrom,pos):return self.allele
class TestEAS5Overlap(unittest.TestCase):
 def test_positive_reference_and_qtl_effect_allele_invariants(self):
  cache={}
  status,variant,strand=a.canonical_qtl("1",20,"A","G",StubLifter(),StubFasta("A"),cache)
  self.assertEqual((status,variant,strand),("GRCH37_FASTA_REF_VERIFIED","1:20:A:G","+"))
  self.assertEqual(len(cache),1)
 def test_reverse_strand_complement_but_never_beta_sign_swapped(self):
  status,variant,strand=a.canonical_qtl("1",20,"T","C",
      StubLifter(strand="-"),StubFasta("A"),{})
  self.assertEqual((status,variant,strand),("GRCH37_FASTA_REF_VERIFIED","1:20:A:G","-"))
 def test_mismatched_ref_and_multimapping_are_unresolved(self):
  st,_,_=a.canonical_qtl("1",20,"A","G",StubLifter(),StubFasta("G"),{})
  self.assertEqual(st,"GRCH37_FASTA_REF_DISAGREES")
  st,_,_=a.canonical_qtl("1",20,"A","G",StubLifter(duplicate=True),StubFasta("A"),{})
  self.assertEqual(st,"LIFTOVER_NOT_UNIQUE")
 def test_gwas_ref_alt_orientation(self):
  row={"chr":"1","pos":"20","variant_id":"1:20:A:G",
       "effect_allele":"A","other_allele":"G","beta":"0.12",
       "se":"0.04","p":"0.003","eaf":"0.72"}
  status,record=a.checked_gwas_row(row,"EAS_AIS",StubFasta("A"))
  self.assertEqual(status,"VALID_GWAS_SOURCE_REF_ALT_EFFECT")
  self.assertAlmostEqual(record["beta_alt"],-.12)
  self.assertAlmostEqual(record["af_alt"],.28)
  row["effect_allele"]="G";row["other_allele"]="A"
  status,record=a.checked_gwas_row(row,"EAS_AIS",StubFasta("A"))
  self.assertAlmostEqual(record["beta_alt"],.12)
  self.assertAlmostEqual(record["af_alt"],.72)
 def test_gwas_ref_verification_fails_closed(self):
  row={"chr":"1","pos":"20","variant_id":"1:20:A:G",
       "effect_allele":"G","other_allele":"A","beta":"0.12",
       "se":"0.04","p":"0.003","eaf":"0.72"}
  status,_=a.checked_gwas_row(row,"EAS_AIS",StubFasta("G"))
  self.assertEqual(status,"GWAS_GRCH37_FASTA_REF_MISMATCH")
 def test_reverse_chain_build_header_gate(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"hg38ToHg19.over.chain"
   p.write_text("chain 100 chr1 248956422 + 100 200 chr1 249250621 + 100 200 1\n")
   self.assertTrue(a.validate_reverse_chain(p))
   p.write_text("chain 100 chr1 249250621 + 100 200 chr1 248956422 + 100 200 1\n")
   with self.assertRaises(ValueError):a.validate_reverse_chain(p)
   forward=Path(d)/"hg19ToHg38.over.chain"
   forward.write_text(p.read_text())
   with self.assertRaises(ValueError):a.validate_reverse_chain(forward)
 def test_fai_wrapped_sequence_random_access_and_cache(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"one.fa"
   seq="AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT"
   p.write_bytes((">1\n"+seq[:10]+"\n"+seq[10:20]+"\n"+seq[20:]+"\n").encode())
   # Build FAI for 10bp per physical row - all lines have exactly 10 bp.
   # Last row length 28 is not valid FASTA fixed wrapping: use 50bp exact.
   seq="ACGTACGTAC"*5
   p.write_bytes((">1\n"+"\n".join(seq[i:i+10] for i in range(0,50,10))+"\n").encode())
   (Path(str(p)+".fai")).write_text("1\t50\t3\t10\t11\n")
   fasta=a.Fasta37(p)
   for pos in (1,10,11,20,21,40,50,5,25,5):
    self.assertEqual(fasta.get("1",pos),seq[pos-1])
   self.assertIsNone(fasta.get("1",51))
   fasta.close()
if __name__=="__main__":unittest.main()
