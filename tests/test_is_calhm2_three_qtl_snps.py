import csv,tempfile,unittest
from pathlib import Path
from scripts.is.audit_is_calhm2_three_qtl_snps import create,FIELDS,POPS,LEADS
class TestCALHM2ThreeQTL(unittest.TestCase):
 def test_exact_three_and_lead_tagging(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);ev=root/"evidence.tsv";ld=root/"ld.tsv";out=root/"result"
   hdr=["locus","symbol","BBJ_exact_Fasta37","GIGA_EAS_AIS_exact_Fasta37","palindromic_review",*FIELDS]
   with ev.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=hdr,delimiter="\t");w.writeheader()
    for i in range(3):
     r=dict(locus="BBJ_IS_L002",symbol="CALHM2",BBJ_exact_Fasta37="True",
      GIGA_EAS_AIS_exact_Fasta37="True",palindromic_review="False",
      SNP_GRCh37=f"10:{105400000+i}:A:G",SNP_GRCh38=f"10:{105500000+i}:A:G",
      QTD000161_nominal_p="0.00001",QTD000161_beta="0.2",QTD000161_se="0.05",
      BBJ_p="0.003",BBJ_beta_ALT="0.04",GIGA_EAS_AIS_p="0.02",
      GIGA_EAS_AIS_beta_ALT="0.03")
     w.writerow(r)
   with ld.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["population","lead","variant_b37","signed_r","r2","r2_ge_0p8"],delimiter="\t");w.writeheader()
    for pop in POPS:
     for lead in LEADS:
      for i in range(3):
       w.writerow(dict(population=pop,lead=lead,variant_b37=f"10:{105400000+i}:A:G",
        signed_r=".95" if lead=="BBJ" else=".1",r2=".9025" if lead=="BBJ" else=".01",
        r2_ge_0p8=str(lead=="BBJ")))
   result=create(ev,ld,out)
   self.assertEqual(len(result),3)
   self.assertTrue(all(r["JPT104_BBJ_tag_r2_ge_0p8"]=="True" for r in result))
   self.assertTrue(all(r["EAS504_GIGASTROKE_tag_r2_ge_0p8"]=="False" for r in result))
   self.assertTrue((out/"CALHM2_THREE_QTL_SNP_MANIFEST.json").exists())
 def test_empty_or_unexpected_snp_count_rejected(self):
  with tempfile.TemporaryDirectory() as t:
   ev=Path(t)/"e.tsv";ld=Path(t)/"ld.tsv"
   ev.write_text("locus\tsymbol\tSNP_GRCh37\n")
   ld.write_text("population\tlead\tvariant_b37\n")
   with self.assertRaises(ValueError):create(ev,ld,Path(t)/"o")
