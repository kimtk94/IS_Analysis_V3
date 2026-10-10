"""Synthetic source-stamped KCPS2 original SAIGE ALCO_AMOUNT rs671 parser."""
import csv,gzip,hashlib,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from extract_kcps2_korean_rs671_alcohol import extract,EXPECTED,NAME

class KoreanSource(unittest.TestCase):
    def create(self,root,allele1="G",allele2="A",pos=112241766):
        x={k:"0" for k in EXPECTED}
        x.update(CHR="12",POS=str(pos),MarkerID="12:112241766:G:A",
           Allele1=allele1,Allele2=allele2,
           AC_Allele2="400",AF_Allele2=".225",
           imputationInfo=".98",BETA="-.59" if allele2=="A" else ".59",
           SE=".03",Tstat="-20",var="1",**{"p.value":"1e-180"},
           N="150000",INFO=".98",neglogP="180",Z="-19.7")
        path=root/NAME
        with gzip.open(path,"wt") as f:
            w=csv.DictWriter(f,delimiter="\t",fieldnames=list(EXPECTED))
            w.writeheader();w.writerow(x)
        h=hashlib.sha256(path.read_bytes()).hexdigest()
        (root/(NAME+".member_crc_verified.json")).write_text(json.dumps({
          "status":"SELECTIVE_ZIP_MEMBER_CRC_VERIFIED_ARCHIVE_MD5_NOT_COMPUTED",
          "ZIP_member_size_bytes":path.stat().st_size,
          "inner_gzip_sha256":h}))
        return path
    def test_rs671_korean_replication_when_allele2_is_A(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.create(root)
            x=extract(root,root/"out",minimum_rows=1)
            self.assertEqual(x["rs671_source_build_at_position"],"GRCh37_hg19")
            self.assertAlmostEqual(x["rs671_ALCO_AMOUNT_beta_A"],-.59)
            self.assertFalse(x["full_KCPS2_archive_MD5_verified"])
            self.assertFalse(x["causal_mediation_estimated"])
    def test_ref_effect_flips_sign(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.create(root,allele1="A",allele2="G")
            x=extract(root,root/"out",1)
            self.assertAlmostEqual(x["rs671_ALCO_AMOUNT_beta_A"],-.59)
    def test_source_integrity_failure_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);path=self.create(root)
            with path.open("ab") as f:f.write(b"x")
            with self.assertRaisesRegex(ValueError,"incomplete"):
                extract(root,root/"out",1)
    def test_rs671_grch38_candidate_position(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.create(root,pos=111803962)
            x=extract(root,root/"out",1)
            self.assertEqual(x["rs671_source_build_at_position"],"GRCh38_hg38")
if __name__=="__main__":unittest.main()
