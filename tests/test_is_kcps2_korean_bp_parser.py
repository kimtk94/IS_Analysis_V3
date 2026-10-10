"""KCPS2 Korean SBP/DBP rs671 original-source allele and integrity tests."""
import csv,gzip,hashlib,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from extract_kcps2_korean_rs671_bp import extract
from extract_kcps2_korean_rs671_alcohol import EXPECTED
def src(root,trait,allele2="A",pos=112241766):
    name=f"SAIGE_{trait.upper()}_INFO.txt.gz"
    r={key:"0" for key in EXPECTED}
    r.update(CHR="12",POS=str(pos),MarkerID=f"12:{pos}:G:A",
        Allele1="G" if allele2=="A" else "A",Allele2=allele2,
        AC_Allele2="30000",AF_Allele2=".155",
        imputationInfo=".98",BETA="-.08" if allele2=="A" else ".08",
        SE=".007",Tstat="-11",var="1",**{"p.value":"1e-20"},
        N="130000",INFO=".98",neglogP="20",Z="-11")
    p=root/name
    with gzip.open(p,"wt") as f:
        w=csv.DictWriter(f,fieldnames=list(EXPECTED),delimiter="\t")
        w.writeheader();w.writerow(r)
    (root/(name+".member_crc_verified.json")).write_text(json.dumps({
        "archive_member_name":name,
        "status":"SELECTIVE_ZIP_MEMBER_CRC_VERIFIED_ARCHIVE_MD5_NOT_COMPUTED",
        "ZIP_member_size_bytes":p.stat().st_size,
        "inner_gzip_sha256":hashlib.sha256(p.read_bytes()).hexdigest()}))
class BP(unittest.TestCase):
    def test_sbp_harmonization(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);src(root,"sbp")
            r=extract(root,root/"out","sbp",minrows=1)
            self.assertEqual(r["rs671_effect_A_beta"],-.08)
            self.assertFalse(r["mmHg_interpretation_permitted"])
            self.assertFalse(r["entire_original_ZIP_MD5_verified"])
    def test_dbp_ref_effect_flip(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);src(root,"dbp",allele2="G")
            r=extract(root,root/"out","dbp",1)
            self.assertEqual(r["rs671_effect_A_beta"],-.08)
    def test_fail_changed_zip_inner_data(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);src(root,"sbp")
            with (root/"SAIGE_SBP_INFO.txt.gz").open("ab") as f:f.write(b"x")
            with self.assertRaisesRegex(ValueError,"integrity"):
                extract(root,root/"out","sbp",1)
if __name__=="__main__":unittest.main()
