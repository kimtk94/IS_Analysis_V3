"""Synthetic guardrails: 7-region PLINK2 clump MAF sensitivity audit."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_alcohol_eas_regional_clump_sensitivity import run,REGIONS,file_for
def save(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
class RegionalClumpAudits(unittest.TestCase):
    def test_7_region_correct_sentinel_and_no_causal_promotion(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/"audit";root.mkdir()
            ref=Path(td)/"ref";ref.mkdir()
            src=[];stats={}
            for c,v in REGIONS.items():
                srcchr=[];af=[]
                for name,pos in v:
                    variant=f"{c}:{pos}:A:G"
                    srcchr.append(dict(ID=variant,chr=c,pos=pos,
                       P_original="0" if name=="ALDH2" else "1e-10",
                       P_clamped_from_zero="1" if name=="ALDH2" else "0"))
                    af.append({"ID":variant,"ALT_FREQS":.3})
                src.extend(srcchr)
                stats[c]={"clump_input_GWS":len(srcchr),"source_p_lt_5e8":len(srcchr)}
                (ref/f"chr{c}.region_source_audit.json").write_text(json.dumps({
                   "validated":True,"plink_sample_count":504}))
                save(ref/f"chr{c}.EAS504.regions.afreq",af)
                for maf in (0,.01,.05):
                    entries=[{"ID":r["ID"],"POS":r["pos"],
                       "P":"1e-300" if r["P_clamped_from_zero"]=="1" else "1e-10",
                       "TOTAL":"0","SP2":"."} for r in srcchr]
                    save(root/file_for(c,maf),entries)
                if c in ("4","12"):
                    for r2 in (.2,.5):save(root/file_for(c,.01,r2),[{
                        "ID":srcchr[0]["ID"],"POS":srcchr[0]["pos"],
                        "P":"1e-10","TOTAL":"0","SP2":"."}])
            save(root/"IS_ALCOHOL_SOURCE_REF_MATCHED_REGIONAL_ASSOCIATIONS.tsv",src)
            (root/"IS_ALCOHOL_REGIONAL_GWAS_REFERENCE_COVERAGE.json").write_text(json.dumps({
                "source_original_total_rows":7676853,
                "p1_index_threshold":5e-8,"p2_secondary_threshold":1e-4,
                "chrom_source_stats":stats,"p_zero_clamped_for_PLINK":1}))
            r=run(root,ref)
            self.assertEqual(r["base_clump_count_MAF_unrestricted"],7)
            self.assertEqual(r["sensitivity_clump_count_MAF_ge_0p05"],7)
            self.assertEqual(r["sentinel_to_index_or_member"]["ALDH2"]["source_p_numeric_underflow"],True)
            self.assertFalse(r["one_clump_one_causal_signal_assumption_justified"])
            self.assertEqual(r["causal_eligible_MR_IVs"],0)
            (root/"chr12.EAS504.alcohol_maf_0.05_r2_0p1.clumps").unlink()
            with self.assertRaises(FileNotFoundError):
                run(root,ref)
if __name__=="__main__":unittest.main()
