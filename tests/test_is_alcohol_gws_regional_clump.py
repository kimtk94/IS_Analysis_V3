"""End-to-end synthetic regional GWAS association-to-1KG REF/ALT clump input."""
import csv,gzip,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_alcohol_gws_regional_clump import make,same_pair,is_palindrome
from extract_alcohol_eas_region_reference_for_clump import SENTINELS
ALLELES={("2",27730940):("T","C"),("4",39413780):("A","G"),
 ("4",100239319):("T","C"),("9",38395928):("T","C"),
 ("9",75461066):("T","C"),("12",106750302):("A","G"),
 ("12",112241766):("G","A")}
HDR=["SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","HetP","N"]
class RegionClumpTest(unittest.TestCase):
    def fixture(self,base):
        ref=base/"ref";ref.mkdir()
        for ch,poss in SENTINELS.items():
            (ref/f"chr{ch}.region_source_audit.json").write_text(json.dumps({
                "EAS_samples":504,"plink_sample_count":504,"validated":True}))
            with (ref/f"chr{ch}.EAS504.regions.pvar").open("w") as f:
                f.write('##meta=VCF4\n#CHROM\tPOS\tID\tREF\tALT\n')
                for k in range(110):
                    f.write(f"{ch}\t{200000000+k}\t{ch}:{200000000+k}:G:A\tG\tA\n")
                for pos in poss:
                    a,b=ALLELES[(ch,pos)]
                    if (ch,pos)==("4",39413780):a,b=b,a
                    f.write(f"{ch}\t{pos}\t{ch}:{pos}:{a}:{b}\t{a}\t{b}\n")
        path=base/"gwas.gz"
        with gzip.open(path,"wt") as f:
            w=csv.DictWriter(f,fieldnames=HDR,delimiter="\t")
            w.writeheader()
            for (ch,pos),(a,b) in ALLELES.items():
                w.writerow(dict(SNP=f"chr{ch}_{pos}_{a}_{b}",CHR=ch,POS=pos,
                    EA=b,NEA=a,EAF=".4",BETA="-.4",SE=".04",
                    P="0" if pos==112241766 else "1e-10",HetP=".2",N="154570"))
        return ref,path
    def test_all_gws_sentinel_reference_ids_present_and_clamped_zero(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);ref,src=self.fixture(base)
            r=make(base/"out",ref,src,source_min_rows=7)
            self.assertEqual(r["source_original_total_rows"],7)
            self.assertEqual(r["p_zero_clamped_for_PLINK"],1)
            self.assertEqual(sum(x["clump_input_rows_after_dedup"] for x in r["chrom_source_stats"].values()),7)
            with (base/"out"/"chr12.alcohol_regional.clump_input.tsv").open() as f:
                vals=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(len(vals),2)
            self.assertEqual(min(float(x["P"]) for x in vals),1e-300)
    def test_harmonization_and_palindrome_gates(self):
        self.assertTrue(same_pair("A","G","G","A"))
        self.assertFalse(same_pair("A","G","T","A"))
        self.assertTrue(is_palindrome("A","T"))
        self.assertFalse(is_palindrome("A","G"))
if __name__=="__main__":unittest.main()
