"""Non-GTEx cross-ancestry QTL source candidate integrity and SNP liftover."""
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_g0022_non_gtex_qtl_sources import build
from audit_g0022_non_gtex_qtl_variant_coverage import mapping
class DummyChain:
    def __init__(self,reverse=False):self.reverse=reverse
    def convert_coordinate(self,chrom,pos):
        if self.reverse:
            return [("chr12",99,"+",1)] if pos==199 else []
        return [("chr12",199,"+",1)] if pos==99 else []
class SourceTests(unittest.TestCase):
    def test_six_exact_accessions_and_no_ancestry_assumption(self):
        with tempfile.TemporaryDirectory() as td:
            folder=Path(td);out=folder/"out"
            src=folder/"catalogue.tsv"
            records=[]
            ids=("QTD000609","QTD000620","QTD000434",
                 "QTD000051","QTD000021","QTD000110")
            for i,ds in enumerate(ids):
                records.append({"dataset_id":ds,"study_id":f"QTS{i}",
                    "study_label":f"Study{i}","sample_group":"blood",
                    "sample_size":"300","quant_method":"ge",
                    "ftp_path":f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS{i}/{ds}/{ds}.all.tsv.gz"})
            with src.open("w") as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
                w.writeheader();w.writerows(records)
            rows=build(src,out)
            self.assertEqual(len(rows),6)
            self.assertEqual(set(x["dataset_id"] for x in rows),set(ids))
            self.assertTrue(all(x["ancestry_verified"]=="NO" for x in rows))
            records[0]["quant_method"]="tx"
            with src.open("w") as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
                w.writeheader();w.writerows(records)
            with self.assertRaisesRegex(ValueError,"Not gene-level"):
                build(src,out)
    def test_liftover_checked_reciprocal(self):
        input=[{"variant_id":"12:100:A:G","z":"2","ref":"A","alt":"G"}]
        r=mapping(input,DummyChain(),DummyChain(True))
        self.assertEqual(r[(200,"A","G")]["variant_id"],"12:100:A:G")
        class BadReverse:
            def convert_coordinate(self,chrom,pos):return []
        self.assertEqual(mapping(input,DummyChain(),BadReverse()),{})
if __name__=="__main__":unittest.main()
