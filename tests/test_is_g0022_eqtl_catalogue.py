"""eQTL Catalogue 19-column nominal QTL input and exact allele liftover."""
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_g0022_eqtl_catalogue_aorta_coloc import (
  EQTLCATALOGUE_COLUMNS,parse_nominal,read_qtl,numeric)

class StubLift:
    def convert_coordinate(self,chrom,pos):
        return [("chr12",200,"+",800)] if chrom=="chr12" and pos==100 else []
def qtlrow(geneid="ENSG00000111275",ref="A",alt="G",pos="101",beta=".22"):
    d={k:"NA" for k in EQTLCATALOGUE_COLUMNS}
    d.update({"molecular_trait_id":geneid,"chromosome":"12","position":pos,
      "ref":ref,"alt":alt,"variant":f"chr12_{pos}_{ref}_{alt}",
      "ma_samples":"200","maf":".35","pvalue":"1e-7","beta":beta,
      "se":".06","type":"SNP","ac":"250","an":"774",
      "molecular_trait_object_id":geneid,"gene_id":geneid,
      "median_tpm":"10","rsid":"rs1"})
    return "\t".join(d[x] for x in EQTLCATALOGUE_COLUMNS)
class QTLPreparationTests(unittest.TestCase):
    def test_schema_19_and_malformed(self):
        self.assertEqual(len(EQTLCATALOGUE_COLUMNS),19)
        self.assertEqual(parse_nominal(qtlrow())["maf"],".35")
        with self.assertRaisesRegex(ValueError,"19-column"):
            parse_nominal("a\tb")
        self.assertIsNone(numeric("NA"))
    def test_exact_ref_alt_duplicate_rsid_and_no_misleading_flip(self):
        import gzip
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/"a.gz"
            with gzip.open(path,"wt") as f:
                f.write(qtlrow()+"\n")
                f.write(qtlrow()+"\n") # same variant, different rsid would be equivalent
                f.write(qtlrow(geneid="ENSG9999")+"\n")
            src={"12:201:A:G":{"variant_id":"12:201:A:G","ref":"A",
                    "alt":"G","alt_effect_beta":".2","se":".03","p":"1e-8",
                    "alt_effect_eaf":".35"}}
            out,qc=read_qtl(path,StubLift(),{"ENSG00000111275":"ALDH2"},src)
            self.assertEqual(qc["matched_unique_gene_snp"],1)
            self.assertEqual(qc["same_variant_multiple_rsid_duplicate"],1)
            self.assertEqual(len(out["ENSG00000111275"]),1)
            self.assertEqual(out["ENSG00000111275"]["12:201:A:G"]["allele_match"],
                "EXACT_REF_ALT_AFTER_LIFTOVER")
    def test_wrong_allele_does_not_match(self):
        import gzip
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/"a.gz"
            with gzip.open(path,"wt") as f:f.write(qtlrow(ref="A",alt="C")+"\n")
            src={"12:201:A:G":{}}
            out,qc=read_qtl(path,StubLift(),{"ENSG00000111275":"ALDH2"},src)
            self.assertEqual(qc["no_gwas_allele_pair"],1)
            self.assertEqual(dict(out),{})
if __name__=="__main__":unittest.main()
