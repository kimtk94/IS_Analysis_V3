"""G0022 source-verified four-tissue eQTL harmonization and assay coverage gates."""
import csv,gzip,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from prepare_g0022_eqtl_catalogue_multitissue import process_study,SOURCES
from prepare_g0022_eqtl_catalogue_aorta_coloc import EQTLCATALOGUE_COLUMNS
from audit_g0022_qtl_cs_assay_coverage import orient_alleles

def make_qtl(pos,ref,alt,gene,beta=".1"):
    data={col:"NA" for col in EQTLCATALOGUE_COLUMNS}
    data.update(molecular_trait_id=gene,chromosome="12",position=str(pos),
        ref=ref,alt=alt,variant=f"chr12_{pos}_{ref}_{alt}",
        ma_samples="100",maf=".3",pvalue="1e-5",beta=beta,
        se=".04",type="SNP",ac="120",an="400",
        molecular_trait_object_id=gene,gene_id=gene)
    return "\t".join(data[c] for c in EQTLCATALOGUE_COLUMNS)+"\n"
class DummyLift:
    def convert_coordinate(self,c,pos):
        return [("chr12",200,"+",987)] if c=="chr12" and pos==100 else []
class MultiQTLTests(unittest.TestCase):
    def test_forward_negative_strand_effect_alleles(self):
        self.assertEqual(orient_alleles("12:201:A:G","-"),("T","C"))
        self.assertEqual(orient_alleles("12:201:A:G","+"),("A","G"))
    def test_tissue_source_must_match_target_and_not_falsely_confirm_causality(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"nominal.gz"
            with gzip.open(p,"wt") as f:
                f.write(make_qtl(101,"A","G","ENSG00000111275"))
                f.write(make_qtl(101,"A","G","ENSG00000111275"))
                f.write(make_qtl(101,"A","C","ENSG00000111275"))
            g={"12:201:A:G":{"ref":"A","alt":"G","alt_effect_beta":".2",
                 "se":".03","p":"1e-8","alt_effect_eaf":".3"}}
            out=Path(td)/"results"
            result=process_study("QTD000136",p,
                {"ENSG00000111275":"ALDH2"},g,{"12:201:A:G"},
                DummyLift(),out,min_source_bytes=0)
            self.assertEqual(result["target_gene_snp_matches"],1)
            self.assertEqual(result["qc"]["duplicate_same_gene_snp"],1)
            self.assertEqual(result["GWAS_CS_SNPs_matched_in_8_target_genes"],
                ["12:201:A:G"])
            # Test that matching variant is still not a formal colocalization.
            self.assertEqual(result["status"],
                "ALLELE_HARMONIZATION_COMPLETE_COLOC_UNVALIDATED")
            with (out/"G0022_MULTITISSUE_QTL_INPUT_AUDIT.tsv").open() as f:
                audit=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(audit[0]["status"],
                "MISSING_GWS_OR_CS_AND_OR_QTL_SIGNAL_NOT_NEGATIVE")
    def test_gene_variant_absent_not_negative(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"nominal.gz"
            with gzip.open(p,"wt") as f:
                f.write(make_qtl(101,"A","C","ENSG00000111275"))
            result=process_study("QTD000171",p,
                {"ENSG00000111275":"ALDH2"},
                {"12:201:A:G":{"ref":"A","alt":"G","alt_effect_beta":".2",
                  "se":".03","p":"1e-8","alt_effect_eaf":".3"}},
                {"12:201:A:G"},DummyLift(),Path(td)/"out",min_source_bytes=0)
            self.assertEqual(result["target_gene_snp_matches"],0)
            self.assertEqual(result["GWAS_CS_SNPs_in_QTL_extract_any_gene"],[])
            self.assertEqual(result["genes_ready_for_diagnostic"],0)
if __name__=="__main__":unittest.main()
