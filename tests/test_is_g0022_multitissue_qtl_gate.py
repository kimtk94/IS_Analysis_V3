"""Fail-closed multi-tissue publication gate: missing CS coverage != eQTL negative."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_multitissue_qtl_gate import audit
from prepare_g0022_eqtl_catalogue_multitissue import SOURCES
def write(path,items):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]),delimiter="\t")
        w.writeheader();w.writerows(items)
def dumps(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value))
class PublicationGate(unittest.TestCase):
    def fixtures(self,root):
        dumps(root/"G0022_MULTITISSUE_ABF_SUMMARY.json",{
            "analyzed_studies":4,"gene_tissue_tests":32,
            "computed_posteriors":24,
            "posterior_results_eligible_for_publication":0})
        dumps(root/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE_SUMMARY.json",{
            "dataset_variant_tests":16,
            "by_coverage_status":{"NOT_TESTED_VARIANT_POSITION_ABSENT":16}})
        data=[]
        for key,(tissue,n) in SOURCES.items():
            dumps(root/key/"G0022_MULTITISSUE_QTL_INPUT_SUMMARY.json",{
              "dataset_id":key,"catalogue_sample_size":n,
              "GWAS_CS_SNPs_in_QTL_extract_any_gene":[],
              "genes_ready_for_diagnostic":0,
              "region_grch38":"12:111650000-112370000",
              "qc":{"source_records":1500},"raw_sha256":"a"*64,
              "target_gene_snp_matches":800})
            for gene in ("ACAD10","ALDH2","NAA25","HECTD4","BRAP","PTPN11","IFT81","ATP2A2"):
                data.append({"dataset_id":key,"gene":gene,
                   "matched_snps":100,"PP.H4":".99",
                   "gate":"BLOCKED_NO_SHARED_GENOMEWIDE_SIGNIFICANT_AIS_VARIANT",
                   "conclusion":"NOT_VALIDATED_COLOCALIZATION"})
        write(root/"G0022_GTEX_MULTITISSUE_ABF_DIAGNOSTIC.tsv",data)
        return data
    def test_even_pph4_099_never_promoted_without_cs(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);rows=self.fixtures(root)
            verdict=audit(root)
            self.assertEqual(verdict["validated_causal_genes"],0)
            self.assertEqual(verdict["gtEx_gene_tissue_comparisons"],32)
            rows[0]["conclusion"]="VALIDATED_COLOCALIZATION"
            write(root/"G0022_GTEX_MULTITISSUE_ABF_DIAGNOSTIC.tsv",rows)
            with self.assertRaisesRegex(RuntimeError,"causal colocalization"):
                audit(root)
    def test_missing_source_coverage_denominator_fails(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);self.fixtures(root)
            dumps(root/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE_SUMMARY.json",{
              "dataset_variant_tests":15,
              "by_coverage_status":{"NOT_TESTED_VARIANT_POSITION_ABSENT":15}})
            with self.assertRaisesRegex(RuntimeError,"source.*coverage|Source"):
                audit(root)
if __name__=="__main__":unittest.main()
