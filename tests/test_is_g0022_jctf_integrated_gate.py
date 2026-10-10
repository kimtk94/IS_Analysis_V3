"""Integrated G0022 QTL evidence never silently promotes marginal QTL to coloc."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_jctf_integrated_evidence import build
def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj))
def tab(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t")
        w.writeheader();w.writerows(data)
class IntegratedGateTest(unittest.TestCase):
    def fixtures(self,base):
        alt=base/"alternate_qtl_sources_v1";job=base/"jctf_japan_omics_variants_v1"
        gt=base/"eqtl_catalogue_multitissue_v1"
        save(base/"G0022_FULL_AIS_PUBLICATION_GATE_SUMMARY.json",{
          "causal_gene_publication_ready":False,"causal_signal_publication_ready":False})
        save(alt/"G0022_NON_GTEX_QTL_CS_COVERAGE_SUMMARY.json",{
          "cs_variant_dataset_tests":24,"source_test_status":{"POSITION_ABSENT":24}})
        save(gt/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE_SUMMARY.json",{
          "by_coverage_status":{"NOT_TESTED_VARIANT_POSITION_ABSENT":16}})
        save(base/"G0022_AIS_RS671_CONDITIONAL_DIAGNOSTIC_SUMMARY.json",{
          "conditional_gws_p_lt_5e8":0,"reference_n":504})
        save(job/"G0022_JCTF_4SNP_QTL_EVIDENCE_SUMMARY.json",{
          "ALDH2_eQTL_variant_records":4,"ALDH2_pQTL_variant_records":4,
          "colocalization_numerically_calculated":False,"validated_causal_stroke_genes":0})
        variants=["12:112241766:G:A","12:112168009:G:A",
                  "12:112468206:C:T","12:112736118:A:G"]
        tab(base/"G0022_AIS_FOUR_CS_RS671_LD.tsv",
            [dict(variant_id=v,annotation="ALDH2 - positional",
                  r2_to_rs671=1.0) for v in variants])
        rs=[]
        for v in variants:
            for x in ("eQTL","pQTL"):
                rs.append(dict(GWAS_variant_GRCh37=v,gene_symbol="ALDH2",
                     qtl_category=x,gwas_ais_full_locus_pip=.25,
                     jctf_allele_association_p=1e-8,jctf_susie_variant_pip=.3))
        tab(job/"G0022_JCTF_4SNP_JAPANESE_EQTL_PQTL_EVIDENCE.tsv",rs)
        return base,alt,job,gt
    def test_fully_present_jctf_evidence_is_not_causal_coloc(self):
        with tempfile.TemporaryDirectory() as td:
            base,alt,job,gt=self.fixtures(Path(td))
            result=build(base,alt,job,gt)
            self.assertFalse(result["molecular_colocalization_complete"])
            self.assertFalse(result["causal_ALDH2_stroke_mechanism_resolved"])
            self.assertEqual(result["East_Asian_JCTF_QTL_variants_with_ALDH2_eQTL"],4)
            x=json.loads((job/"G0022_JCTF_4SNP_QTL_EVIDENCE_SUMMARY.json").read_text())
            x["validated_causal_stroke_genes"]=1
            save(job/"G0022_JCTF_4SNP_QTL_EVIDENCE_SUMMARY.json",x)
            with self.assertRaisesRegex(RuntimeError,"causal"):
                build(base,alt,job,gt)
if __name__=="__main__":unittest.main()
