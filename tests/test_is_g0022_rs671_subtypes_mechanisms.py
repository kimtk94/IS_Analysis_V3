"""G0022 rs671 subtypes, molecular gene diversity and mediation blockers."""
import csv,gzip,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_rs671_stroke_subtypes import run as stroke,STUDIES,CS
from audit_g0022_rs671_mediation_readiness import run as mediation
from audit_g0022_jctf_all_gene_qtl import build as genes

def put(path,rows,compress=False):
    path.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if compress else open
    with op(path,"wt") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
class NextStageTests(unittest.TestCase):
    def fixture_source(self):
        out=[]
        for ds in STUDIES:
            for var in CS:
                chrom,pos,ref,alt=var.split(":")
                reverse=(ds=="GCST90104548")
                out.append({"group_id":"IS_XDATA_G0022","dataset":ds,
                    "variant_id":var,"chr":chrom,"pos":pos,"ref":ref,"alt":alt,
                    "build":"GRCh37","ancestry":"Japanese" if ds=="BBJ" else "EAS",
                    "qc_status":"MATCH_REF_EFFECT" if reverse else "MATCH_ALT_EFFECT",
                    "effect_allele":ref if reverse else alt,
                    "beta":".10" if reverse else "-.10",
                    "alt_effect_beta":"-.10","eaf":".8" if reverse else ".2",
                    "alt_effect_eaf":".2","se":".03","p":".001"})
        return out
    def test_24_allele_oriented_rows_and_flipped_study(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td);source=t/"a.gz";put(source,self.fixture_source(),True)
            r=stroke(source,t/"out")
            self.assertEqual(r["rows"],24)
            self.assertTrue(r["all_24_study_variant_effects_protective"])
            self.assertFalse(r["phenotype_independence_verified"])
    def test_wrong_eaf_or_missing_study_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td);data=self.fixture_source();data[0]["alt_effect_eaf"]=".9"
            source=t/"a.gz";put(source,data,True)
            with self.assertRaisesRegex(ValueError,"beta/EAF"):
                stroke(source,t/"out")
            put(source,data[:-1],True)
            with self.assertRaisesRegex(ValueError,"Incomplete"):
                stroke(source,t/"out")
    def test_mediation_gates_never_promote_single_instrument(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td)
            (t/"G0022_RS671_STROKE_SUBTYPE_SUMMARY.json").write_text(json.dumps({
                "rows":24,"all_24_study_variant_effects_protective":True,
                "rs671_AIS_OR":.86,"rs671_phenotypes":6}))
            (t/"G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json").write_text(json.dumps({
                "molecular_colocalization_complete":False,
                "causal_ALDH2_stroke_mechanism_resolved":False,
                "East_Asian_JCTF_QTL_variants_with_ALDH2_eQTL":4,
                "East_Asian_JCTF_QTL_variants_with_ALDH2_pQTL":4}))
            (t/"G0022_AIS_RS671_CONDITIONAL_DIAGNOSTIC_SUMMARY.json").write_text(json.dumps({
                "reference_n":504}))
            r=mediation(t)
            self.assertFalse(r["eligible_single_rs671_MR_Egger"])
            self.assertFalse(r["eligible_multivariable_mediation"])
            self.assertEqual(r["valid_mediation_effects_computed"],0)
    def test_jctf_all_gene_evidence_not_narrowed(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td);p=t/"input.tsv"
            data=[]
            for i in range(72):
                data.append({"gene_id":f"ENSG{i%19:011d}","gene_symbol":f"G{i%19}",
                    "qtl_category":"eQTL" if i%19!=18 else "pQTL",
                    "jctf_allele_association_p":".005",
                    "jctf_susie_variant_pip":".1",
                    "jctf_finemap_variant_pip":".05",
                    "GWAS_variant_GRCh37":f"12:{i}:A:G"})
            put(p,data)
            r=genes(p,t/"out")
            self.assertEqual(r["source_rows"],72)
            self.assertEqual(r["unique_genes"],19)
            self.assertEqual(r["validated_AIS_gene_mediation"],0)
if __name__=="__main__":unittest.main()
