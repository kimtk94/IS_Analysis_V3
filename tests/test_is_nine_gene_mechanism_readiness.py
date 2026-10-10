"""Guard QTL H4 vs healthy control expression, NEURL legacy IDs and 2225 scope."""
import importlib.util,sys,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/is/build_is_nine_gene_mechanism_readiness.py"
spec=importlib.util.spec_from_file_location("nine_gene_ctx",P)
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m
spec.loader.exec_module(m)

def fixture():
    stable={gene:f"ENSG0000014{i:06d}" for i,gene in enumerate(m.NINE)}
    stable["NEURL1"]="ENSG00000107954"
    genes=[]
    for i,gene in enumerate(m.NINE):
        symbol="NEURL" if gene=="NEURL1" else gene
        genes.append({"gene_id_stable":stable[gene],
                      "canonical_GENCODE_v19_gene_symbol":symbol,
                      "n_positional_regions":"1",
                      "new_DIRECT_prior_grid_assays":"646" if gene=="ALDH2" else "0"})
    genes+=[{"gene_id_stable":f"ENSG0000090{i:06d}",
            "canonical_GENCODE_v19_gene_symbol":f"OTHER{i}",
            "n_positional_regions":"1","new_DIRECT_prior_grid_assays":"0"}
            for i in range(2216)]
    cell=[];features=[]
    for gene in m.NINE:
        present=gene not in ("FGF5","C4orf22")
        cell.append({"gene":gene,"present":"TRUE" if present else "FALSE",
                     "top_celltype":"Endothelial cells" if present else "",
                     "top_broad_class":"Endothelial" if present else "",
                     "top_detection_fraction":"0.21" if present else "NA",
                     "top_pseudobulk_cpm":"23.4" if present else "NA",
                     "total_detected_cells":"200" if present else "0"})
        features.append({"gene":gene,"feature_status":"FEATURE_PRESENT" if present else "FEATURE_NOT_IN_MATRIX"})
    grid=[]
    for i in range(646):
        for p in [1e-6,3e-6,1e-5,3e-5,1e-4]:
            grid.append({"gene_base":stable["ALDH2"],"locus":"BBJ_IS_L003",
                         "dataset_key":f"GTEx_V8__Brain_{i}","p12":str(p),
                         "PP_H4":"0.3","PP_H3":"0.4"})
    manifest={"donor_analysis_unit":"Patient",
              "reference_type":"ADULT_CONTROL_TEMPORAL_LOBE_VASCULAR_PERIVASCULAR"}
    return {"genes":genes,"cell":cell,"features":features,
            "donor":[],"grid":grid,"human_manifest":manifest}

class TestNineGeneReadiness(unittest.TestCase):
    def test_no_stroke_dge_and_missing_features_are_not_negative(self):
        records=m.audited(**fixture())
        self.assertEqual(len(records),9)
        self.assertTrue(all(x["stroke_patient_control_DGE"]=="NOT_PERFORMED" for x in records))
        self.assertTrue(all(x["candidate_causality"]=="NOT_ESTABLISHED" for x in records))
        self.assertEqual(sum(x["healthy_adult_reference_feature_status"]=="FEATURE_NOT_IN_MATRIX" for x in records),2)
        self.assertEqual(next(x for x in records if x["historic_symbol"]=="NEURL1")["current_GENCODE_v19_symbol"],"NEURL")
    def test_bad_feature_reference_fails(self):
        inp=fixture()
        inp["cell"][0]["present"]="TRUE"
        with self.assertRaises(ValueError):m.audited(**inp)
    def test_incorrect_stable_id_alias_fails(self):
        inp=fixture()
        neurl=next(x for x in inp["genes"] if x["canonical_GENCODE_v19_gene_symbol"]=="NEURL")
        neurl["gene_id_stable"]="ENSG00000000000"
        with self.assertRaises(ValueError):m.audited(**inp)
if __name__=="__main__":
    unittest.main()
