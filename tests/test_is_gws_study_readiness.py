"""Study-level readiness: never silently promote allele overlap into fine-mapping."""
import csv,gzip,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from build_gws_study_readiness import process,FIELDS

def save(path,fields,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=fields)
        w.writeheader();w.writerows(rows)
class ReadinessTest(unittest.TestCase):
    def test_lead_status_and_no_finemap_promotion(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            groups=[];reports=[];variants=[]
            for i in range(7):
                tag=f"G{i}"
                groups.append(dict(group_id=tag,has_gws="1",chr="2",
                    lead_variant=f"2:{100+i}:A:G",phenotypes="AIS"))
                reports.append(dict(dataset="BBJ",group_id=tag,
                    window_variants="1",exported_matching_pairs="1"))
                row={k:"" for k in FIELDS}
                row.update(dict(dataset="BBJ",phenotype="AIS",group_id=tag,chr="2",
                    pos=str(100+i),variant_id=f"2:{100+i}:A:G",ref="A",alt="G",
                    effect_allele="G",other_allele="A",beta="0.1",se="0.02",
                    p="1e-9",eaf="0.25",alt_effect_beta="0.1",
                    alt_effect_eaf="0.25",qc_status="MATCH_ALT_EFFECT",
                    ancestry="EAS",build="GRCh37"))
                variants.append(row)
            save(root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv",list(groups[0]),groups)
            save(root/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv",list(reports[0]),reports)
            with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","wt") as f:
                w=csv.DictWriter(f,delimiter="\t",fieldnames=FIELDS)
                w.writeheader();w.writerows(variants)
            result=process(root)
            self.assertEqual(result["with_oriented_lead"],7)
            self.assertEqual(result["ready_for_multi_signal_finemap"],0)
            with (root/"IS_GWS_STUDY_INPUT_READINESS.tsv").open() as f:
                records=list(csv.DictReader(f,delimiter="\t"))
            self.assertTrue(all(r["multi_signal_finemap_ready"]=="NO_LD_MATRIX_OR_SIGNAL_QC" for r in records))
if __name__=="__main__":unittest.main()
