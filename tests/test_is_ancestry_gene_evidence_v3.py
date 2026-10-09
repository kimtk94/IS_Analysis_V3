"""Strict source/locus boundary between EAS legacy coloc and EUR candidate discovery."""
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from build_ancestry_gene_evidence_v3 import build

def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(data[0]))
        w.writeheader();w.writerows(data)

class EvidenceV3(unittest.TestCase):
    def test_no_legacy_coloc_leak_into_eur(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)/"base";out=Path(td)/"v2";base.mkdir();out.mkdir()
            manifest=[{"group_id":f"IS_XDATA_G{i:04d}","ancestry":"EAS_AND_JAPANESE",
                "source":"BBJ_GIGASTROKE_EAS","has_gws":"1","lead_variant":"4:100:A:T",
                "lead_p":"1e-8"} for i in range(30)]
            manifest.append({"group_id":"IS_EUR_AIS_G0001","ancestry":"EUR",
                "source":"GCST90104540","has_gws":"1",
                "lead_variant":"4:100:A:T","lead_p":"1e-8"})
            save(out/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv",manifest)
            genes=[]
            for g in ("IS_XDATA_G0000","IS_EUR_AIS_G0001"):
                genes.append({"group_id":g,"gene_id":"ENSG123.3","gene_symbol":"GENE1",
                  "biotype":"protein_coding","lead_distance_bp":"100"})
            save(out/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv",genes)
            original=[{"group_id":"IS_XDATA_G0000","gene_id":"ENSG123.3",
              "coloc_status":"LOCUS_AND_GENE_ID_MATCHED","coloc_max_pp_h4":"0.95",
              "coloc_original_locus":"BBJ_IS_L001"}]
            save(base/"IS_LOCUS_GENE_EVIDENCE_V2.tsv",original)
            save(base/"LEGACY_ANCHOR_CANDIDATES.tsv",[{"gene":"GENE1"},{"gene":"NEURL1"}])
            info=build(out,base)
            self.assertEqual(info["candidate_evidence_rows"],3)
            self.assertEqual(info["legacy_coloc_region_pairs"],1)
            self.assertEqual(info["eur_region_pairs_inappropriately_reusing_coloc"],0)
            with (out/"IS_ANCESTRY_GENE_EVIDENCE_V3.tsv").open() as f:
                rows=list(csv.DictReader(f,delimiter="\t"))
            a,b,c=rows
            self.assertEqual(a["legacy_coloc_pp_h4"],"0.95")
            self.assertEqual(b["legacy_coloc_pp_h4"],"")
            self.assertEqual(b["eqtl_status"],"NOT_TESTED")
            self.assertEqual(c["gene_symbol"],"NEURL1")
            self.assertEqual(c["mol_qtl_task"],"KEEP_AS_EXTERNAL_ANCHOR")
if __name__=="__main__":unittest.main()
