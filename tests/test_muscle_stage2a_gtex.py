import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "muscle_stage2a_gtex.py"

spec = importlib.util.spec_from_file_location("muscle_stage2a", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class MuscleStage2ATests(unittest.TestCase):
    def test_extract_records_v2_data_envelope(self):
        self.assertEqual(
            mod.extract_records({"data": [{"snpId": "rs1"}]}),
            [{"snpId": "rs1"}],
        )

    def test_select_grch38_mapping_prefers_chromosome(self):
        payload = {
            "mappings": [
                {"assembly_name": "GRCh38", "coord_system": "scaffold",
                 "seq_region_name": "PATCH", "start": 1},
                {"assembly_name": "GRCh38", "coord_system": "chromosome",
                 "seq_region_name": "3", "start": 100, "allele_string": "A/G"},
            ]
        }
        m = mod.select_grch38_mapping(payload)
        self.assertEqual(m["seq_region_name"], "3")
        self.assertEqual(m["start"], 100)

    def test_split_source_genes(self):
        self.assertEqual(mod.split_source_genes("SPTLC3;ISM1"), ["SPTLC3", "ISM1"])

    def test_candidate_gene_scoring(self):
        stage1 = [{
            "study_id": "S1", "rsid": "rs1", "gene": "GENEA", "p": "1e-6",
            "functional_evidence_points": "3",
            "stage1_priority": "MEDIUM_FUNCTIONAL_PRIORITY",
        }]
        eq = [{"study_id": "S1", "rsid": "rs1", "gene_symbol": "GENEA", "gencode_id": ""}]
        sq = [{"study_id": "S1", "rsid": "rs1", "gene_symbol": "GENEB", "gencode_id": ""}]
        rows = {r["gene"]: r for r in mod.build_candidate_genes(stage1, eq, sq)}
        self.assertEqual(rows["GENEA"]["stage2a_evidence_score"], 7)
        self.assertEqual(rows["GENEB"]["stage2a_evidence_score"], 5)

    def test_run_with_mocked_apis(self):
        header = [
            "study_id", "rsid", "chr", "pos", "gene", "p",
            "recalculated_tier", "functional_evidence_points", "stage1_priority"
        ]
        row = [
            "S1", "rs123", "3", "100", "GENEA", "1e-6",
            "B_suggestive", "3", "MEDIUM_FUNCTIONAL_PRIORITY"
        ]

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            inp = tmp / "stage1.tsv"
            inp.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")
            out = tmp / "out"

            ens = {
                "mappings": [{
                    "assembly_name": "GRCh38", "coord_system": "chromosome",
                    "seq_region_name": "3", "start": 100, "allele_string": "A/G"
                }]
            }
            gv = {"data": [{"snpId": "rs123", "variantId": "chr3_100_A_G_b38"}]}
            eq = {"data": [{
                "snpId": "rs123", "variantId": "chr3_100_A_G_b38",
                "geneSymbol": "GENEA", "gencodeId": "ENSG1.1",
                "tissueSiteDetailId": "Muscle_Skeletal",
                "datasetId": "gtex_v10", "pval": 1e-10, "slope": 0.2
            }]}
            sq = {"data": []}

            with patch.object(mod, "query_ensembl", return_value=(ens, "OK")), \
                 patch.object(mod, "query_gtex_variant", return_value=(gv, "OK")), \
                 patch.object(mod, "query_gtex_qtl", side_effect=[(eq, "OK"), (sq, "OK")]):
                summary = mod.run(inp, out, pause=0)

            self.assertEqual(summary["n_input_variants"], 1)
            self.assertEqual(summary["n_ensembl_grch38_resolved"], 1)
            self.assertEqual(summary["n_source_position_matches_grch38"], 1)
            self.assertEqual(summary["n_variants_with_significant_eqtl"], 1)
            self.assertEqual(summary["n_api_error_variants"], 0)
            self.assertTrue((out / "MUSCLE_STAGE2A_CANDIDATE_GENES.tsv").exists())
            self.assertTrue((out / "raw/ensembl/rs123.json").exists())


if __name__ == "__main__":
    unittest.main()
