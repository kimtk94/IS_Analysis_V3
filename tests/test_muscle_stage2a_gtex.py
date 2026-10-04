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

    def test_parse_b37_variant(self):
        self.assertEqual(
            mod.parse_gtex_b37_variant_id("3_77476889_C_T_b37"),
            ("3", 77476889),
        )
        self.assertEqual(
            mod.source_matches_any_b37("3", 77476889, ["3_77476889_C_T_b37"]),
            "1",
        )

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

    def test_run_queries_qtl_by_variant_id(self):
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
                    "seq_region_name": "3", "start": 200, "allele_string": "A/G"
                }]
            }
            gv = {"data": [{
                "snpId": "rs123",
                "variantId": "chr3_200_A_G_b38",
                "b37VariantId": "3_100_A_G_b37",
            }]}
            eq = {"data": [{
                "snpId": "rs123", "variantId": "chr3_200_A_G_b38",
                "geneSymbol": "GENEA", "gencodeId": "ENSG1.1",
                "tissueSiteDetailId": "Muscle_Skeletal",
                "datasetId": "gtex_v10", "pValue": 1e-10, "nes": 0.2
            }]}
            sq = {"data": []}

            calls = []
            def fake_qtl(variant_id, dataset, tissue, qtl):
                calls.append((variant_id, dataset, tissue, qtl))
                return (eq, "OK") if qtl == "eqtl" else (sq, "OK")

            with patch.object(mod, "query_ensembl", return_value=(ens, "OK")), \
                 patch.object(mod, "query_gtex_variant", return_value=(gv, "OK")), \
                 patch.object(mod, "query_gtex_qtl", side_effect=fake_qtl):
                summary = mod.run(inp, out, pause=0)

            self.assertEqual(summary["version"], "2A.2")
            self.assertEqual(summary["n_input_variants"], 1)
            self.assertEqual(summary["n_source_position_matches_grch38"], 0)
            self.assertEqual(summary["n_source_position_matches_gtex_b37"], 1)
            self.assertEqual(summary["n_variants_with_significant_eqtl"], 1)
            self.assertEqual(summary["n_api_error_variants"], 0)
            self.assertEqual(calls[0][0], "chr3_200_A_G_b38")
            self.assertEqual(calls[1][0], "chr3_200_A_G_b38")

    def test_unresolved_gtex_variant_is_not_api_error(self):
        header = [
            "study_id", "rsid", "chr", "pos", "gene", "p",
            "recalculated_tier", "functional_evidence_points", "stage1_priority"
        ]
        row = [
            "S1", "rs404", "1", "100", "GENEX", "1e-6",
            "B_suggestive", "1", "EXPLORATORY"
        ]

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            inp = tmp / "stage1.tsv"
            inp.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")
            out = tmp / "out"

            ens = {"mappings": [{
                "assembly_name": "GRCh38", "coord_system": "chromosome",
                "seq_region_name": "1", "start": 120, "allele_string": "A/G"
            }]}

            with patch.object(mod, "query_ensembl", return_value=(ens, "OK")), \
                 patch.object(mod, "query_gtex_variant", return_value=({"data": []}, "OK")), \
                 patch.object(mod, "query_gtex_qtl") as qtl_mock:
                summary = mod.run(inp, out, pause=0)

            qtl_mock.assert_not_called()
            self.assertEqual(summary["n_gtex_variant_unresolved"], 1)
            self.assertEqual(summary["gtex_variant_unresolved_rsids"], ["rs404"])
            self.assertEqual(summary["n_api_error_variants"], 0)


if __name__ == "__main__":
    unittest.main()
