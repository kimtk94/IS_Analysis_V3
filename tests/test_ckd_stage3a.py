#!/usr/bin/env python3
import importlib.util
import unittest
from pathlib import Path

P = Path(__file__).resolve().parents[1] / "scripts" / "run_ckd_stage3a_kidney_evidence.py"
spec = importlib.util.spec_from_file_location("m", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class Stage3ATests(unittest.TestCase):
    def test_valid_xlsx_guard(self):
        import tempfile, zipfile
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.xlsx"
            p.write_bytes(b"<html>not xlsx</html>")
            self.assertFalse(m.is_valid_xlsx(p))
            with zipfile.ZipFile(p, "w") as z:
                z.writestr("[Content_Types].xml", "<Types/>")
                z.writestr("xl/worksheets/sheet1.xml", "<worksheet/>")
                z.writestr("padding.bin", b"0" * 120000)
            self.assertTrue(m.is_valid_xlsx(p))

    def test_gene_match_exact_and_multivalue(self):
        genes = {"F12", "GSTA3", "ACP1"}
        self.assertEqual(m.exact_gene_in_cell("F12", genes), "F12")
        self.assertEqual(
            m.exact_gene_in_cell("ABC;GSTA3;XYZ", genes), "GSTA3"
        )
        self.assertIsNone(m.exact_gene_in_cell("F12A1", genes))

    def test_table_number_and_evidence(self):
        self.assertEqual(
            m.sheet_table_number("Supplementary Table 3"), 3
        )
        self.assertEqual(
            m.evidence_type_from_table(3),
            "kidney_cis_pqtl_significant",
        )
        self.assertEqual(
            m.evidence_type_from_table(13),
            "kidney_cis_eqtl_significant_same_study",
        )
        self.assertEqual(
            m.evidence_type_from_table(16),
            "kidney_qtl_egfr_colocalization",
        )

    def test_stage2_classification(self):
        self.assertEqual(
            m.classify_stage2({
                "comparison_status": "ABF_and_SuSiE_shared_signal"
            }),
            "core_shared_signal",
        )
        self.assertEqual(
            m.classify_stage2({
                "comparison_status": "SuSiE_nonconverged_or_failed"
            }),
            "susie_unresolved",
        )
        self.assertEqual(
            m.classify_stage2({
                "comparison_status": "ABF_not_confirmed_by_SuSiE",
                "max_PP.H4": "0.77",
            }),
            "abf_supported_susie_weaker",
        )
        self.assertEqual(
            m.classify_stage2({
                "comparison_status": "ABF_not_confirmed_by_SuSiE",
                "max_PP.H4": "3e-14",
            }),
            "abf_susie_discordant",
        )


if __name__ == "__main__":
    unittest.main()
