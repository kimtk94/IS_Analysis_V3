"""No-network tests for molecular-QTL source gates and cached ImmuNexUT parsing."""
import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_is_eas_molecular_qtl_source_gates import SOURCES, assess, gate
from audit_is_immunexut_aldh2_browser import FOUR_CS, parse_page, build


class TestEQTLSourceGates(unittest.TestCase):
    def test_jctf_v2_v3_fail_closed_and_immunexut_eligible_for_screen(self):
        states = {s["dataset"]: s["gate"] for s in assess(SOURCES)}
        self.assertEqual(states["NHA000172"], "BLOCKED_THRESHOLD_SELECTED")
        self.assertEqual(states["NHA000193"], "BLOCKED_THRESHOLD_SELECTED")
        self.assertEqual(states["E-GEAD-420"], "SOURCE_PROMISING_GENOTYPE_COVERAGE_UNKNOWN")
        self.assertEqual(states["ZENODO_21296030"], "BLOCKED_ACCESS")
        self.assertTrue(all(not s["can_publish_causal_mediator"] for s in assess(SOURCES)))

    def test_gate_failure_for_source_claiming_valid_causal_conclusion(self):
        rows = [dict(SOURCES[2], publication_ready=True)]
        with self.assertRaises(ValueError):
            assess(rows)

    def example_html(self, skipped=None, duplicate=False):
        items = []
        for i, (_, rsid, pos) in enumerate(FOUR_CS):
            if rsid == skipped:
                continue
            item = {"gene_symbol":"ALDH2","snp_id":rsid,
                    "snp_url":"eqtlSnpsLink?snp_id="+rsid,
                    "cell_type":"Neu","chromosome":"chr12","position":pos,
                    "eqtl_pval1":4.88802e-9, "eqtl_effect_beta1":-0.207695}
            items.append(json.dumps(item,separators=(",",":")))
            if duplicate and rsid == "rs671":
                items.append(json.dumps(item,separators=(",",":")))
        return "<html><body>var eqtl = ["+",".join(items)+"];</body></html>"

    def test_four_cs_neutrophil_extract(self):
        rows, j = parse_page(self.example_html())
        self.assertEqual(len(rows),4)
        self.assertTrue(j["identical_four_p_and_beta"])
        self.assertEqual(j["new_valid_coloc"],0)
        self.assertFalse(j["gwas_qtl_allele_direction_compared"])
        self.assertEqual({x["qtl_rsid"] for x in rows},{x[1] for x in FOUR_CS})

    def test_missing_or_duplicate_must_fail(self):
        with self.assertRaisesRegex(ValueError,"Missing/duplicate"):
            parse_page(self.example_html(skipped="rs671"))
        with self.assertRaisesRegex(ValueError,"Missing/duplicate"):
            parse_page(self.example_html(duplicate=True))

    def test_no_network_read_only_fixture_and_sha(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t)
            infile=base/"saved.html"
            infile.write_text(self.example_html())
            result=build(argparse.Namespace(source_html=infile,outdir=base/"output",fetch=False,overwrite=False))
            self.assertEqual(len(result["html_sha256"]),64)
            self.assertEqual(len((base/"output/G0022_IMMUNEXUT_ALDH2_FOUR_CS_NEU_EQTL.tsv").read_text().splitlines()),5)
            self.assertFalse(result["validated_causal_gene"])


if __name__ == "__main__":
    unittest.main()
