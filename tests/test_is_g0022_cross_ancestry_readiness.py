"""Fixtures for fail-closed G0022 EAS/EUR GWAS coverage evaluation."""
import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_cross_ancestry_readiness import load_eas, collect, build


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with (gzip.open(path, "wt") if str(path).endswith(".gz") else path.open("w")) as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestG0022CrossAncestryReadiness(unittest.TestCase):
    def fixtures(self, root):
        g = root / "eas.tsv"
        records = []
        for pos, ref, alt, eaf in [(101, "G", "A", 0.20), (105, "C", "T", 0.25),
                                   (110, "A", "G", 0.12), (115, "G", "C", 0.4)]:
            records.append({"variant_id": f"12:{pos}:{ref}:{alt}", "chr": "12",
                "pos": str(pos), "ref": ref, "alt": alt,
                "status": "ALT_SIGNED_CONFIRMED", "reference_alt_af": eaf,
                "gwas_alt_eaf": eaf})
        write_tsv(g, records)
        s = root / "model.json"
        s.write_text(json.dumps({"credible_set_variant_ids": {"L1": ["12:101:G:A", "12:105:C:T"]},
                                 "n_snps": 4, "n_credible_sets": 1, "approximate_n_eff": 100}))
        eur = root / "eur.tsv.gz"
        write_tsv(eur, [
            {"chromosome": "1", "base_pair_location": "101", "effect_allele": "A",
             "other_allele": "G", "beta": "0.1", "standard_error": "0.1",
             "effect_allele_frequency": "0.3", "p_value": "0.25"},
            {"chromosome": "12", "base_pair_location": "110", "effect_allele": "A",
             "other_allele": "G", "beta": "0.4", "standard_error": "0.1",
             "effect_allele_frequency": "0.9", "p_value": "0.001"},
            {"chromosome": "12", "base_pair_location": "115", "effect_allele": "A",
             "other_allele": "T", "beta": "0.1", "standard_error": "0.1",
             "effect_allele_frequency": "0.3", "p_value": "0.2"},
        ])
        return g, s, eur

    def test_source_absence_separate_from_allele_mismatch(self):
        with tempfile.TemporaryDirectory() as t:
            g, s, eur = self.fixtures(Path(t))
            eas, cs, _ = load_eas(g, s)
            rows, summary = collect(eas, cs, eur)
            states = {x["variant_id"]: x["eur_source_coverage"] for x in rows}
            self.assertEqual(summary["eur_matched_eas_credible_set_variants"], 0)
            self.assertEqual(summary["eur_source_gwas_total_rows"], 3)
            self.assertEqual(summary["eur_source_g0022_region_rows"], 2)
            self.assertEqual(summary["eur_missing_positions"], 2)
            self.assertEqual(summary["eur_position_only_or_mismatched"], 1)
            self.assertEqual(states["12:110:A:G"], "EXACT_ALLELE_PAIR_MATCHED")
            self.assertEqual(states["12:115:G:C"], "POSITION_PRESENT_ALLELES_UNMATCHED")
            match = [x for x in rows if x["variant_id"] == "12:110:A:G"][0]
            self.assertEqual(float(match["eur_gwas_alt_eaf"]), 0.1)
            self.assertEqual(float(match["eur_alt_beta"]), -0.4)
            self.assertFalse(summary["ready_for_multi_ancestry_fine_mapping"])

    def test_outputs_and_sha_and_must_remain_blocked(self):
        with tempfile.TemporaryDirectory() as t:
            base = Path(t)
            g, s, eur = self.fixtures(base)
            class Args:
                eas_variants = g
                eas_summary = s
                eur_original = eur
                eur_signed_ld = None
                outdir = base / "output"
            j = build(Args)
            self.assertEqual(j["eur_exact_matched_eas_variants"], 1)
            self.assertEqual(j["eur_source_chromosome_rows"], 2)
            self.assertEqual(len(j["sources"]["eur_original_sha256"]), 64)
            self.assertEqual(len((Args.outdir / "G0022_EAS_EUR_SOURCE_VARIANT_COVERAGE.tsv").read_text().splitlines()), 5)
            self.assertFalse(j["multi_ancestry_direct_cs_ready"])

    def test_missing_cs_or_bad_eas_orientation_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            g, s, eur = self.fixtures(Path(t))
            source = json.loads(s.read_text())
            source["credible_set_variant_ids"]["L1"].append("12:500:C:T")
            s.write_text(json.dumps(source))
            with self.assertRaises(ValueError):
                load_eas(g, s)


if __name__ == "__main__":
    unittest.main()
