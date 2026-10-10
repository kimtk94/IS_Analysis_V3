"""Deterministic fixture tests for genotype-based JPT/EAS rs671 LD sensitivity."""
import csv
import math
import struct
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_jpt_eas_ld_stability import (
    LEAD, load_panel, get_reference_markers, load_dosages, summarize,
)


def write_tsv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestJptEasLD(unittest.TestCase):
    def setup_files(self, root):
        panel = root / "panel.txt"
        samples = ["A", "B", "C", "D", "E", "F", "G", "H"]
        write_tsv(panel, [{"sample": x, "pop": "JPT" if i < 4 else "CHB",
                           "super_pop": "EAS", "gender": "male"} for i, x in enumerate(samples)])
        variants = root / "variants.tsv"
        second = "12:112168009:G:A"
        write_tsv(variants, [{"variant_id": LEAD, "status": "ALT_SIGNED_CONFIRMED"},
                             {"variant_id": second, "status": "ALT_SIGNED_CONFIRMED"}])
        # rows/columns [lead, second]. Real pipeline checks exact dense row-major f64.
        ld = root / "ld.f64"
        ld.write_bytes(struct.pack("<4d", 1, 0.75, 0.75, 1))
        traw = root / "geno.traw"
        header = ["CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"] + [s + "_" + s for s in samples]
        vals = {
            LEAD: [0, 1, 2, 0, 1, 2, 1, 0],
            second: [0, 1, 2, 1, 1, 2, 1, 0],
        }
        with traw.open("w") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
            for vid in [LEAD, second]:
                parts = vid.split(":")
                # value is REF dosage, ALT = 2 - REF
                w.writerow([12, vid, 0, parts[1], parts[2], parts[3]] + vals[vid])
        return panel, variants, ld, traw

    def test_genotype_selection_and_population(self):
        with tempfile.TemporaryDirectory() as t:
            panel, var, ld, traw = self.setup_files(Path(t))
            sample_pop = load_panel(panel)
            selected = get_reference_markers(var, ld, cutoff=.5)
            self.assertEqual(len(selected), 2)
            samples, geno = load_dosages(traw, selected, sample_pop)
            self.assertEqual(len(samples), 8)
            self.assertAlmostEqual(geno[LEAD][0], 2)
            self.assertAlmostEqual(geno[LEAD][2], 0)
            self.assertEqual(int(np.sum(np.isfinite(geno[LEAD]))), 8)

    def test_bad_allele_must_raise(self):
        with tempfile.TemporaryDirectory() as t:
            panel, var, ld, traw = self.setup_files(Path(t))
            original = traw.read_text()
            traw.write_text(original.replace("112168009\tG\tA", "112168009\tT\tA"))
            with self.assertRaises(ValueError):
                load_dosages(traw, get_reference_markers(var, ld), load_panel(panel))

    def test_mismatched_signed_ld_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            panel, var, ld, traw = self.setup_files(Path(t))
            pop = load_panel(panel)
            chosen = get_reference_markers(var, ld)
            samples, geno = load_dosages(traw, chosen, pop)
            # The fixture's LD is deliberately not the actual genotype r for 8 samples.
            with self.assertRaises(ValueError):
                summarize(samples, geno, chosen, pop)


if __name__ == "__main__":
    unittest.main()
