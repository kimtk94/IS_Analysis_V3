"""Synthetic-only unit tests; no live human or GWAS data, no causal claims."""
import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_thesis_evidence_freeze import CKD_GENES, freeze, prepare_ckd, prepare_is


def tab(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_is(root, n=2):
    p = root / "is.tsv"
    rows = [
        {
            "gene_id_stable": f"ENSG{i:011d}",
            "canonical_gene_symbol": f"GENE{i}",
            "region_id": "G0001",
            "new_nominal_QTL_gene_source_status": "NOT_IN_EAS_30_GROUPS",
            "GWAS_and_QTL_alleles_coharmonized": "False",
            "matched_QTL_LD_verified": "False",
            "coloc_recomputed": "False",
            "ischemic_stroke_causal_gene": "NOT_ESTABLISHED",
        } for i in range(1, n + 1)
    ]
    tab(p, list(rows[0]), rows)
    return p


def write_ckd(root):
    p = root / "ckd.tsv"
    rows = [
        {
            "gene_symbol": g,
            "eur_abf_h4": "0.92" if i == 0 else "NA",
            "eur_susie_max_h4": "",
            "eas_egfr_abf_h4": "0.71" if i == 0 else "",
            "eas_bun_abf_h4": "",
            "top_cell_type": "proximal tubule" if i == 0 else "",
            "koges_rsid": "rs123" if i == 0 else "",
            "eur_stage2_class": ("susie_unresolved" if g == "UMOD"
                                 else "abf_susie_discordant" if g == "INHBC"
                                 else "core_shared_signal"),
            "eur_comparison_status": ("SuSiE_nonconverged_or_failed" if g == "UMOD"
                                      else "ABF_not_confirmed_by_SuSiE" if g == "INHBC"
                                      else "ABF_and_SuSiE_shared_signal"),
            # EUR source status alone must flag UMOD; no preassigned technical marker.
            "technical_status": "",
            "dual_egfr_status": "",
            "dual_bun_status": "",
        } for i, g in enumerate(CKD_GENES)
    ]
    tab(p, list(rows[0]), rows)
    return p


class EvidenceFreezeTests(unittest.TestCase):
    def test_is_missing_qtl_is_not_biological_negative(self):
        with tempfile.TemporaryDirectory() as d:
            rows, report = prepare_is(write_is(Path(d)), expected=2)
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(x["causal_claim_status"] == "NOT_ESTABLISHED" for x in rows))
            self.assertTrue(all(x["scientific_audit_status"].endswith("UNVERIFIED") for x in rows))
            self.assertTrue(all("QTL_DONOR_LD_NOT_VERIFIED" in x["limitations"] for x in rows))
            self.assertEqual(report["coloc_recomputed_rows"], 0)

    def test_invalid_is_denominator_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "count"):
                prepare_is(write_is(Path(d)), expected=2225)

    def test_duplicate_is_ensembl_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_is(Path(d))
            s = p.read_text()
            p.write_text(s.replace("ENSG00000000002", "ENSG00000000001"))
            with self.assertRaisesRegex(ValueError, "Invalid/duplicate"):
                prepare_is(p, expected=2)

    def test_region_cross_check_fails_on_missing_gene(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            p = write_is(root)
            q = root / "regions.tsv"
            tab(q, ["gene_id_stable"], [{"gene_id_stable": "ENSG00000000001"}])
            with self.assertRaisesRegex(ValueError, "omits or adds"):
                prepare_is(p, expected=2, region_path=q)

    def test_ckd_nine_gene_contract_and_no_koges_false_pass(self):
        with tempfile.TemporaryDirectory() as d:
            rows, stats = prepare_ckd(write_ckd(Path(d)))
            self.assertEqual(len(rows), 9)
            self.assertEqual(stats["eur_abf_present_rows"], 1)
            self.assertEqual(rows[0]["koges_anchor_present"], "TRUE")
            self.assertEqual(rows[0]["koges_longitudinal_validated"], "FALSE")
            self.assertIn("MULTISIGNAL_TECHNICAL_UNRESOLVED", rows[-1]["limitations"])
            self.assertIn("EUR_MULTISIGNAL_UNRESOLVED", rows[-1]["limitations"])
            inhbc = next(x for x in rows if x["gene_symbol"] == "INHBC")
            self.assertIn("EUR_ABF_SUSIE_DISCORDANT", inhbc["limitations"])

    def test_ckd_duplicate_or_missing_candidate_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_ckd(Path(d))
            s = p.read_text().replace("UMOD\t", "ACP1\t")
            p.write_text(s)
            with self.assertRaisesRegex(ValueError, "exactly the 9"):
                prepare_ckd(p)

    def test_ckd_out_of_range_h4_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_ckd(Path(d))
            p.write_text(p.read_text().replace("0.92", "1.25"))
            with self.assertRaisesRegex(ValueError, "outside"):
                prepare_ckd(p)

    def test_manifest_hash_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            gene, ckd = write_is(root), write_ckd(root)
            out = root / "frozen"
            manifest = freeze(gene, ckd, out, expected_is=2)
            self.assertFalse(manifest["scientific_validation_completed"])
            self.assertEqual(manifest["sources"]["is.tsv"]["sha256"],
                             hashlib.sha256(gene.read_bytes()).hexdigest())
            with (out / "CKD_EVIDENCE_READINESS_FROZEN.tsv").open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(rows), 9)
            self.assertEqual(json.loads((out / "EVIDENCE_FREEZE_PROVENANCE.json").read_text())
                             ["coverage"]["IS"]["gene_rows"], 2)
            with self.assertRaisesRegex(ValueError, "already populated"):
                freeze(gene, ckd, out, expected_is=2)

    def test_missing_input_fails_before_output(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            out = root / "output"
            with self.assertRaisesRegex(ValueError, "Missing input"):
                freeze(root / "not_here.tsv", None, out, expected_is=2)
            self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
