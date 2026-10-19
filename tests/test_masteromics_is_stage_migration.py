import csv
import json
import tempfile
import unittest
from pathlib import Path

from masteromics.is_stage_migration import build


def write_tsv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class IsStageMigrationTest(unittest.TestCase):
    def populate(self, root: Path, baseline: Path, colab_status="READY", converged="TRUE"):
        for rel in [
            "data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz",
            "data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz",
            "results/is/stage4_cross_eas/fine_mapping/BBJ_GIGASTROKE_CS_COMPARISON.tsv",
            "results/is/stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
            "results/is/stage5_functional/phase9c_convergence/COLOC_SUSIE_SIGNAL_PAIRS.tsv",
            "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv",
            "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE225948_TARGET_TESTABILITY.tsv",
            "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_RDS_SHA256.txt",
        ]:
            p = root / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"fixture\n")

        write_tsv(
            root / "results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv",
            ["locus_id","chr","region_start","region_end","lead_variant","role"],
            [
                {"locus_id":"BBJ_IS_L001","chr":"4","region_start":"1","region_end":"2","lead_variant":"4:1:A:G","role":"PRIMARY"},
                {"locus_id":"BBJ_IS_L002","chr":"10","region_start":"1","region_end":"2","lead_variant":"10:1:A:G","role":"PRIMARY"},
                {"locus_id":"BBJ_IS_L003","chr":"12","region_start":"1","region_end":"2","lead_variant":"12:1:A:G","role":"PRIMARY"},
                {"locus_id":"BBJ_IS_L004","chr":"13","region_start":"1","region_end":"2","lead_variant":"13:1:A:G","role":"PRIMARY"},
            ],
        )
        write_tsv(
            root / "results/is/stage3_finemap/japan/bbj/SUSIE_NEFF_MASTER_SUMMARY.tsv",
            ["locus","status","converged","n_eff","nvar","n_cs","top_variant","top_pip","error"],
            [
                {"locus":x,"status":"PASS","converged":converged,"n_eff":"78894","nvar":"100","n_cs":"1","top_variant":"1:1:A:G","top_pip":"0.5","error":""}
                for x in ["BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L003","BBJ_IS_L004"]
            ],
        )
        write_tsv(
            root / "results/is/stage5_functional/phase9d_literature_benchmark/IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv",
            ["gene","locus","role","mechanism_branch"],
            [
                {"gene":"FGF5","locus":"BBJ_IS_L001","role":"CORE","mechanism_branch":"REGULATORY_EXPRESSION_PROTEIN"},
                {"gene":"ALDH2","locus":"BBJ_IS_L003","role":"CORE","mechanism_branch":"CODING_PROTEIN_METABOLIC"},
                {"gene":"SH3PXD2A","locus":"BBJ_IS_L002","role":"CORE","mechanism_branch":"CELL_SPECIFIC_REGULATORY"},
                {"gene":"COL4A2","locus":"BBJ_IS_L004","role":"CORE","mechanism_branch":"VASCULAR_STRUCTURAL_REGULATORY"},
            ],
        )
        write_tsv(
            root / "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_DRIVE_HANDOFF.tsv",
            ["component","status"],
            [
                {"component":"REMOTE_RDS","status":"PRESENT_SIZE_MATCH"},
                {"component":"COLAB_STATUS","status":colab_status},
            ],
        )

        baseline.mkdir(parents=True, exist_ok=True)
        (baseline / "IS_BASELINE_STATUS.json").write_text(json.dumps({
            "schema_version":3,
            "adapter_status":"PASS",
            "checksum_status":"COMPLETE",
            "source_count":10,
        }) + "\n")
        write_tsv(
            baseline / "IS_STAGE_COVERAGE.tsv",
            ["stage","status","source_ids","note"],
            [
                *[
                    {"stage":x,"status":"FROZEN_ARTIFACT_READY","source_ids":"fixture","note":""}
                    for x in [
                        "acquisition","source_qc","normalize","gwas_loci","finemap",
                        "cross_ancestry","molecular_coloc","mechanism","celltype",
                    ]
                ],
                {"stage":"human_annotation","status":"PENDING_AUTHOR_ANNOTATION","source_ids":"fixture","note":""},
                {"stage":"evidence","status":"BLOCKED_UPSTREAM","source_ids":"","note":""},
                {"stage":"report","status":"BLOCKED_UPSTREAM","source_ids":"","note":""},
            ],
        )
        (baseline / "IS_BASELINE_SOURCE_AUDIT.tsv").write_text("dataset_id\tstatus\nfixture\tPASS\n")

    def test_registers_nine_frozen_stages_and_blocks_three(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "legacy"
            baseline = Path(tmp) / "baseline"
            out = Path(tmp) / "out"
            self.populate(root, baseline)
            result = build(root, baseline, out)
            self.assertEqual(result["catalog_stage_count"], 12)
            self.assertEqual(result["migrated_frozen_stages"], 9)
            self.assertEqual(result["blocked_stages"], ["human_annotation","evidence","report"])
            self.assertEqual(result["current_gate"], "human_annotation")
            with (out / "IS_STAGE_MIGRATION.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_stage = {r["stage_id"]: r for r in rows}
            self.assertEqual(by_stage["finemap"]["migration_status"], "MIGRATED_FROZEN")
            self.assertEqual(by_stage["human_annotation"]["migration_status"], "BLOCKED_AUTHOR_ANNOTATION")
            human = json.loads((out / "stages/human_annotation.json").read_text())
            self.assertEqual(human["scientific_status"], "NOT_EXECUTED")
            self.assertIn("AUTHOR_VALIDATED", human["blocker"])

    def test_nonconverged_finemap_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "legacy"
            baseline = Path(tmp) / "baseline"
            self.populate(root, baseline, converged="FALSE")
            with self.assertRaisesRegex(ValueError, "SuSiE convergence changed"):
                build(root, baseline, Path(tmp) / "out")

    def test_not_ready_human_handoff_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "legacy"
            baseline = Path(tmp) / "baseline"
            self.populate(root, baseline, colab_status="NOT_READY")
            with self.assertRaisesRegex(ValueError, "not Colab-ready"):
                build(root, baseline, Path(tmp) / "out")


if __name__ == "__main__":
    unittest.main()
