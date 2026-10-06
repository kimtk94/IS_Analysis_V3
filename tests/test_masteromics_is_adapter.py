import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from masteromics.is_adapter import build_baseline


def write_tsv(path: Path, header, rows, gz=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gz else open
    kwargs = {"mode": "wt", "encoding": "utf-8", "newline": ""} if gz else {"mode": "w", "encoding": "utf-8", "newline": ""}
    with opener(path, **kwargs) as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class IsAdapterTest(unittest.TestCase):
    def populate(self, root: Path):
        bbj = root / "data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"
        write_tsv(
            bbj,
            ["dataset","phenotype","build","chr","pos","effect_allele","other_allele","beta","se","p","eaf","n","variant_id"],
            [{"dataset":"BBJ","phenotype":"ischemic_stroke","build":"GRCh37","chr":"4","pos":"81182554","effect_allele":"C","other_allele":"T","beta":"0.1","se":"0.02","p":"1e-9","eaf":"0.2","n":"174686","variant_id":"4:81182554:T:C"}],
            gz=True,
        )
        giga = root / "data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz"
        write_tsv(
            giga,
            ["dataset","phenotype","ancestry","build","chr","pos","effect_allele","other_allele","beta","se","p","eaf","variant_id"],
            [{"dataset":"GCST90104545","phenotype":"AIS","ancestry":"EAS","build":"GRCh37","chr":"4","pos":"81182554","effect_allele":"C","other_allele":"T","beta":"0.08","se":"0.03","p":"0.01","eaf":"0.2","variant_id":"4:81182554:T:C"}],
            gz=True,
        )
        regions = root / "results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv"
        write_tsv(
            regions,
            ["locus_id","chr","region_start","region_end","lead_variant","role"],
            [
                {"locus_id":"BBJ_IS_L001","chr":"4","region_start":"1","region_end":"2","lead_variant":"4:1:A:G","role":"PRIMARY"},
                {"locus_id":"BBJ_IS_L002","chr":"10","region_start":"1","region_end":"2","lead_variant":"10:1:A:G","role":"PRIMARY"},
                {"locus_id":"BBJ_IS_L003","chr":"12","region_start":"1","region_end":"2","lead_variant":"12:1:A:G","role":"PRIMARY_WIDE_REQUIRES_LD_SPLIT"},
                {"locus_id":"BBJ_IS_L004","chr":"13","region_start":"1","region_end":"2","lead_variant":"13:1:A:G","role":"PRIMARY"},
                {"locus_id":"BBJ_IS_L005","chr":"19","region_start":"1","region_end":"2","lead_variant":"19:1:A:G","role":"EXCLUDE_PRIMARY_RARE_LOW_INFO"},
            ],
        )
        comparison = root / "results/is/stage4_cross_eas/fine_mapping/BBJ_GIGASTROKE_CS_COMPARISON.tsv"
        write_tsv(
            comparison,
            ["phenotype","locus","bbj_top_variant","bbj_top_pip","giga_top_variant","giga_top_pip","cs_jaccard"],
            [{"phenotype":"AIS","locus":"BBJ_IS_L001","bbj_top_variant":"4:1:A:G","bbj_top_pip":"0.5","giga_top_variant":"4:1:A:G","giga_top_pip":"0.4","cs_jaccard":"1"}],
        )
        convergence = root / "results/is/stage5_functional/phase9d_literature_benchmark/IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv"
        header = ["gene","locus","role","mechanism_branch","best_abf_h4","best_susie_h4","best_abf_dataset","best_abf_tissue","interpretation"]
        write_tsv(
            convergence,
            header,
            [
                {"gene":"FGF5","locus":"BBJ_IS_L001","role":"CORE","mechanism_branch":"REGULATORY_EXPRESSION_PROTEIN","best_abf_h4":"0.78","best_susie_h4":"0.75","best_abf_dataset":"GTEx","best_abf_tissue":"Brain","interpretation":"SUPPORTED"},
                {"gene":"ALDH2","locus":"BBJ_IS_L003","role":"CORE","mechanism_branch":"CODING_PROTEIN_METABOLIC","best_abf_h4":"0.08","best_susie_h4":"","best_abf_dataset":"GTEx","best_abf_tissue":"Artery","interpretation":"RS671"},
                {"gene":"SH3PXD2A","locus":"BBJ_IS_L002","role":"CORE","mechanism_branch":"CELL_SPECIFIC_REGULATORY","best_abf_h4":"0.35","best_susie_h4":"","best_abf_dataset":"GTEx","best_abf_tissue":"Artery","interpretation":"CELL_PRIORITY"},
                {"gene":"COL4A2","locus":"BBJ_IS_L004","role":"CORE","mechanism_branch":"VASCULAR_STRUCTURAL_REGULATORY","best_abf_h4":"0.33","best_susie_h4":"","best_abf_dataset":"GTEx","best_abf_tissue":"Brain","interpretation":"VASCULAR_PRIORITY"},
            ],
        )
        mouse = root / "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv"
        write_tsv(
            mouse,
            ["component","status"],
            [
                {"component":"MOUSE_LAYER_CONCLUSION","status":"CELLTYPE_LOCALIZATION_SUPPORTED_TARGET_SPECIFIC_STROKE_DGE_NOT_ESTABLISHED"},
                {"component":"PRIORITIZED_TARGETS_FDR05_PRIMARY","status":"0"},
            ],
        )
        handoff = root / "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_DRIVE_HANDOFF.tsv"
        write_tsv(
            handoff,
            ["component","status"],
            [
                {"component":"REMOTE_PATH","status":"gdrive:fixture.rds.gz"},
                {"component":"COLAB_STATUS","status":"READY"},
            ],
        )
        return bbj

    def test_frozen_baseline_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.populate(root)
            out = root / "out"
            status = build_baseline(root, out)
            self.assertEqual(status["adapter_status"], "PASS")
            self.assertEqual(status["scientific_status"], "FROZEN_RESULTS_NOT_RECOMPUTED")
            self.assertEqual(status["checksum_status"], "COMPLETE")
            self.assertEqual(status["core_genes"], ["FGF5","ALDH2","SH3PXD2A","COL4A2"])
            self.assertEqual(status["primary_loci"], ["BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L003","BBJ_IS_L004"])
            self.assertEqual(status["human_vascular_colab_status"], "READY")
            with (out / "IS_BASELINE_CORE_EVIDENCE.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([x["gene_symbol"] for x in rows], ["FGF5","ALDH2","SH3PXD2A","COL4A2"])
            with (out / "IS_BASELINE_SOURCE_AUDIT.tsv").open() as handle:
                source_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(source_rows), 7)
            self.assertTrue(all(x["checksum_status"] == "PINNED" for x in source_rows))

    def test_missing_required_source_column_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bbj = self.populate(root)
            write_tsv(
                bbj,
                ["dataset","phenotype","build","chr","pos","effect_allele","other_allele","beta","se","p","eaf","variant_id"],
                [{"dataset":"BBJ","phenotype":"ischemic_stroke","build":"GRCh37","chr":"4","pos":"1","effect_allele":"A","other_allele":"G","beta":"0.1","se":"0.02","p":"1e-9","eaf":"0.2","variant_id":"4:1:A:G"}],
                gz=True,
            )
            with self.assertRaisesRegex(ValueError, "missing required columns: n"):
                build_baseline(root, root / "out")


if __name__ == "__main__":
    unittest.main()
