"""Register frozen ischemic-stroke evidence against the MasterOmics IS stage catalog.

This migration layer never recomputes scientific results. It requires a successful
pinned IS baseline audit, validates frozen scientific guardrails, hashes stage source
artifacts, and records which IS stages are migrated versus still scientifically blocked.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

from .architecture import IS_CATALOG

READY = "MIGRATED_FROZEN"
BLOCKED_AUTHOR = "BLOCKED_AUTHOR_ANNOTATION"
BLOCKED_DEP = "BLOCKED_DEPENDENCY"

LEGACY_STAGE_SOURCES = {
    "acquisition": [
        "data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz",
        "data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz",
    ],
    "normalize": [
        "data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz",
        "data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz",
    ],
    "gwas_loci": [
        "results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv",
    ],
    "finemap": [
        "results/is/stage3_finemap/japan/bbj/SUSIE_NEFF_MASTER_SUMMARY.tsv",
    ],
    "cross_ancestry": [
        "results/is/stage4_cross_eas/fine_mapping/BBJ_GIGASTROKE_CS_COMPARISON.tsv",
    ],
    "molecular_coloc": [
        "results/is/stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
        "results/is/stage5_functional/phase9c_convergence/COLOC_SUSIE_SIGNAL_PAIRS.tsv",
    ],
    "mechanism": [
        "results/is/stage5_functional/phase9d_literature_benchmark/IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv",
    ],
    "celltype": [
        "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv",
        "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE225948_TARGET_TESTABILITY.tsv",
    ],
    "human_annotation": [
        "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_DRIVE_HANDOFF.tsv",
        "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_RDS_SHA256.txt",
    ],
}

EXPECTED_LOCI = {"BBJ_IS_L001", "BBJ_IS_L002", "BBJ_IS_L003", "BBJ_IS_L004"}
EXPECTED_CORE = {"FGF5", "ALDH2", "SH3PXD2A", "COL4A2"}


def _sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="	"))


def _artifact(path: Path):
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError("Missing/empty frozen stage source: " + str(path))
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def validate_baseline(baseline_dir: Path):
    status_path = baseline_dir / "IS_BASELINE_STATUS.json"
    coverage_path = baseline_dir / "IS_STAGE_COVERAGE.tsv"
    source_path = baseline_dir / "IS_BASELINE_SOURCE_AUDIT.tsv"
    for path in [status_path, coverage_path, source_path]:
        if not path.is_file():
            raise ValueError("Missing IS baseline adapter artifact: " + str(path))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("adapter_status") != "PASS" or status.get("checksum_status") != "COMPLETE":
        raise ValueError("IS baseline adapter has not passed pinned validation")
    if status.get("source_count") != 10:
        raise ValueError("IS baseline source count changed")
    cov = {r["stage"]: r["status"] for r in _rows(coverage_path)}
    expected_ready = {
        "acquisition","source_qc","normalize","gwas_loci","finemap",
        "cross_ancestry","molecular_coloc","mechanism","celltype",
    }
    if {k for k,v in cov.items() if v == "FROZEN_ARTIFACT_READY"} != expected_ready:
        raise ValueError("IS frozen stage coverage changed")
    if cov.get("human_annotation") != "PENDING_AUTHOR_ANNOTATION":
        raise ValueError("IS human annotation gate changed")
    return status


def validate_frozen_science(root: Path):
    finemap = _rows(root / LEGACY_STAGE_SOURCES["finemap"][0])
    observed = {r["locus"]: r for r in finemap}
    if set(observed) != EXPECTED_LOCI:
        raise ValueError("Frozen BBJ fine-map locus set changed")
    bad = [
        locus for locus,row in observed.items()
        if row.get("status") != "PASS" or row.get("converged", "").upper() != "TRUE"
    ]
    if bad:
        raise ValueError("Frozen BBJ SuSiE convergence changed: " + ",".join(sorted(bad)))

    convergence = _rows(root / LEGACY_STAGE_SOURCES["mechanism"][0])
    core = {r["gene"] for r in convergence if r.get("role") == "CORE"}
    missing_core = EXPECTED_CORE - core
    if missing_core:
        raise ValueError("Required frozen core mechanism genes missing: " + ",".join(sorted(missing_core)))

    handoff = _rows(root / LEGACY_STAGE_SOURCES["human_annotation"][0])
    handoff_map = {r["component"]: r["status"] for r in handoff}
    if handoff_map.get("REMOTE_RDS") != "PRESENT_SIZE_MATCH":
        raise ValueError("Human vascular RDS handoff is not size-matched")
    if handoff_map.get("COLAB_STATUS") != "READY":
        raise ValueError("Human vascular author annotation is not Colab-ready")


def build(root: Path, baseline_dir: Path, outdir: Path):
    root = root.resolve()
    baseline_dir = baseline_dir.resolve()
    outdir = outdir.resolve()
    baseline_status = validate_baseline(baseline_dir)
    validate_frozen_science(root)

    stage_ids = [x[0] for x in IS_CATALOG]
    if stage_ids != [
        "acquisition","source_qc","normalize","gwas_loci","finemap",
        "cross_ancestry","molecular_coloc","mechanism","celltype",
        "human_annotation","evidence","report",
    ]:
        raise ValueError("MasterOmics IS stage catalog changed")

    rows = []
    payloads = {}
    for stage_id in stage_ids:
        if stage_id in {"evidence","report"}:
            status = BLOCKED_DEP
        elif stage_id == "human_annotation":
            status = BLOCKED_AUTHOR
        else:
            status = READY

        paths = []
        if stage_id == "source_qc":
            paths = [
                baseline_dir / "IS_BASELINE_SOURCE_AUDIT.tsv",
                baseline_dir / "IS_STAGE_COVERAGE.tsv",
                baseline_dir / "IS_BASELINE_STATUS.json",
            ]
        else:
            paths = [root / rel for rel in LEGACY_STAGE_SOURCES.get(stage_id, [])]
        artifacts = [_artifact(path) for path in paths]

        payload = {
            "schema_version": 1,
            "project": "ischemic_stroke",
            "stage_id": stage_id,
            "migration_status": status,
            "scientific_status": "FROZEN_RESULTS_NOT_RECOMPUTED" if status == READY else "NOT_EXECUTED",
            "source_artifacts": artifacts,
        }
        if stage_id == "human_annotation":
            payload["blocker"] = "AUTHOR_VALIDATED_GSE256493_HUMAN_VASCULAR_ANNOTATION_REQUIRED"
            payload["next_action"] = "RUN_COLAB_HUMAN_AUTHOR_ANNOTATION"
        elif stage_id == "evidence":
            payload["blocker"] = "HUMAN_ANNOTATION_REQUIRED"
        elif stage_id == "report":
            payload["blocker"] = "FINAL_EVIDENCE_REQUIRED"
        payloads[stage_id] = payload
        rows.append({
            "stage_id": stage_id,
            "migration_status": status,
            "scientific_status": payload["scientific_status"],
            "n_source_artifacts": len(artifacts),
            "blocker": payload.get("blocker", ""),
        })

    outdir.mkdir(parents=True, exist_ok=True)
    stage_dir = outdir / "stages"
    stage_dir.mkdir(parents=True, exist_ok=True)
    for stage_id,payload in payloads.items():
        (stage_dir / f"{stage_id}.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )

    with (outdir / "IS_STAGE_MIGRATION.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    ready = sum(r["migration_status"] == READY for r in rows)
    manifest = {
        "schema_version": 1,
        "project": "ischemic_stroke",
        "baseline_schema_version": baseline_status.get("schema_version"),
        "catalog_stage_count": len(rows),
        "migrated_frozen_stages": ready,
        "blocked_stages": [r["stage_id"] for r in rows if r["migration_status"] != READY],
        "current_gate": "human_annotation",
        "current_blocker": "AUTHOR_VALIDATED_GSE256493_HUMAN_VASCULAR_ANNOTATION_REQUIRED",
        "masteromics_binding_status": "MIGRATION_PROVENANCE_REGISTERED_NOT_EXECUTABLE",
        "scientific_guardrail": (
            "Migration registers frozen legacy evidence; it does not recompute or "
            "upgrade scientific claims."
        ),
    }
    (outdir / "IS_STAGE_MIGRATION.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy_root", type=Path)
    parser.add_argument("baseline_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args(argv)
    try:
        print(json.dumps(build(args.legacy_root, args.baseline_dir, args.outdir), indent=2))
        return 0
    except (OSError, ValueError, KeyError, csv.Error, json.JSONDecodeError) as exc:
        print("IS stage migration error: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
