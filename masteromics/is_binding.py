"""Bind frozen ischemic-stroke migration records into a MasterOmics blueprint.

The binder connects only stages already registered as MIGRATED_FROZEN. It leaves
human_annotation, evidence and report unbound so a frozen migration can never be
mistaken for a complete newly executed stroke analysis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .architecture import IS_CATALOG, inspect, validate
from .engine import atomic_json, sha

READY = "MIGRATED_FROZEN"
EXPECTED_BLOCKED = ["human_annotation", "evidence", "report"]
OUTPUT_KEYS = [
    "schema_version",
    "project",
    "stage_id",
    "execution_mode",
    "migration_status",
    "scientific_status",
    "parity_status",
    "migration_record_sha256",
    "source_artifacts",
]


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _catalog_ids():
    return [x[0] for x in IS_CATALOG]


def load_migration(migration_dir: Path):
    migration_dir = migration_dir.resolve()
    manifest_path = migration_dir / "IS_STAGE_MIGRATION.json"
    table_path = migration_dir / "IS_STAGE_MIGRATION.tsv"
    if not manifest_path.is_file() or not table_path.is_file():
        raise ValueError("Missing IS stage-migration manifest/table")
    manifest = _read_json(manifest_path)
    if manifest.get("project") != "ischemic_stroke":
        raise ValueError("Migration project is not ischemic_stroke")
    if manifest.get("catalog_stage_count") != len(IS_CATALOG):
        raise ValueError("Migration catalog stage count changed")
    if manifest.get("migrated_frozen_stages") != 9:
        raise ValueError("Expected exactly nine frozen-ready IS stages")
    if manifest.get("blocked_stages") != EXPECTED_BLOCKED:
        raise ValueError("Unexpected IS blocked-stage set/order")

    payloads = {}
    for stage_id in _catalog_ids():
        path = migration_dir / "stages" / f"{stage_id}.json"
        if not path.is_file():
            raise ValueError("Missing IS stage migration record: " + str(path))
        payload = _read_json(path)
        if payload.get("project") != "ischemic_stroke" or payload.get("stage_id") != stage_id:
            raise ValueError("Invalid IS stage migration identity: " + stage_id)
        payloads[stage_id] = (path, payload)
    return manifest_path, manifest, payloads


def _verify_artifact(artifact: dict):
    path = Path(artifact["path"])
    if not path.is_file():
        raise ValueError("Missing frozen source artifact: " + str(path))
    size = path.stat().st_size
    if size != int(artifact["size_bytes"]):
        raise ValueError(f"Frozen source size changed: {path}")
    digest = sha(path)
    if digest != artifact["sha256"]:
        raise ValueError(f"Frozen source SHA256 changed: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": size,
        "sha256": digest,
    }


def emit_frozen_stage(stage_json: Path, output_json: Path):
    stage_json = stage_json.resolve()
    payload = _read_json(stage_json)
    if payload.get("project") != "ischemic_stroke":
        raise ValueError("Frozen stage project mismatch")
    if payload.get("migration_status") != READY:
        raise ValueError("Only MIGRATED_FROZEN stages can be emitted")
    if payload.get("scientific_status") != "FROZEN_RESULTS_NOT_RECOMPUTED":
        raise ValueError("Frozen stage scientific status changed")
    artifacts = [_verify_artifact(x) for x in payload.get("source_artifacts", [])]
    if not artifacts:
        raise ValueError("Frozen stage has no source artifacts")
    out = {
        "schema_version": 1,
        "project": "ischemic_stroke",
        "stage_id": payload["stage_id"],
        "execution_mode": "FROZEN_PASSTHROUGH",
        "migration_status": READY,
        "scientific_status": "FROZEN_RESULTS_NOT_RECOMPUTED",
        "parity_status": "SOURCE_IDENTITY_VERIFIED",
        "migration_record_sha256": sha(stage_json),
        "source_artifacts": artifacts,
    }
    atomic_json(output_json, out)
    return out


def bind_frozen_prefix(project_config: Path, migration_dir: Path, output_config: Path):
    project_config = project_config.resolve()
    migration_dir = migration_dir.resolve()
    output_config = output_config.resolve()
    if output_config.exists():
        raise ValueError("Bound project config already exists; refusing overwrite")
    cfg = _read_json(project_config)
    if cfg.get("project") != "ischemic_stroke":
        raise ValueError("Frozen IS binder requires an ischemic_stroke blueprint")
    validate(cfg)

    manifest_path, manifest, payloads = load_migration(migration_dir)
    ready_ids = {
        stage_id
        for stage_id, (_, payload) in payloads.items()
        if payload.get("migration_status") == READY
    }
    expected_ready = set(_catalog_ids()) - set(EXPECTED_BLOCKED)
    if ready_ids != expected_ready:
        raise ValueError("Frozen migration-ready stage set changed")

    policy = cfg.setdefault("analysis_policy", {})
    policy.update({
        "genome_build": "GRCh37",
        "replication_definition": "GIGASTROKE EAS frozen credible-set and direction replication",
        "palindromic_policy": "Frozen legacy canonical allele policy; no re-harmonization in migration",
        "testing_family": "Frozen BBJ primary loci L001-L004",
        "locus_definition": "Frozen GRCh37 BBJ_IS_L001-L004 regions",
        "discovery_ancestry": "Japanese",
        "outcome_validation_ancestry": "EAS",
    })
    cfg["data_review_status"] = "FROZEN_BASELINE_AUDITED"
    cfg["scope"] = "frozen_migration_prefix_bound_human_annotation_pending"
    cfg["migration_registry"] = {
        "path": str(manifest_path),
        "sha256": sha(manifest_path),
        "status": manifest.get("masteromics_binding_status"),
    }

    workspace = Path(cfg["workspace"]).resolve()
    bound = []
    for stage in cfg["stages"]:
        stage_id = stage["id"]
        if stage_id not in ready_ids:
            stage["binding"] = None
            continue
        stage_json, payload = payloads[stage_id]
        source_paths = [str(Path(x["path"]).resolve()) for x in payload["source_artifacts"]]
        output_path = workspace / "results" / "frozen_migration" / f"{stage_id}.json"
        stage["binding"] = {
            "argv": [
                "{python}",
                "-m",
                "masteromics",
                "is-bind",
                "emit",
                str(stage_json),
                "@out0",
            ],
            "inputs": [str(stage_json)] + source_paths,
            "outputs": [{
                "path": str(output_path),
                "kind": "json",
                "keys": OUTPUT_KEYS,
            }],
            "timeout_seconds": 900,
        }
        bound.append(stage_id)

    assessment = inspect(cfg)
    if assessment["unbound_stages"] != EXPECTED_BLOCKED:
        raise ValueError("Bound IS blueprint must leave only final three stages unbound")
    output_config.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output_config, cfg)
    return {
        "project": "ischemic_stroke",
        "bound_stages": bound,
        "bound_stage_count": len(bound),
        "unbound_stages": assessment["unbound_stages"],
        "execution_status": assessment["execution_status"],
        "scientific_status": "FROZEN_RESULTS_NOT_RECOMPUTED",
        "output_config": str(output_config),
    }



def materialize_frozen_prefix(migration_dir: Path, output_dir: Path):
    migration_dir = migration_dir.resolve()
    output_dir = output_dir.resolve()
    _, manifest, payloads = load_migration(migration_dir)
    ready_ids = [
        stage_id for stage_id in _catalog_ids()
        if payloads[stage_id][1].get("migration_status") == READY
    ]
    expected_ready = [x for x in _catalog_ids() if x not in EXPECTED_BLOCKED]
    if ready_ids != expected_ready:
        raise ValueError("Frozen migration-ready stage order changed")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for stage_id in ready_ids:
        target = output_dir / f"{stage_id}.json"
        if target.exists():
            raise ValueError("Frozen passthrough output already exists; refusing overwrite: " + str(target))
        stage_path, _ = payloads[stage_id]
        emitted = emit_frozen_stage(stage_path, target)
        outputs.append({
            "stage_id": stage_id,
            "path": str(target),
            "parity_status": emitted["parity_status"],
        })
    summary = {
        "schema_version": 1,
        "project": "ischemic_stroke",
        "execution_mode": "FROZEN_PREFIX_PROVENANCE_VERIFICATION",
        "scientific_status": "FROZEN_RESULTS_NOT_RECOMPUTED",
        "verified_stage_count": len(outputs),
        "verified_stages": [x["stage_id"] for x in outputs],
        "blocked_stages": EXPECTED_BLOCKED,
        "full_blueprint_status": "BLOCKED_UNBOUND_ADAPTERS",
        "migration_manifest_status": manifest.get("masteromics_binding_status"),
        "outputs": outputs,
    }
    atomic_json(output_dir / "FROZEN_PREFIX_VERIFICATION.json", summary)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("bind")
    p.add_argument("project_config", type=Path)
    p.add_argument("migration_dir", type=Path)
    p.add_argument("output_config", type=Path)
    p = sub.add_parser("emit")
    p.add_argument("stage_json", type=Path)
    p.add_argument("output_json", type=Path)
    p = sub.add_parser("materialize")
    p.add_argument("migration_dir", type=Path)
    p.add_argument("output_dir", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.action == "bind":
            result = bind_frozen_prefix(args.project_config, args.migration_dir, args.output_config)
        elif args.action == "materialize":
            result = materialize_frozen_prefix(args.migration_dir, args.output_dir)
        else:
            result = emit_frozen_stage(args.stage_json, args.output_json)
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print("IS frozen binding error: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
