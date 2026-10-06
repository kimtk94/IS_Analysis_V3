"""Read-only migration adapter for the frozen ischemic-stroke baseline.

The adapter is configuration-driven and does not rerun stroke analyses. It audits
the existing BBJ/GIGASTROKE and functional-validation artifacts, pins provenance,
and emits baseline artifacts for later explicit MasterOmics stage bindings.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "projects" / "ischemic_stroke.example.json"
DEFAULT_PINS = REPO_ROOT / "projects" / "ischemic_stroke.baseline_pins.json"


def load_config(path: Path):
    cfg = json.loads(path.read_text(encoding="utf-8"))
    if cfg.get("schema_version") != 1:
        raise ValueError("Unsupported IS migration config schema")
    if cfg.get("project") != "ischemic_stroke":
        raise ValueError("IS migration config project must be ischemic_stroke")
    if cfg.get("mode") != "frozen_locus_first_migration":
        raise ValueError("IS migration config must use frozen_locus_first_migration mode")
    design = cfg.get("design", {})
    if design.get("discovery_model") != "BBJ_EAS_LOCUS_FIRST":
        raise ValueError("Current canonical IS design must be BBJ_EAS_LOCUS_FIRST")
    if design.get("discovery_build") != "GRCh37" or design.get("replication_build") != "GRCh37":
        raise ValueError("Frozen BBJ/GIGASTROKE baseline is GRCh37")
    loci = design.get("primary_loci", [])
    genes = design.get("core_genes", [])
    if loci != ["BBJ_IS_L001", "BBJ_IS_L002", "BBJ_IS_L003", "BBJ_IS_L004"]:
        raise ValueError("Canonical IS primary loci contract changed")
    if genes != ["FGF5", "ALDH2", "SH3PXD2A", "COL4A2"]:
        raise ValueError("Canonical IS core mechanism genes contract changed")
    branches = design.get("core_mechanism_branches", {})
    if set(branches) != set(genes):
        raise ValueError("Core mechanism branch map must cover exactly the core genes")
    sources = cfg.get("sources", [])
    if not sources:
        raise ValueError("IS migration config requires source contracts")
    ids = [x.get("id") for x in sources]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate IS source contract id")
    required_source_ids = {
        "bbj_is_canonical",
        "gigastroke_eas_ais",
        "bbj_finemap_regions",
        "bbj_susie_neff_summary",
        "bbj_gigastroke_cs_comparison",
        "molecular_coloc_abf",
        "molecular_coloc_susie",
        "functional_convergence",
        "mouse_functional_freeze",
        "human_vascular_handoff",
    }
    if set(ids) != required_source_ids:
        raise ValueError("Canonical IS source contract set changed")
    for source in sources:
        if not source.get("path") or not source.get("role") or not source.get("required_columns"):
            raise ValueError("Each IS source needs path, role, and required_columns")
        if Path(source["path"]).is_absolute():
            raise ValueError("IS migration source paths must be relative to legacy_root")
    return cfg


def load_pins(path: Path, cfg: dict):
    pins = json.loads(path.read_text(encoding="utf-8"))
    if pins.get("schema_version") != 1 or pins.get("project") != "ischemic_stroke":
        raise ValueError("Unsupported IS baseline pin manifest")
    source_ids = {x["id"] for x in cfg["sources"]}
    pin_ids = set(pins.get("sources", {}))
    if pin_ids != source_ids:
        raise ValueError("IS baseline pin manifest source set changed")
    for source_id, spec in pins["sources"].items():
        if not isinstance(spec.get("size_bytes"), int) or spec["size_bytes"] <= 0:
            raise ValueError("Invalid pinned size for " + source_id)
        digest = spec.get("sha256", "")
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest.lower()):
            raise ValueError("Invalid pinned SHA256 for " + source_id)
    return pins


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _header(path: Path):
    with _open_text(path) as handle:
        first = handle.readline().rstrip("\r\n")
    return first.split("\t") if first else []


def _sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _rows(path: Path):
    with _open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _kv(path: Path):
    return {r["component"]: r["status"] for r in _rows(path)}


def _source_map(cfg):
    return {x["id"]: x for x in cfg["sources"]}


def audit_sources(
    root: Path,
    cfg: dict,
    hash_large: bool = False,
    large_threshold: int = 50 * 1024 * 1024,
    pins: dict | None = None,
):
    records = []
    pin_sources = pins.get("sources", {}) if pins else {}
    for spec in cfg["sources"]:
        path = root / spec["path"]
        if not path.is_file():
            raise ValueError("Missing IS baseline source: " + str(path))
        size = path.stat().st_size
        if size <= 0:
            raise ValueError("Empty IS baseline source: " + str(path))
        header = _header(path)
        required = spec["required_columns"]
        missing = sorted(set(required) - set(header))
        if missing:
            raise ValueError(spec["id"] + " missing required columns: " + ",".join(missing))

        pin = pin_sources.get(spec["id"])
        if pin and size != pin["size_bytes"]:
            raise ValueError(
                f"{spec['id']} size changed: observed={size} pinned={pin['size_bytes']}"
            )
        do_hash = bool(pin) or hash_large or size <= large_threshold
        digest = _sha256(path) if do_hash else ""
        if pin and digest != pin["sha256"]:
            raise ValueError(spec["id"] + " SHA256 changed from frozen baseline")
        checksum_status = (
            "VERIFIED_PIN" if pin else "PINNED" if do_hash else "DEFERRED_LARGE_FILE"
        )
        records.append({
            "dataset_id": spec["id"],
            "role": spec["role"],
            "ancestry": spec.get("ancestry", ""),
            "build": spec.get("build", ""),
            "path": str(path.resolve()),
            "size_bytes": size,
            "sha256": digest,
            "checksum_status": checksum_status,
            "required_columns": ";".join(required),
            "observed_columns": ";".join(header),
        })
    return records


def _write_tsv(path: Path, rows: list[dict]):
    if not rows:
        raise ValueError("Refusing to write empty adapter artifact")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_baseline(
    root: Path,
    outdir: Path,
    config_path: Path = DEFAULT_CONFIG,
    hash_large: bool = False,
    pins_path: Path | None = None,
):
    root = root.resolve()
    outdir = outdir.resolve()
    config_path = config_path.resolve()
    cfg = load_config(config_path)
    pins = load_pins(pins_path.resolve(), cfg) if pins_path else None
    sources = _source_map(cfg)
    source_rows = audit_sources(root, cfg, hash_large=hash_large, pins=pins)
    design = cfg["design"]

    finemap = _rows(root / sources["bbj_finemap_regions"]["path"])
    primary = {r["locus_id"]: r for r in finemap if r.get("role", "").startswith("PRIMARY")}
    missing_loci = [x for x in design["primary_loci"] if x not in primary]
    if missing_loci:
        raise ValueError("Missing primary BBJ loci: " + ",".join(missing_loci))

    finemap_summary = _rows(root / sources["bbj_susie_neff_summary"]["path"])
    finemap_by_locus = {r["locus"]: r for r in finemap_summary}
    missing_finemap = [x for x in design["primary_loci"] if x not in finemap_by_locus]
    if missing_finemap:
        raise ValueError("Missing BBJ SuSiE summaries: " + ",".join(missing_finemap))
    bad_finemap = [
        x for x in design["primary_loci"]
        if finemap_by_locus[x].get("status") != "PASS"
        or finemap_by_locus[x].get("converged", "").upper() != "TRUE"
    ]
    if bad_finemap:
        raise ValueError("Nonconverged BBJ SuSiE loci: " + ",".join(bad_finemap))

    coloc_abf = _rows(root / sources["molecular_coloc_abf"]["path"])
    coloc_susie = _rows(root / sources["molecular_coloc_susie"]["path"])
    abf_core = {
        r.get("gene_symbol") for r in coloc_abf
        if r.get("status") == "PASS" and r.get("gene_symbol")
    }
    susie_core = {r.get("gene_symbol") for r in coloc_susie if r.get("gene_symbol")}
    missing_abf_core = [x for x in design["core_genes"] if x not in abf_core]
    missing_susie_core = [x for x in design["core_genes"] if x not in susie_core]
    if missing_abf_core:
        raise ValueError("Core genes missing from ABF coloc: " + ",".join(missing_abf_core))
    if missing_susie_core:
        raise ValueError("Core genes missing from SuSiE coloc: " + ",".join(missing_susie_core))

    convergence = _rows(root / sources["functional_convergence"]["path"])
    core = {r["gene"]: r for r in convergence if r.get("role") == "CORE"}
    missing_genes = [x for x in design["core_genes"] if x not in core]
    if missing_genes:
        raise ValueError("Missing core IS mechanism genes: " + ",".join(missing_genes))

    evidence_rows = []
    for gene in design["core_genes"]:
        r = core[gene]
        expected_branch = design["core_mechanism_branches"][gene]
        if r.get("mechanism_branch") != expected_branch:
            raise ValueError(f"{gene} mechanism branch changed: {r.get('mechanism_branch')}")
        evidence_rows.append({
            "gene_symbol": gene,
            "locus_id": r["locus"],
            "mechanism_branch": r["mechanism_branch"],
            "best_abf_h4": r.get("best_abf_h4", ""),
            "best_susie_h4": r.get("best_susie_h4", ""),
            "best_abf_dataset": r.get("best_abf_dataset", ""),
            "best_abf_tissue": r.get("best_abf_tissue", ""),
            "interpretation": r.get("interpretation", ""),
            "evidence_status": "FROZEN_BASELINE",
        })

    mouse = _kv(root / sources["mouse_functional_freeze"]["path"])
    handoff = _kv(root / sources["human_vascular_handoff"]["path"])
    expected = cfg["expected_status"]
    if mouse.get("MOUSE_LAYER_CONCLUSION") != expected["mouse_layer_conclusion"]:
        raise ValueError("Unexpected mouse functional-layer conclusion")
    if handoff.get("COLAB_STATUS") != expected["human_vascular_colab_status"]:
        raise ValueError("Human vascular handoff is not Colab-ready")

    coverage_rows = [
        {"stage":"acquisition","status":"FROZEN_ARTIFACT_READY","source_ids":";".join(x["id"] for x in cfg["sources"]),"note":"10 pinned baseline sources are present"},
        {"stage":"source_qc","status":"FROZEN_ARTIFACT_READY","source_ids":"baseline_pins","note":"size/SHA256 verification is enforced by the production adapter"},
        {"stage":"normalize","status":"FROZEN_ARTIFACT_READY","source_ids":"bbj_is_canonical;gigastroke_eas_ais","note":"canonical GRCh37 GWAS inputs are frozen"},
        {"stage":"gwas_loci","status":"FROZEN_ARTIFACT_READY","source_ids":"bbj_finemap_regions","note":"BBJ primary loci L001-L004 are present"},
        {"stage":"finemap","status":"FROZEN_ARTIFACT_READY","source_ids":"bbj_susie_neff_summary","note":"all four primary BBJ SuSiE fits are PASS and converged"},
        {"stage":"cross_ancestry","status":"FROZEN_ARTIFACT_READY","source_ids":"bbj_gigastroke_cs_comparison","note":"BBJ-GIGASTROKE credible-set comparison is frozen"},
        {"stage":"molecular_coloc","status":"FROZEN_ARTIFACT_READY","source_ids":"molecular_coloc_abf;molecular_coloc_susie","note":"all four core genes are represented in ABF and SuSiE molecular evidence"},
        {"stage":"mechanism","status":"FROZEN_ARTIFACT_READY","source_ids":"functional_convergence","note":"four core mechanism branches are frozen"},
        {"stage":"celltype","status":"FROZEN_ARTIFACT_READY","source_ids":"mouse_functional_freeze","note":"mouse cell-type layer is frozen with stroke-DGE limitation explicit"},
        {"stage":"human_annotation","status":"PENDING_AUTHOR_ANNOTATION","source_ids":"human_vascular_handoff","note":"GSE256493 RDS handoff is READY, but author annotation is not yet frozen"},
        {"stage":"evidence","status":"BLOCKED_UPSTREAM","source_ids":"","note":"requires completed human_annotation"},
        {"stage":"report","status":"BLOCKED_UPSTREAM","source_ids":"","note":"requires completed integrated evidence"},
    ]

    outdir.mkdir(parents=True, exist_ok=True)
    source_out = outdir / "IS_BASELINE_SOURCE_AUDIT.tsv"
    evidence_out = outdir / "IS_BASELINE_CORE_EVIDENCE.tsv"
    coverage_out = outdir / "IS_STAGE_COVERAGE.tsv"
    status_out = outdir / "IS_BASELINE_STATUS.json"
    _write_tsv(source_out, source_rows)
    _write_tsv(evidence_out, evidence_rows)
    _write_tsv(coverage_out, coverage_rows)

    checksum_complete = all(r["checksum_status"] in {"PINNED", "VERIFIED_PIN"} for r in source_rows)
    frozen_ready = sum(r["status"] == "FROZEN_ARTIFACT_READY" for r in coverage_rows)
    pending = [r["stage"] for r in coverage_rows if r["status"] != "FROZEN_ARTIFACT_READY"]
    status = {
        "schema_version": 3,
        "project": cfg["project"],
        "mode": cfg["mode"],
        "adapter_status": "PASS",
        "scientific_status": cfg["scientific_status"],
        "masteromics_binding_status": cfg["execution_status"],
        "checksum_status": "COMPLETE" if checksum_complete else "PARTIAL_LARGE_FILES_DEFERRED",
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "pin_manifest_path": str(pins_path.resolve()) if pins_path else None,
        "pin_manifest_sha256": _sha256(pins_path.resolve()) if pins_path else None,
        "source_count": len(source_rows),
        "stage_coverage": {
            "total": len(coverage_rows),
            "frozen_artifact_ready": frozen_ready,
            "not_ready_stages": pending,
        },
        "discovery_model": design["discovery_model"],
        "discovery_ancestry": design["discovery_ancestry"],
        "discovery_build": design["discovery_build"],
        "replication_model": design["replication_model"],
        "replication_ancestry": design["replication_ancestry"],
        "replication_build": design["replication_build"],
        "primary_loci": design["primary_loci"],
        "core_genes": design["core_genes"],
        "mouse_layer_conclusion": mouse.get("MOUSE_LAYER_CONCLUSION"),
        "human_vascular_colab_status": handoff.get("COLAB_STATUS"),
        "human_vascular_remote_path": handoff.get("REMOTE_PATH"),
        "artifacts": {
            "source_audit": str(source_out),
            "core_evidence": str(evidence_out),
            "stage_coverage": str(coverage_out),
        },
        "next_gate": cfg["next_gate"],
    }
    status_out.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy_root", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--pins", type=Path, default=DEFAULT_PINS)
    parser.add_argument("--hash-large", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = build_baseline(
            args.legacy_root,
            args.outdir,
            config_path=args.config,
            hash_large=args.hash_large,
            pins_path=args.pins,
        )
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, ValueError, KeyError, csv.Error, json.JSONDecodeError) as exc:
        print("IS adapter error: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
