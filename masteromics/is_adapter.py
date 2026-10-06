"""Read-only adapter for the frozen ischemic-stroke analysis baseline.

This module does not rerun stroke analyses. It audits the existing BBJ/GIGASTROKE
and functional-validation artifacts, pins their provenance, and emits canonical
baseline audit/evidence artifacts for later MasterOmics stage bindings.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

SOURCE_CONTRACTS = [
    {
        "id": "bbj_is_canonical",
        "rel": "data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz",
        "role": "discovery_gwas",
        "ancestry": "Japanese",
        "build": "GRCh37",
        "required": ["dataset","phenotype","build","chr","pos","effect_allele","other_allele","beta","se","p","eaf","n","variant_id"],
    },
    {
        "id": "gigastroke_eas_ais",
        "rel": "data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz",
        "role": "replication_gwas",
        "ancestry": "EAS",
        "build": "GRCh37",
        "required": ["dataset","phenotype","ancestry","build","chr","pos","effect_allele","other_allele","beta","se","p","eaf","variant_id"],
    },
    {
        "id": "bbj_finemap_regions",
        "rel": "results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv",
        "role": "regional_contract",
        "required": ["locus_id","chr","region_start","region_end","lead_variant","role"],
    },
    {
        "id": "bbj_gigastroke_cs_comparison",
        "rel": "results/is/stage4_cross_eas/fine_mapping/BBJ_GIGASTROKE_CS_COMPARISON.tsv",
        "role": "cross_ancestry_finemap",
        "required": ["phenotype","locus","bbj_top_variant","bbj_top_pip","giga_top_variant","giga_top_pip","cs_jaccard"],
    },
    {
        "id": "functional_convergence",
        "rel": "results/is/stage5_functional/phase9d_literature_benchmark/IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv",
        "role": "functional_evidence",
        "required": ["gene","locus","role","mechanism_branch","best_abf_h4","best_susie_h4","interpretation"],
    },
    {
        "id": "mouse_functional_freeze",
        "rel": "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv",
        "role": "single_cell_freeze",
        "required": ["component","status"],
    },
    {
        "id": "human_vascular_handoff",
        "rel": "results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1/GSE256493_DRIVE_HANDOFF.tsv",
        "role": "human_single_cell_handoff",
        "required": ["component","status"],
    },
]

CORE_GENES = ["FGF5","ALDH2","SH3PXD2A","COL4A2"]
PRIMARY_LOCI = ["BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L003","BBJ_IS_L004"]


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


def audit_sources(root: Path, hash_large: bool = False, large_threshold: int = 50 * 1024 * 1024):
    records = []
    for spec in SOURCE_CONTRACTS:
        path = root / spec["rel"]
        if not path.is_file():
            raise ValueError("Missing IS baseline source: " + str(path))
        if path.stat().st_size <= 0:
            raise ValueError("Empty IS baseline source: " + str(path))
        header = _header(path)
        missing = sorted(set(spec["required"]) - set(header))
        if missing:
            raise ValueError(spec["id"] + " missing required columns: " + ",".join(missing))
        do_hash = hash_large or path.stat().st_size <= large_threshold
        records.append({
            "dataset_id": spec["id"],
            "role": spec["role"],
            "ancestry": spec.get("ancestry", ""),
            "build": spec.get("build", ""),
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path) if do_hash else "",
            "checksum_status": "PINNED" if do_hash else "DEFERRED_LARGE_FILE",
            "required_columns": ";".join(spec["required"]),
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


def build_baseline(root: Path, outdir: Path, hash_large: bool = False):
    root = root.resolve()
    outdir = outdir.resolve()
    source_rows = audit_sources(root, hash_large=hash_large)

    finemap_path = root / next(s["rel"] for s in SOURCE_CONTRACTS if s["id"] == "bbj_finemap_regions")
    finemap = _rows(finemap_path)
    primary = {r["locus_id"]: r for r in finemap if r.get("role", "").startswith("PRIMARY")}
    missing_loci = [x for x in PRIMARY_LOCI if x not in primary]
    if missing_loci:
        raise ValueError("Missing primary BBJ loci: " + ",".join(missing_loci))

    convergence_path = root / next(s["rel"] for s in SOURCE_CONTRACTS if s["id"] == "functional_convergence")
    convergence = _rows(convergence_path)
    core = {r["gene"]: r for r in convergence if r.get("role") == "CORE"}
    missing_genes = [x for x in CORE_GENES if x not in core]
    if missing_genes:
        raise ValueError("Missing core IS mechanism genes: " + ",".join(missing_genes))

    evidence_rows = []
    for gene in CORE_GENES:
        r = core[gene]
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

    mouse_path = root / next(s["rel"] for s in SOURCE_CONTRACTS if s["id"] == "mouse_functional_freeze")
    handoff_path = root / next(s["rel"] for s in SOURCE_CONTRACTS if s["id"] == "human_vascular_handoff")
    mouse = _kv(mouse_path)
    handoff = _kv(handoff_path)
    if mouse.get("MOUSE_LAYER_CONCLUSION") != "CELLTYPE_LOCALIZATION_SUPPORTED_TARGET_SPECIFIC_STROKE_DGE_NOT_ESTABLISHED":
        raise ValueError("Unexpected mouse functional-layer conclusion")
    if handoff.get("COLAB_STATUS") != "READY":
        raise ValueError("Human vascular handoff is not Colab-ready")

    outdir.mkdir(parents=True, exist_ok=True)
    source_out = outdir / "IS_BASELINE_SOURCE_AUDIT.tsv"
    evidence_out = outdir / "IS_BASELINE_CORE_EVIDENCE.tsv"
    status_out = outdir / "IS_BASELINE_STATUS.json"
    _write_tsv(source_out, source_rows)
    _write_tsv(evidence_out, evidence_rows)

    checksum_complete = all(r["checksum_status"] == "PINNED" for r in source_rows)
    status = {
        "schema_version": 1,
        "project": "ischemic_stroke",
        "adapter_status": "PASS",
        "scientific_status": "FROZEN_RESULTS_NOT_RECOMPUTED",
        "masteromics_binding_status": "ADAPTER_AUDIT_ONLY",
        "checksum_status": "COMPLETE" if checksum_complete else "PARTIAL_LARGE_FILES_DEFERRED",
        "source_count": len(source_rows),
        "primary_loci": PRIMARY_LOCI,
        "core_genes": CORE_GENES,
        "mouse_layer_conclusion": mouse.get("MOUSE_LAYER_CONCLUSION"),
        "human_vascular_colab_status": handoff.get("COLAB_STATUS"),
        "human_vascular_remote_path": handoff.get("REMOTE_PATH"),
        "artifacts": {
            "source_audit": str(source_out),
            "core_evidence": str(evidence_out),
        },
        "next_gate": "IMPLEMENT_EXPLICIT_MASTEROMICS_IS_STAGE_BINDINGS",
    }
    status_out.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy_root", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--hash-large", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = build_baseline(args.legacy_root, args.outdir, hash_large=args.hash_large)
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        print("IS adapter error: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
