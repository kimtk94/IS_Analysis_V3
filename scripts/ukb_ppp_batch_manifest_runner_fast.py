#!/usr/bin/env python3
"""Create paired EUR/EAS batches and optionally execute local preparation.

The runner is conservative: execution requires ``--run`` and source URLs are
recorded but never downloaded by a dry run. Local file paths and file:// URLs are
supported for deterministic testing; HTTP downloads use urllib.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REQUIRED = {"gene", "ancestry", "source_url"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        rows = [{k: (v or "").strip() for k, v in row.items()} for row in reader]
    for row in rows:
        row["gene"] = row["gene"].upper()
        row["ancestry"] = row["ancestry"].upper()
    return rows


def paired_batches(rows: list[dict[str, str]], size: int) -> list[list[str]]:
    represented: dict[str, set[str]] = {}
    for row in rows:
        represented.setdefault(row["gene"], set()).add(row["ancestry"])
    genes = sorted(g for g, ancestries in represented.items() if {"EUR", "EAS"} <= ancestries)
    return [genes[i : i + size] for i in range(0, len(genes), size)]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stage(row: dict[str, str], destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = row["source_url"]
    if url.startswith("file://"):
        shutil.copy2(url[7:], destination)
    elif Path(url).is_file():
        shutil.copy2(url, destination)
    else:
        urllib.request.urlretrieve(url, destination)
    expected = row.get("sha256", "").lower()
    if expected and sha256(destination) != expected:
        destination.unlink(missing_ok=True)
        raise ValueError(f"checksum mismatch for {row['gene']} {row['ancestry']}")
    expected_size = row.get("size_bytes", "")
    if expected_size and destination.stat().st_size != int(expected_size):
        destination.unlink(missing_ok=True)
        raise ValueError(f"size mismatch for {row['gene']} {row['ancestry']}")
    return destination


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser()
    for name in ("base", "qc-dir", "outdir", "standardized-dir", "instrument-dir",
                 "download-manifest", "gene-coordinate-file"):
        p.add_argument(f"--{name}", required=True, type=Path)
    p.add_argument("--batch-size", type=int, default=15)
    p.add_argument("--focus-gene")
    p.add_argument("--focus-max-bytes", type=int, default=20_000_000)
    p.add_argument("--other-max-file-lines", type=int, default=1_000)
    p.add_argument("--p-value-threshold", type=float, default=5e-8)
    p.add_argument("--f-statistic-threshold", type=float, default=10.0)
    p.add_argument("--cis-window-bp", type=int, default=1_000_000)
    p.add_argument("--prepare-script", type=Path, default=Path(__file__).with_name("01_prepare_exposure_fast.R"))
    execution = p.add_mutually_exclusive_group()
    execution.add_argument("--run", action="store_true",
                           help="download the selected batch and run the preparation script")
    execution.add_argument("--download-only", action="store_true",
                           help="download the selected batch without running the preparation script")
    p.add_argument("--stop-on-error", action="store_true")
    args = p.parse_args(argv)
    if args.batch_size < 1 or args.focus_max_bytes < 1 or args.other_max_file_lines < 2:
        p.error("batch size/limits must be positive (line limit must include header and data)")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = read_tsv(args.download_manifest)
    batches = paired_batches(rows, args.batch_size)
    focus = args.focus_gene.upper() if args.focus_gene else None
    selected = [(i + 1, b) for i, b in enumerate(batches) if not focus or focus in b]
    if focus and len(selected) != 1:
        raise SystemExit(f"focus gene {focus} must occur in exactly one paired batch; found {len(selected)}")
    args.qc_dir.mkdir(parents=True, exist_ok=True)
    for child in ("downloads", "processing_logs", "raw_cleanup"):
        (args.qc_dir / child).mkdir(exist_ok=True)
    batch_rows = [{"batch_id": f"batch_{i:03d}", "gene_count": len(b), "genes": ",".join(b)}
                  for i, b in enumerate(batches, 1)]
    write_tsv(args.qc_dir / "batch_manifest.tsv", ["batch_id", "gene_count", "genes"], batch_rows)
    selected_keys = {(gene, anc) for _, genes in selected for gene in genes for anc in ("EUR", "EAS")}
    source_rows = [r for r in rows if (r["gene"], r["ancestry"]) in selected_keys]
    plan = []
    for batch_number, genes in selected:
        for row in source_rows:
            if row["gene"] in genes:
                limit_type = "bytes" if row["gene"] == focus else "lines"
                limit = args.focus_max_bytes if limit_type == "bytes" else args.other_max_file_lines
                plan.append({"batch_id": f"batch_{batch_number:03d}", "gene": row["gene"],
                             "ancestry": row["ancestry"], "source_url": row["source_url"],
                             "limit_type": limit_type, "limit": limit, "status": "planned"})
    write_tsv(args.qc_dir / "execution_plan.tsv",
              ["batch_id", "gene", "ancestry", "source_url", "limit_type", "limit", "status"], plan)
    base_parents = args.base.resolve().parents
    inferred_work_root = base_parents[3] if len(base_parents) > 3 else args.base.resolve().parent
    execute_downloads = args.run or args.download_only
    metadata = {"started_at": datetime.now(timezone.utc).isoformat(), "run": args.run,
                "download_only": args.download_only,
                "batch_size": args.batch_size, "focus_gene": focus,
                "focus_max_bytes": args.focus_max_bytes,
                "other_max_file_lines": args.other_max_file_lines,
                "download_manifest_sha256": sha256(args.download_manifest),
                "gene_coordinates_sha256": sha256(args.gene_coordinate_file),
                "code_root": str(Path.cwd()), "work_root": str(inferred_work_root),
                "p_value_threshold": args.p_value_threshold,
                "f_statistic_threshold": args.f_statistic_threshold,
                "cis_window_bp": args.cis_window_bp,
                "raw_cleanup": False}
    (args.qc_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if not execute_downloads:
        print(f"Dry run: wrote {len(plan)} source tasks to {args.qc_dir}")
        return 0
    failures = 0
    source_lookup = {(row["gene"], row["ancestry"]): row for row in source_rows}
    for task in plan:
        row = source_lookup[(task["gene"], task["ancestry"])]
        suffix = Path(row["source_url"]).name or f"{row['gene']}.tar"
        try:
            archive = stage(row, args.base / row["ancestry"] / suffix)
            if args.download_only:
                task["status"] = "downloaded"
                continue
            command = ["Rscript", str(args.prepare_script), "--archive", str(archive),
                       "--gene", row["gene"], "--ancestry", row["ancestry"],
                       "--coordinates", str(args.gene_coordinate_file),
                       "--standardized-dir", str(args.standardized_dir),
                       "--instrument-dir", str(args.instrument_dir), "--legacy-dir", str(args.outdir),
                       "--limit-type", task["limit_type"], "--limit", str(task["limit"]),
                       "--p-value-threshold", str(args.p_value_threshold),
                       "--f-statistic-threshold", str(args.f_statistic_threshold),
                       "--cis-window-bp", str(args.cis_window_bp)]
            subprocess.run(command, check=True)
            task["status"] = "complete"
        except Exception as exc:
            failures += 1
            task["status"] = f"failed: {exc}"
            if args.stop_on_error:
                break
    write_tsv(args.qc_dir / "batch_progress.tsv",
              ["batch_id", "gene", "ancestry", "source_url", "limit_type", "limit", "status"], plan)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
