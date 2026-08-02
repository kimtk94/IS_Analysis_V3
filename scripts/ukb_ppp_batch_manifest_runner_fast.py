#!/usr/bin/env python3
"""Create paired EUR/EAS batches and optionally execute local preparation.

The runner is conservative: execution requires ``--run`` and source URLs are
recorded but never downloaded by a dry run. Local file paths and file:// URLs are
supported for deterministic testing; HTTP downloads use urllib.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

MANIFEST_ALIASES = {
    "gene": ("gene", "gene_symbol"),
    "ancestry": ("ancestry",),
    "source_url": ("source_url", "url"),
    "source_file": ("source_file",),
    "size_bytes": ("size_bytes", "expected_size_bytes"),
    "sha256": ("sha256",),
    "md5": ("md5",),
    "synapse_id": ("synapse_id",),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = {canonical for canonical in ("gene", "ancestry", "source_url")
                   if not any(alias in fields for alias in MANIFEST_ALIASES[canonical])}
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        rows = []
        for raw in reader:
            raw = {k: (v or "").strip() for k, v in raw.items()}
            row = dict(raw)
            for canonical, aliases in MANIFEST_ALIASES.items():
                row[canonical] = next((raw[name] for name in aliases if raw.get(name)), "")
            rows.append(row)
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


def add_source_keys(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Return manifest rows with a stable key that is unique within the manifest."""
    keyed_rows = []
    occurrences: dict[str, int] = {}
    for row in rows:
        if row.get("synapse_id"):
            base_key = f"synapse:{row['synapse_id']}"
        else:
            identity = "\0".join(row.get(field, "") for field in
                                  ("gene", "ancestry", "source_file", "source_url"))
            base_key = f"source:{hashlib.sha256(identity.encode()).hexdigest()}"
        occurrences[base_key] = occurrences.get(base_key, 0) + 1
        source_key = base_key
        if occurrences[base_key] > 1:
            source_key = f"{base_key}#{occurrences[base_key]}"
        keyed_rows.append({**row, "source_key": source_key})
    return keyed_rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_digests(path: Path, names: tuple[str, ...] | list[str]) -> dict[str, str]:
    """Calculate the requested supported digests in a single read pass."""
    unsupported = set(names) - {"sha256", "md5"}
    if unsupported:
        raise ValueError(f"unsupported digest(s): {', '.join(sorted(unsupported))}")
    digests = {
        name: hashlib.sha256() if name == "sha256" else hashlib.md5(usedforsecurity=False)
        for name in dict.fromkeys(names)
    }
    if not digests:
        return {}
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            for digest in digests.values():
                digest.update(chunk)
    return {name: digest.hexdigest() for name, digest in digests.items()}


def sha256(path: Path) -> str:
    return file_digests(path, ["sha256"])["sha256"]


def md5(path: Path) -> str:
    return file_digests(path, ["md5"])["md5"]


ARCHIVE_SUFFIXES = (".tar.gz", ".tar.bz2", ".tar.xz", ".tgz", ".tbz2", ".txz",
                    ".zip", ".tar", ".gz", ".bz2", ".xz")
DATA_SUFFIXES = (".tsv", ".txt", ".csv")


def safe_source_stem(name: str) -> str:
    """Return a filesystem-safe assay identifier, including for compound archives."""
    stem = Path(name).name
    lowered = stem.lower()
    for suffix in ARCHIVE_SUFFIXES:
        if lowered.endswith(suffix):
            stem = stem[:-len(suffix)]
            lowered = lowered[:-len(suffix)]
            break
    for suffix in DATA_SUFFIXES:
        if lowered.endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._-")


def test_sample_name(row: dict[str, str]) -> str:
    """Build a sample name that preserves the manifest source identity."""
    source_name = row.get("source_file", "") or Path(row.get("source_url", "")).name
    identifier = safe_source_stem(source_name) if source_name else ""
    if not identifier:
        identifier = safe_source_stem(row.get("synapse_id", ""))
    if not identifier:
        raise ValueError(f"cannot derive sample name for {row['gene']} {row['ancestry']}")
    return f"{identifier}.tsv"


def validation_error(row: dict[str, str], path: Path) -> str | None:
    """Return the first manifest validation failure for *path*, if any."""
    expected_size = row.get("size_bytes", "")
    if expected_size and path.stat().st_size != int(expected_size):
        return "size"
    expected = {name: row.get(name, "").lower() for name in ("sha256", "md5")
                if row.get(name, "")}
    actual = file_digests(path, list(expected))
    if expected.get("sha256") and actual["sha256"] != expected["sha256"]:
        return "checksum"
    if expected.get("md5") and actual["md5"] != expected["md5"]:
        return "MD5"
    return None


def stage(row: dict[str, str], destination: Path, *, reuse_unverified: bool = False) -> Path:
    """Download a source, reusing an existing destination only when it is safe.

    Existing regular files must pass every validation value supplied by the
    manifest.  With no supplied values, the safe default is to redownload;
    callers must explicitly set ``reuse_unverified`` to reuse such a file.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    validations = [row.get(name, "") for name in ("size_bytes", "sha256", "md5")]
    if destination.is_file():
        if any(validations):
            failure = validation_error(row, destination)
            if failure is None:
                return destination
            print(f"Existing destination {failure} mismatch; redownloading {destination}",
                  file=sys.stderr)
        elif reuse_unverified:
            return destination
        else:
            print(f"Existing destination has no manifest validation; redownloading {destination}",
                  file=sys.stderr)
        destination.unlink()

    url = row["source_url"]
    if url.startswith("file://"):
        shutil.copy2(url[7:], destination)
    elif Path(url).is_file():
        shutil.copy2(url, destination)
    elif row.get("synapse_id"):
        if not shutil.which("synapse"):
            raise RuntimeError("Synapse CLI is required: install synapseclient and run 'synapse login'")
        subprocess.run(["synapse", "get", row["synapse_id"],
                        "--downloadLocation", str(destination.parent)], check=True)
        downloaded = destination.parent / (row.get("source_file") or destination.name)
        if downloaded != destination:
            shutil.move(downloaded, destination)
        if not destination.is_file():
            raise FileNotFoundError(f"Synapse download did not create {destination}")
    else:
        urllib.request.urlretrieve(url, destination)
    failure = validation_error(row, destination)
    if failure:
        destination.unlink(missing_ok=True)
        raise ValueError(f"{failure} mismatch for {row['gene']} {row['ancestry']}")
    return destination


def write_test_sample(source: Path, destination: Path, limit_type: str, limit: int) -> Path:
    """Write a bounded, uncompressed sample from a plain or archived text file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with contextlib.ExitStack() as stack:
        if tarfile.is_tarfile(source):
            archive = stack.enter_context(tarfile.open(source, mode="r:*"))
            member = next((item for item in archive if item.isfile()), None)
            if member is None:
                raise ValueError(f"{source}: archive contains no regular member")
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f"{source}: cannot read archive member {member.name}")
            handle = stack.enter_context(stream)
            if member.name.lower().endswith(".gz"):
                handle = stack.enter_context(gzip.GzipFile(fileobj=handle, mode="rb"))
        else:
            with source.open("rb") as probe:
                gzip_magic = probe.read(2) == b"\x1f\x8b"
            if gzip_magic or source.name.lower().endswith(".gz"):
                handle = stack.enter_context(gzip.open(source, "rb"))
            else:
                handle = stack.enter_context(source.open("rb"))

        try:
            output = stack.enter_context(destination.open("xb"))
        except FileExistsError as exc:
            raise FileExistsError(
                f"refusing to overwrite existing test sample: {destination}"
            ) from exc

        complete_lines = 0
        bytes_written = 0
        try:
            if limit_type == "lines":
                for _ in range(limit):
                    line = handle.readline()
                    if not line:
                        break
                    output.write(line)
                    complete_lines += 1
                    bytes_written += len(line)
            else:
                while bytes_written < limit:
                    line = handle.readline(limit - bytes_written + 1)
                    if not line or len(line) > limit - bytes_written:
                        break
                    output.write(line)
                    complete_lines += 1
                    bytes_written += len(line)

            if complete_lines < 2:
                raise ValueError(f"{source}: sample contains no complete data rows")
            assert bytes_written <= limit or limit_type == "lines"
        except Exception:
            output.close()
            destination.unlink(missing_ok=True)
            raise
    return destination


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser()
    for name in ("base", "qc-dir", "outdir", "standardized-dir", "instrument-dir",
                 "download-manifest", "gene-coordinate-file"):
        p.add_argument(f"--{name}", required=True, type=Path)
    p.add_argument("--batch-size", type=int, default=15)
    p.add_argument("--focus-gene")
    focus_limit = p.add_mutually_exclusive_group()
    focus_limit.add_argument(
        "--focus-max-file-lines", type=int, default=500_000,
        help="maximum decompressed lines for the focus gene, including the header",
    )
    focus_limit.add_argument(
        "--focus-max-bytes", type=int,
        help="legacy alternative: maximum decompressed bytes for the focus gene",
    )
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
    execution.add_argument("--test-data-only", action="store_true",
                           help="reuse existing raw files when available, otherwise download the "
                                "selected batch, and write bounded TSV test samples")
    p.add_argument("--stop-on-error", action="store_true")
    p.add_argument("--reuse-unverified", action="store_true",
                   help="reuse existing downloads when the manifest has no size or checksum")
    args = p.parse_args(argv)
    focus_limit_value = (args.focus_max_bytes if args.focus_max_bytes is not None
                         else args.focus_max_file_lines)
    if args.batch_size < 1 or focus_limit_value < 2 or args.other_max_file_lines < 2:
        p.error("batch size/limits must be positive (line limit must include header and data)")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = add_source_keys(read_tsv(args.download_manifest))
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
    selected_genes = {gene for _, genes in selected for gene in genes}
    source_rows = [row for row in rows if row["gene"] in selected_genes]
    plan = []
    for batch_number, genes in selected:
        for row in source_rows:
            if row["gene"] in genes:
                is_focus = row["gene"] == focus
                limit_type = "bytes" if is_focus and args.focus_max_bytes is not None else "lines"
                if is_focus:
                    limit = (args.focus_max_bytes if args.focus_max_bytes is not None
                             else args.focus_max_file_lines)
                else:
                    limit = args.other_max_file_lines
                plan.append({"batch_id": f"batch_{batch_number:03d}", "gene": row["gene"],
                             "ancestry": row["ancestry"], "source_key": row["source_key"],
                             "source_file": row.get("source_file", ""),
                             "synapse_id": row.get("synapse_id", ""),
                             "source_url": row["source_url"],
                             "size_bytes": row.get("size_bytes", ""),
                             "sha256": row.get("sha256", ""), "md5": row.get("md5", ""),
                             "limit_type": limit_type, "limit": limit, "status": "planned"})
    plan_fields = ["batch_id", "gene", "ancestry", "source_key", "source_file",
                   "synapse_id", "source_url", "size_bytes", "sha256", "md5",
                   "limit_type", "limit", "status"]
    write_tsv(args.qc_dir / "execution_plan.tsv",
              plan_fields, plan)
    base_parents = args.base.resolve().parents
    inferred_work_root = base_parents[3] if len(base_parents) > 3 else args.base.resolve().parent
    execute_downloads = args.run or args.download_only or args.test_data_only
    metadata = {"started_at": datetime.now(timezone.utc).isoformat(), "run": args.run,
                "download_only": args.download_only,
                "test_data_only": args.test_data_only,
                "reuse_unverified": args.reuse_unverified,
                "batch_size": args.batch_size, "focus_gene": focus,
                "focus_max_file_lines": args.focus_max_file_lines,
                "focus_max_bytes": args.focus_max_bytes,
                "other_max_file_lines": args.other_max_file_lines,
                "download_manifest_sha256": file_digests(
                    args.download_manifest, ["sha256"])["sha256"],
                "gene_coordinates_sha256": file_digests(
                    args.gene_coordinate_file, ["sha256"])["sha256"],
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
    test_sources: dict[str, Path] = {}
    test_samples: dict[str, Path] = {}
    if args.test_data_only:
        test_samples = {
            task["source_key"]: (args.outdir / "test_data" / task["ancestry"] /
                                 test_sample_name(task))
            for task in plan
        }
        if all(sample.is_file() for sample in test_samples.values()):
            for task in plan:
                task["status"] = "sampled"
            write_tsv(args.qc_dir / "batch_progress.tsv", plan_fields, plan)
            return 0

        # Complete the raw-file preflight for the whole selected batch before
        # writing any missing samples.  This prevents a failed late download
        # from creating additional test data for only part of a batch.
        for task in plan:
            suffix = (task.get("source_file") or Path(task["source_url"]).name
                      or f"{task['gene']}.tar")
            try:
                test_sources[task["source_key"]] = stage(
                    task,
                    args.base / task["ancestry"] / suffix,
                    reuse_unverified=True,
                )
                task["status"] = "raw_ready"
            except Exception as exc:
                failures += 1
                task["status"] = f"failed: {exc}"
                if args.stop_on_error:
                    break
        if failures:
            write_tsv(args.qc_dir / "batch_progress.tsv", plan_fields, plan)
            return 1

    for task in plan:
        # The plan retains all source metadata, so execution cannot select a
        # different assay when gene and ancestry are shared by multiple rows.
        row = task
        suffix = row.get("source_file") or Path(row["source_url"]).name or f"{row['gene']}.tar"
        try:
            if args.test_data_only:
                archive = test_sources[row["source_key"]]
            else:
                archive = stage(
                    row,
                    args.base / row["ancestry"] / suffix,
                    reuse_unverified=args.reuse_unverified,
                )
            if args.download_only:
                task["status"] = "downloaded"
                continue
            if args.test_data_only:
                sample = test_samples[row["source_key"]]
                if not sample.is_file():
                    write_test_sample(archive, sample, task["limit_type"], int(task["limit"]))
                task["status"] = "sampled"
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
              plan_fields, plan)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
