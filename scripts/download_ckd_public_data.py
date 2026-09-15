#!/usr/bin/env python3
"""Download and verify public CKD thesis datasets from the reviewed manifest.

The downloader deliberately separates public bulk GWAS downloads from restricted
KoGES data and from UKB-PPP full per-protein archives. UKB-PPP bulk cis-region
files should only be fetched for proteins that survive the initial instrument MR
screen, using data/metadata/ukb_ppp_download_manifest.tsv.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data/metadata/ckd_public_download_manifest.tsv"


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def gzip_test(path: Path) -> None:
    if not path.name.endswith(".gz"):
        return
    with gzip.open(path, "rb") as fh:
        while fh.read(1024 * 1024):
            pass


def curl_download(url: str, output: Path) -> None:
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError("curl is required for resumable production downloads")
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [curl, "-L", "--fail", "--retry", "5", "--retry-all-errors", "-C", "-", "-o", str(output), url]
    subprocess.run(cmd, check=True)


def resolve_figshare(article_url: str) -> list[dict]:
    req = urllib.request.Request(article_url, headers={"User-Agent": "IS-Analysis-V3/CKD-data-acquisition"})
    with urllib.request.urlopen(req, timeout=60) as response:
        payload = json.load(response)
    files = payload.get("files", [])
    return [{
        "name": f.get("name", ""),
        "download_url": f.get("download_url", ""),
        "size": f.get("size"),
        "id": f.get("id"),
        "computed_md5": f.get("computed_md5", ""),
    } for f in files]


def write_metadata(row: dict[str, str], output: Path, resolved_url: str) -> None:
    stat = output.stat()
    meta = {
        "data_id": row["data_id"],
        "role": row["role"],
        "ancestry": row["ancestry"],
        "phenotype": row["phenotype"],
        "source": row["source"],
        "source_url": resolved_url,
        "file": str(output),
        "bytes": stat.st_size,
        "sha256": sha256(output),
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "genome_build": row["genome_build"],
        "sample_n": row["sample_n"],
    }
    output.with_suffix(output.suffix + ".download.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )


def select_rows(rows: list[dict[str, str]], data_ids: list[str], public_core: bool) -> list[dict[str, str]]:
    if data_ids:
        wanted = set(data_ids)
        selected = [r for r in rows if r["data_id"] in wanted]
        missing = wanted - {r["data_id"] for r in selected}
        if missing:
            raise SystemExit(f"Unknown data_id(s): {', '.join(sorted(missing))}")
        return selected
    if public_core:
        return [r for r in rows if r["access"] == "public" and r["priority"] in {"1", "2"}]
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--data-id", action="append", default=[])
    parser.add_argument("--public-core", action="store_true", help="priority 1/2 public rows only")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--figshare-list", action="store_true", help="resolve Figshare article files and print them")
    parser.add_argument("--figshare-file-regex", default="", help="regex selecting Figshare file name(s) to download")
    args = parser.parse_args()

    rows = select_rows(load_manifest(args.manifest), args.data_id, args.public_core)
    if args.list:
        for r in rows:
            print("\t".join([r["data_id"], r["role"], r["ancestry"], r["phenotype"], r["access"], r["url"]]))
        return 0

    failures = []
    for row in rows:
        data_id = row["data_id"]
        access = row["access"]
        if access != "public":
            print(f"SKIP {data_id}: access={access}")
            continue
        source_dir = args.dest / row["role"] / row["ancestry"].replace("/", "_") / data_id

        try:
            if row["source"] == "Figshare" and row["file_name"] == "AUTO":
                files = resolve_figshare(row["url"])
                print(f"FIGSHARE {data_id}: {len(files)} file(s)")
                for f in files:
                    print(f"  {f['id']}\t{f['size']}\t{f['name']}\t{f['download_url']}")
                if args.figshare_list or not args.figshare_file_regex:
                    continue
                rx = re.compile(args.figshare_file_regex, re.I)
                matches = [f for f in files if rx.search(f["name"])]
                if not matches:
                    raise RuntimeError(f"No Figshare files match: {args.figshare_file_regex}")
                for f in matches:
                    out = source_dir / f["name"]
                    print(f"DOWNLOAD {data_id}: {f['name']}")
                    if not args.dry_run:
                        curl_download(f["download_url"], out)
                        gzip_test(out)
                        write_metadata(row, out, f["download_url"])
                continue

            out = source_dir / row["file_name"]
            print(f"DOWNLOAD {data_id}: {row['url']} -> {out}")
            if args.dry_run:
                continue
            curl_download(row["url"], out)
            gzip_test(out)
            write_metadata(row, out, row["url"])
            print(f"OK {data_id}: {out.stat().st_size} bytes sha256={sha256(out)}")
        except Exception as exc:  # continue to report all acquisition failures
            failures.append((data_id, str(exc)))
            print(f"FAIL {data_id}: {exc}", file=sys.stderr)

    if failures:
        print("\nFAILED DATASETS", file=sys.stderr)
        for data_id, error in failures:
            print(f"- {data_id}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
