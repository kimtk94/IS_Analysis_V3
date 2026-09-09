#!/usr/bin/env python3
"""Fail fast when the server cannot safely stage another analysis batch."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

GIB = 1024 ** 3


def dir_size(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file() and not item.is_symlink():
                total += item.stat().st_size
        except FileNotFoundError:
            # A worker may atomically move/delete a staging file while we scan.
            continue
    return total


def gib(value: int) -> float:
    return value / GIB


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", default="/srv/is-analysis")
    parser.add_argument("--staging", default="/srv/is-analysis/data/staging")
    parser.add_argument("--min-free-gib", type=float, default=40.0)
    parser.add_argument("--max-staging-gib", type=float, default=45.0)
    args = parser.parse_args()

    root = Path(args.path)
    staging = Path(args.staging)
    root.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)

    usage = shutil.disk_usage(root)
    staging_bytes = dir_size(staging)

    print(
        f"storage total={gib(usage.total):.1f}GiB "
        f"used={gib(usage.used):.1f}GiB free={gib(usage.free):.1f}GiB "
        f"staging={gib(staging_bytes):.1f}GiB"
    )

    errors: list[str] = []
    if gib(usage.free) < args.min_free_gib:
        errors.append(
            f"free space {gib(usage.free):.1f}GiB is below floor "
            f"{args.min_free_gib:.1f}GiB"
        )
    if gib(staging_bytes) > args.max_staging_gib:
        errors.append(
            f"staging usage {gib(staging_bytes):.1f}GiB exceeds cap "
            f"{args.max_staging_gib:.1f}GiB"
        )

    if errors:
        for error in errors:
            print(f"STORAGE_GUARD_FAIL: {error}", file=sys.stderr)
        return 2

    print("STORAGE_GUARD_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
