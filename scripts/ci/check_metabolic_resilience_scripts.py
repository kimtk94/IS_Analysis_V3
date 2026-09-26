#!/usr/bin/env python3
"""Static and scientific-contract CI checks for metabolic-resilience code."""

from __future__ import annotations

import ast
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT_DIR = ROOT / "scripts" / "metabolic_resilience"

ERREXIT_RE = re.compile(
    r"(?m)^\s*set\s+(?:-[A-Za-z]*e[A-Za-z]*|-o\s+errexit)\b"
)
SECRET_RE = re.compile(
    r"(?i)(?:api[_-]?key|secret|token|password)\s*=\s*['\"][^'\"]{8,}['\"]"
)
HARDCODED_ROOT_RE = re.compile(
    r"(?m)^\s*ROOT=['\"]/srv/is-analysis['\"]\s*$"
)
PY_START_RE = re.compile(
    r"^\s*python3?\s+-\s+<<['\"]?([A-Za-z_][A-Za-z0-9_]*)['\"]?\s*$"
)


def check_bash_syntax(path: pathlib.Path) -> list[str]:
    proc = subprocess.run(
        ["bash", "-n", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode == 0:
        return []
    return [f"{path}: bash -n failed: {proc.stderr.strip()}"]


def extract_python_heredocs(text: str) -> list[tuple[int, str]]:
    lines = text.splitlines()
    blocks: list[tuple[int, str]] = []
    i = 0

    while i < len(lines):
        m = PY_START_RE.match(lines[i])
        if not m:
            i += 1
            continue

        delimiter = m.group(1)
        start = i + 2
        buf: list[str] = []
        i += 1

        while i < len(lines) and lines[i].strip() != delimiter:
            buf.append(lines[i])
            i += 1

        blocks.append((start, "\n".join(buf) + "\n"))
        i += 1

    return blocks


def check_embedded_python(path: pathlib.Path, text: str) -> list[str]:
    errors: list[str] = []

    for start_line, code in extract_python_heredocs(text):
        try:
            ast.parse(code, filename=f"{path}:heredoc@{start_line}")
        except SyntaxError as exc:
            line = start_line + (exc.lineno or 1) - 1
            errors.append(
                f"{path}: embedded Python syntax error near shell line "
                f"{line}: {exc.msg}"
            )

    return errors


def check_python(path: pathlib.Path) -> list[str]:
    try:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        return []
    except SyntaxError as exc:
        return [f"{path}:{exc.lineno}: Python syntax error: {exc.msg}"]


def check_r(path: pathlib.Path) -> list[str]:
    if not shutil.which("Rscript"):
        return [f"{path}: Rscript unavailable on CI runner"]

    proc = subprocess.run(
        ["Rscript", "-e", f"parse(file={str(path)!r})"],
        text=True,
        capture_output=True,
        check=False,
    )

    if proc.returncode == 0:
        return []

    return [f"{path}: R parse failed: {proc.stderr.strip()}"]


def policy_checks(path: pathlib.Path, text: str) -> list[str]:
    errors: list[str] = []

    if path.suffix == ".sh" and ERREXIT_RE.search(text):
        errors.append(f"{path}: prohibited errexit/set -e usage")

    if SECRET_RE.search(text):
        errors.append(f"{path}: possible hard-coded credential literal")

    if HARDCODED_ROOT_RE.search(text):
        errors.append(
            f"{path}: hard-coded ROOT=/srv/is-analysis; "
            "use configurable IS_ANALYSIS_ROOT with /srv/is-analysis as default"
        )

    return errors


def require(path: pathlib.Path, pattern: str, description: str) -> list[str]:
    if not path.exists():
        return [f"missing scientific-contract file: {path}"]

    text = path.read_text(encoding="utf-8", errors="replace")

    if re.search(pattern, text, flags=re.I | re.M | re.S):
        return []

    return [f"{path}: missing locked contract: {description}"]


def scientific_contract_checks() -> list[str]:
    v2 = SCRIPT_DIR / "drive_v2"
    errors: list[str] = []

    errors += require(
        v2 / "43_stage3d3_mr.R",
        r"if\s*\(\s*k\s*>=\s*3\s*\).*?WEIGHTED_MEDIAN",
        "weighted median only when k >= 3",
    )
    errors += require(
        v2 / "43_stage3d3_mr.R",
        r"if\s*\(\s*k\s*>=\s*10\s*\).*?MR_EGGER",
        "MR-Egger only when k >= 10",
    )
    errors += require(
        v2 / "50_stage3e0_coloc_preflight.sh",
        r"T2D.*HOLD.*case proportion",
        "T2D coloc HOLD until verified case proportion",
    )
    errors += require(
        v2 / "52_stage3e2_coloc_abf.R",
        r"1e-6.*1e-5.*1e-4",
        "coloc p12 sensitivity at 1e-6, 1e-5, 1e-4",
    )
    errors += require(
        v2 / "32_stage3c1_ld_clump.py",
        r"0\.01",
        "LD clumping r2 threshold 0.01",
    )
    errors += require(
        v2 / "31_stage3c0r_build_1000g_eur_references.sh",
        r"integrated_v5b",
        "1000G Phase 3 v5b reference",
    )

    return errors


def main() -> int:
    files = sorted(
        p for p in SCRIPT_DIR.rglob("*")
        if p.is_file() and p.suffix in {".sh", ".py", ".R"}
    )

    if not files:
        print(f"[FAIL] no source files found under {SCRIPT_DIR}", file=sys.stderr)
        return 1

    failures: list[str] = []

    print(f"[INFO] checking {len(files)} metabolic-resilience source files")

    for path in files:
        text = path.read_text(encoding="utf-8", errors="replace")
        print(f"[CHECK] {path.relative_to(ROOT)}")

        if path.suffix == ".sh":
            failures.extend(check_bash_syntax(path))
            failures.extend(check_embedded_python(path, text))
        elif path.suffix == ".py":
            failures.extend(check_python(path))
        elif path.suffix == ".R":
            failures.extend(check_r(path))

        failures.extend(policy_checks(path, text))

    failures.extend(scientific_contract_checks())

    if failures:
        print("\n[FAIL] metabolic-resilience CI findings:", file=sys.stderr)
        for item in failures:
            print(f"- {item}", file=sys.stderr)
        return 1

    print("\n[PASS] syntax, policy, and scientific contracts passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
