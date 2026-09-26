#!/usr/bin/env python3
"""Static CI checks for metabolic-resilience server scripts."""

from __future__ import annotations

import ast
import pathlib
import re
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
HARDCODED_ROOT_RE = re.compile(r"(?m)^\s*ROOT=['\"]/srv/is-analysis['\"]\s*$")
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
    return [f"{path}: bash -n failed:\n{proc.stderr.strip()}"]


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


def check_python_heredocs(path: pathlib.Path, text: str) -> list[str]:
    errors: list[str] = []

    for start_line, code in extract_python_heredocs(text):
        try:
            ast.parse(code, filename=f"{path}:heredoc@{start_line}")
        except SyntaxError as exc:
            shell_line = start_line + (exc.lineno or 1) - 1
            errors.append(
                f"{path}: embedded Python syntax error near shell line "
                f"{shell_line}: {exc.msg}"
            )

    return errors


def policy_checks(path: pathlib.Path, text: str) -> list[str]:
    errors: list[str] = []

    if ERREXIT_RE.search(text):
        errors.append(f"{path}: prohibited errexit/set -e usage")

    if SECRET_RE.search(text):
        errors.append(f"{path}: possible hard-coded credential literal")

    if HARDCODED_ROOT_RE.search(text):
        errors.append(
            f"{path}: hard-coded ROOT=/srv/is-analysis; "
            "use a configurable IS_ANALYSIS_ROOT with /srv/is-analysis as the default"
        )

    return errors


def main() -> int:
    scripts = sorted(SCRIPT_DIR.glob("*.sh"))

    if not scripts:
        print(f"[INFO] no metabolic-resilience shell scripts yet under {SCRIPT_DIR}")
        return 0

    failures: list[str] = []

    print(f"[INFO] checking {len(scripts)} scripts")

    for path in scripts:
        text = path.read_text(encoding="utf-8")
        print(f"[CHECK] {path.relative_to(ROOT)}")

        failures.extend(check_bash_syntax(path))
        failures.extend(check_python_heredocs(path, text))
        failures.extend(policy_checks(path, text))

    if failures:
        print("\n[FAIL] metabolic-resilience CI findings:", file=sys.stderr)
        for item in failures:
            print(f"- {item}", file=sys.stderr)
        return 1

    print("\n[PASS] all metabolic-resilience static checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
