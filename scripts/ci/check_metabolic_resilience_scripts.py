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

    errors += require(
        v2 / "62_stage4c_select_eur_candidates.py",
        r"MR_FDR\s*=\s*0\.05.*COLOC_H4\s*=\s*0\.80",
        "EUR candidate gate requires MR FDR 0.05 and coloc H4 0.80",
    )
    errors += require(
        v2 / "62_stage4c_select_eur_candidates.py",
        r"len\(domains\)\s*>=\s*2",
        "multi-domain EUR priority requires at least 2 metabolic domains",
    )
    errors += require(
        v2 / "64_stage4e_cross_ancestry_summarize.py",
        r"CONCORDANT_SIGNIFICANT.*CONCORDANT_NONSIGNIFICANT.*OPPOSITE_DIRECTION.*NOT_TESTABLE",
        "EAS replication keeps four interpretation classes",
    )
    errors += require(
        v2 / "71_stage5b_koges_feasibility_gate.py",
        r"n_gen_long\s*>=\s*4000",
        "KoGES genotype-longitudinal GO gate requires N >= 4000",
    )
    errors += require(
        v2 / "71_stage5b_koges_feasibility_gate.py",
        r"physical_activity_usable_waves.*<\s*3.*USE_BASELINE_PA_STRATIFICATION",
        "KoGES PA fallback uses baseline stratification when repeated PA <3 waves",
    )
    errors += require(
        v2 / "73_stage5d_write_analysis_contract.py",
        r"incident MetS.*Cox proportional hazards.*Repeated continuous Metabolic Burden Index.*linear mixed model",
        "Stage5 contract locks incident MetS primary and repeated MBI secondary",
    )
    errors += require(
        v2 / "46_stage3d6_multiple_testing_overlap.py",
        r"T2D_VALIDATION.*BH-FDR within prespecified metabolic domain",
        "T2D is a separate validation family and metabolic multiplicity is prespecified",
    )
    errors += require(
        v2 / "47_stage3d7_steiger_sensitivity.py",
        r"HOLD_BINARY_LIABILITY_ASSUMPTIONS",
        "binary-outcome Steiger remains HOLD without liability-scale assumptions",
    )
    errors += require(
        v2 / "54_stage3e4_multisignal_coloc_gate.py",
        r"RUN_MULTISIGNAL_COLOC.*OPTIONAL_MULTISIGNAL_SENSITIVITY.*ABF_SINGLE_SIGNAL_ACCEPTABLE",
        "multi-signal loci are explicitly routed to conditional/SuSiE sensitivity",
    )
    errors += require(
        v2 / "55_stage3f0_variant_artifact_audit.py",
        r"EXCLUDE_AND_RERUN_MR_COLOC",
        "protein-altering/coding/splice variants trigger exclusion sensitivity",
    )
    errors += require(
        v2 / "65_stage4f_eas_ld_reference_plan.py",
        r"eur_ld_allowed_primary.*0",
        "EUR LD is prohibited as the primary EAS LD reference",
    )
    errors += require(
        v2 / "74_stage5e_koges_longitudinal_models.R",
        r"coxph.*LMM_MBI_SLOPE",
        "KoGES primary Cox and secondary longitudinal MBI models are implemented",
    )
    errors += require(
        v2 / "75_stage5f_koges_pa_interaction.R",
        r"physical_activity.*time_years",
        "KoGES physical-activity interaction model includes longitudinal interaction",
    )
    errors += require(
        v2 / "80_stage6a_functional_annotation_scaffold.py",
        r"must not override MR/coloc evidence",
        "functional annotation remains supportive rather than causal evidence",
    )
    errors += require(
        v2 / "91_master_analysis_gate.py",
        r"(?:PASS_CORE_READY.*HOLD_CORE_INCOMPLETE|HOLD_CORE_INCOMPLETE.*PASS_CORE_READY)",
        "master manuscript gate blocks incomplete core causal claims",
    )
    errors += require(
        v2 / "92_environment_snapshot.sh",
        r"pip freeze.*sessionInfo",
        "reproducibility snapshot captures Python and R environments",
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
