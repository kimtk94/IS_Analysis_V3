#!/usr/bin/env python3
"""Run self-contained regressions or explicitly selected historical source checks.

The four historical suites authenticate an immutable, separately checked-out
source snapshot. They are integration checks, not prerequisites for testing a
fresh source checkout. Neither mode downloads data or changes source hash pins.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
LEGACY_MODULES = frozenset({
    "test_is_24_original_lead_review.py",
    "test_is_gws_eas_24_original_leads.py",
    "test_is_p3a_36_group_recovery.py",
    "test_is_p3a_chr10_legacy_original_signature.py",
})
LEGACY_INPUTS = (
    "IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv",
    "IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv",
    "IS_EAS5_98_BH05_FDR_EGENE_FULLCIS_GWAS_ALLELE_READINESS.tsv",
    "IS_EAS5_SAME_SNP_ASSOCIATION_GENE_TISSUE_QC.tsv",
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("fixtures", "legacy", "all"),
                        default="fixtures")
    args = parser.parse_args(argv)
    os.chdir(REPO_ROOT)

    if args.suite in {"legacy", "all"}:
        root = REPO_ROOT / "legacy" / "docs" / "data" / "is"
        missing = [str(root / name) for name in LEGACY_INPUTS
                   if not (root / name).is_file()]
        if missing:
            print("Historical source validation was requested, but inputs are missing:",
                  file=sys.stderr)
            print("\n".join(missing), file=sys.stderr)
            print("Prepare the pinned legacy checkout documented in README.md. "
                  "Use --suite fixtures for self-contained code regressions.",
                  file=sys.stderr)
            return 2

    test_dir = REPO_ROOT / "tests"
    sys.path.insert(0, str(REPO_ROOT))
    sys.path.insert(0, str(test_dir))
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    selected = []
    for path in sorted(test_dir.glob("test_*.py")):
        historical = path.name in LEGACY_MODULES
        if args.suite == "fixtures" and historical:
            continue
        if args.suite == "legacy" and not historical:
            continue
        selected.append(path.name)
        suite.addTests(loader.discover(str(test_dir), pattern=path.name))

    print(f"Repository test suite: {args.suite}; "
          f"modules={len(selected)}; cases={suite.countTestCases()}", flush=True)
    if args.suite == "fixtures":
        print("Historical source integration checks are separate: "
              "bash scripts/codex_smoke_test.sh --with-legacy", flush=True)
    if not suite.countTestCases():
        print("No tests selected; refusing to report success.", file=sys.stderr)
        return 2
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
