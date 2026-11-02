#!/usr/bin/env python3
"""Validate actual *previously generated* P5 gene-level edgeR artifact.

No source raw 10x download and no donor-level data. Checks only the frozen
aggregate full-tested gene table, preserving a source checksum. This is NOT
reproduction of the original GSE313736 authors' unavailable 80-gene set.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import tempfile
from pathlib import Path

EXPECTED = {
    "tested": 16050,
    "fdr05": 13717,
    "high_total": 83,
    "high_up": 67,
    "high_down": 16,
}
CONTROL_SIGNS = {"MYOC": 1, "FKBP5": 1, "FAM107A": 1, "ANGPTL7": 1}


def source_hash(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def inspect(path, *, enforce_expected):
    ids = set()
    n = fdr_n = high_up = high_down = 0
    controls = {}
    symbols = {}
    with gzip.open(path, "rt", encoding="utf8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        wanted = {"gene_id", "gene_symbol", "logFC", "FDR"}
        if not wanted.issubset(reader.fieldnames or ()):
            raise ValueError("REQUIRED_EDGER_COLUMNS_MISSING")
        for row in reader:
            gene = (row["gene_id"] or "").split(".")[0]
            if not gene.startswith("ENSG") or not gene[4:].isdigit():
                raise ValueError("BAD_GENE_ID")
            if gene in ids:
                raise ValueError("DUPLICATE_ENSEMBL_ID")
            ids.add(gene)
            symbol = row["gene_symbol"]
            fc, q = float(row["logFC"]), float(row["FDR"])
            if not math.isfinite(fc) or not math.isfinite(q) or not 0 <= q <= 1:
                raise ValueError("BAD_FINITE_LOGFC_OR_FDR")
            n += 1
            fdr_n += q < 0.05
            high = q < 0.05 and abs(fc) >= 3
            high_up += bool(high and fc > 0)
            high_down += bool(high and fc < 0)
            if symbol in CONTROL_SIGNS:
                if symbol in controls:
                    raise ValueError("DUPLICATED_CONTROL_GENE_SYMBOL")
                controls[symbol] = {"gene_id": gene, "log2FC": fc, "FDR": q}
            symbols[symbol] = symbols.get(symbol, 0) + 1
    if n == 0:
        raise ValueError("EMPTY_GENE_RESULT")
    summary = {
        "tested_genes": n, "fdr_lt_0_05": fdr_n,
        "high_effect_rule": "FDR < 0.05 AND abs(log2FC) >= 3.0",
        "high_effect_total": high_up + high_down,
        "high_effect_up": high_up, "high_effect_down": high_down,
        "controls": controls,
        "p5_source_gene_table_sha256": source_hash(path),
        "author_original_80_gene_set_reproduced": False,
        "author_original_80_status": "HOLD_NO_AUTHENTICATED_AUTHOR_GENE_TABLE",
        "not_independent_biological_replication": True,
        "drug_vs_library_preparation_unresolved": True,
    }
    if enforce_expected:
        actual = {
            "tested": n, "fdr05": fdr_n, "high_total": high_up + high_down,
            "high_up": high_up, "high_down": high_down,
        }
        if actual != EXPECTED:
            raise ValueError("FROZEN_P5_STATISTIC_DRIFT " + json.dumps(actual))
        if any(g not in controls or controls[g]["log2FC"] * sign <= 0
               or controls[g]["FDR"] >= 0.05 for g, sign in CONTROL_SIGNS.items()):
            raise ValueError("P5_POSITIVE_CONTROL_DRIFT")
    return summary


def fixture_self_test():
    with tempfile.TemporaryDirectory() as work:
        p = Path(work) / "fixture.tsv.gz"
        rows = [
            ("ENSG00000000001", "MYOC", 4.0, 0.01),
            ("ENSG00000000002", "FKBP5", 4.0, 0.01),
            ("ENSG00000000003", "NEG", -3.1, 0.01),
            ("ENSG00000000004", "ZERO", 0.0, 0.1),
        ]
        def write(items):
            with gzip.open(p, "wt", encoding="utf8", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["gene_id", "gene_symbol", "logFC", "FDR"])
                w.writerows(items)
        write(rows)
        q = inspect(p, enforce_expected=False)
        assert (q["tested_genes"], q["fdr_lt_0_05"],
                q["high_effect_up"], q["high_effect_down"]) == (4, 3, 2, 1)
        write(rows + [rows[0]])
        try:
            inspect(p, enforce_expected=False)
            raise AssertionError("Duplicate not rejected")
        except ValueError as e:
            assert "DUPLICATE" in str(e)
        write(rows[:3] + [("ENSG00000000004", "ZERO", 0.0, 1.7)])
        try:
            inspect(p, enforce_expected=False)
            raise AssertionError("Invalid FDR not rejected")
        except ValueError as e:
            assert "BAD_FINITE" in str(e)
        write(rows)
        try:
            inspect(p, enforce_expected=True)
            raise AssertionError("Wrong frozen counts not rejected")
        except ValueError as e:
            assert "FROZEN_P5_STATISTIC_DRIFT" in str(e)
    print("P5_REFERENCE_FREEZE_SELF_TEST_PASS", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-dir", type=Path)
    p.add_argument("--out", type=Path)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        fixture_self_test()
        return
    if not a.artifact_dir or not a.out:
        p.error("--artifact-dir and --out required outside self-test")
    allfiles = list(a.artifact_dir.rglob("p5_all_gene_results.tsv.gz"))
    if len(allfiles) != 1:
        raise SystemExit("EXPECTED_ONE_P5_AGGREGATE_FILE " + str(len(allfiles)))
    result = inspect(allfiles[0], enforce_expected=True)
    summaries = list(a.artifact_dir.rglob("p5_edger_summary.json"))
    if len(summaries) != 1:
        raise SystemExit("EXPECTED_ONE_P5_SUMMARY_JSON " + str(len(summaries)))
    q = json.loads(summaries[0].read_text(encoding="utf8"))
    if (q.get("biological_pairs") != 57
        or q.get("tested_genes") != result["tested_genes"]
        or q.get("fdr_lt_0_05") != result["fdr_lt_0_05"]):
        raise SystemExit("P5_SUMMARY_VERSUS_GENE_TABLE_INCONSISTENT")
    result.update({
        "status": "PASS_OUR_57DONOR_AGGREGATE_REAUDIT_ONLY",
        "source_github_actions_run_id": 41074586446,
        "source_artifact_id": 12987425927,
        "source_artifact_name": "ocular-p5-edger-inference",
        "biological_donor_pairs": 57,
        "original_author_80_gene_ids_available": False,
    })
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "p5_57donor_83high_effect_reference_freeze.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    print("P5_REFERENCE_FREEZE_PASS", "tested", result["tested_genes"],
          "fdr", result["fdr_lt_0_05"],
          "high", result["high_effect_total"], "up", result["high_effect_up"],
          "down", result["high_effect_down"], flush=True)


if __name__ == "__main__":
    main()
