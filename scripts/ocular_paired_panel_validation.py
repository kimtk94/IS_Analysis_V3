#!/usr/bin/env python3
"""Evaluate frozen ocular panels in the donor-paired SRP319815 STAR-count release."""
import argparse
import csv
import hashlib
import io
import json
import math
import re
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from ocular_frozen_241_external import PIN, read_panel, sha

BUNDLE_SHA = "f197a59505b87c08f2df8552cc1732e2cadd1a6c529975b35c80d9ae5173d421"
ENA_SHA = "c60a84d47f003e5050f466dd05eb80131279192a11ff9a3efd03c4f1657ac20a"
ENA_FIELDS = ("run_accession", "secondary_study_accession", "study_accession",
              "experiment_accession", "sample_accession", "sample_alias", "library_name", "experiment_title")
PANELS = ("all_241", "ECM_16")


def ena_design(path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = sorted(csv.DictReader(handle, delimiter="\t"), key=lambda row: row["run_accession"])
    payload = "".join("\t".join(row[key] for key in ENA_FIELDS) + "\n" for row in rows)
    if hashlib.sha256(payload.encode()).hexdigest() != ENA_SHA:
        raise ValueError("ENA study or pairing evidence differs from reviewed metadata")
    design = defaultdict(dict)
    for row in rows:
        donor, arm = re.fullmatch(r"(\d{2}-\d{2})(DEX|ETH)", row["sample_alias"]).groups()
        design[donor][arm] = row
    if len(design) != 8 or any(set(arms) != {"DEX", "ETH"} for arms in design.values()):
        raise ValueError("Expected eight ENA pairs")
    return design


def gene_counts(handle, panel, runs):
    reader = csv.reader(handle, delimiter="\t")
    label, *columns = next(reader)
    if label != "GeneID" or len(set(columns)) != len(columns) or set(columns) != set(runs):
        raise ValueError("Count matrix and run metadata disagree")
    totals = [0] * len(columns)
    genes = set()
    selected = {}
    for gene, *fields in reader:
        if gene in genes or not re.fullmatch(r"ENSG\d+", gene) or len(fields) != len(columns):
            raise ValueError("Malformed or nonunique Ensembl count row")
        values = list(map(int, fields))
        if min(values) < 0:
            raise ValueError("Negative gene count")
        genes.add(gene)
        totals = [total + count for total, count in zip(totals, values)]
        if gene in panel:
            selected[gene] = values
    if not all(totals) or not genes:
        raise ValueError("Empty expression library")
    log_cpm = {
        gene: {run: math.log2(1 + value * 1_000_000 / total)
               for run, value, total in zip(columns, values, totals)}
        for gene, values in selected.items()
    }
    return log_cpm, len(genes)


def exact_sign(scores):
    positive = sum(x > 0 for x in scores)
    negative = sum(x < 0 for x in scores)
    n = positive + negative
    if not scores:
        p_value = None
    elif n == 0:
        p_value = 1.0
    else:
        tail = sum(math.comb(n, k) for k in range(min(positive, negative) + 1))
        p_value = min(1.0, 2 * tail / 2**n)
    return {"pairs": len(scores), "positive": positive, "negative": negative,
            "ties": len(scores) - n, "median_signed_log2_difference": statistics.median(scores) if scores else None,
            "two_sided_p": p_value}


def summarize_scores(scores):
    results = {key: exact_sign(scores[key]) for key in PANELS}
    running = 0.0
    ordered = sorted(PANELS, key=lambda k: (results[k]["two_sided_p"] is None, results[k]["two_sided_p"]))
    for i, key in enumerate(ordered):
        p = results[key]["two_sided_p"]
        if p is None:
            results[key]["holm_p"] = None
        else:
            running = max(running, min(1.0, (len(PANELS) - i) * p))
            results[key]["holm_p"] = running
    return results


def evaluate(panel, design, runs, expression):
    groups = {"all_241": list(expression), "ECM_16": [g for g in expression if panel[g]["ecm_member"]]}
    if any(not genes for genes in groups.values()):
        raise ValueError("No mapped genes for a required panel")
    scores = {name: defaultdict(list) for name in ("primary", "all_usable")}
    pair_coverage = []
    for donor, arms in sorted(design.items()):
        dex, vehicle = (arms[k]["run_accession"] for k in ("DEX", "ETH"))
        qc = {arm: runs[row["run_accession"]]["QC_summary"] if row["run_accession"] in runs
              else "not_in_count_release" for arm, row in arms.items()}
        complete = dex in runs and vehicle in runs
        primary = complete and all(flag in {"PASS", "WARN"} for flag in qc.values())
        pair_coverage.append({"cell_line": donor, "source_group": arms["DEX"]["experiment_title"],
                              "arm_qc": qc, "complete": complete, "in_primary": primary})
        if not complete:
            continue
        for key, genes in groups.items():
            score = statistics.mean(
                (expression[g][dex] - expression[g][vehicle]) * (1 if panel[g]["acute_log2fc"] > 0 else -1)
                for g in genes)
            scores["all_usable"][key].append(score)
            if primary:
                scores["primary"][key].append(score)
    return {"mapped": {key: len(genes) for key, genes in groups.items()},
            "pair_coverage": pair_coverage,
            "primary_group_counts": dict(Counter(p["source_group"] for p in pair_coverage if p["in_primary"])),
            **{name: summarize_scores(values) for name, values in scores.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("panel", "ena", "bundle", "out"):
        parser.add_argument("--" + key, required=True, type=Path)
    args = parser.parse_args()
    panel = read_panel(args.panel)
    design = ena_design(args.ena)
    if sha(args.bundle) != BUNDLE_SHA:
        raise ValueError("DEE2 project bundle checksum mismatch")
    ena_runs = {r["run_accession"]: r for arms in design.values() for r in arms.values()}
    with zipfile.ZipFile(args.bundle) as z:
        with z.open("SRP319815_NA/MetadataSummary.tsv") as handle:
            metadata = list(csv.DictReader(io.TextIOWrapper(handle), delimiter="\t"))
        runs = {r["SRR_accession"]: r for r in metadata}
        if len(runs) != len(metadata):
            raise ValueError("Duplicate run in DEE2 metadata")
        for run, r in runs.items():
            if (run not in ena_runs or r["SRP_accession"] != "SRP319815" or
                    r["SRX_accession"] != ena_runs[run]["experiment_accession"] or
                    r["QC_summary"] not in {"PASS", "WARN", "FAIL"}):
                raise ValueError("Unexpected experiment, project or QC for " + run)
        with z.open("SRP319815_NA/GeneCountMatrix.tsv") as handle:
            expression, num_genes = gene_counts(io.TextIOWrapper(handle), panel, runs)
    result = {"study": "PRJNA729873", "bundle_sha256": BUNDLE_SHA, "ena_projection_sha256": ENA_SHA,
              "panel_sha256": PIN, "source_runs": len(ena_runs), "counted_runs": len(runs),
              "genome_wide_genes": num_genes, "frozen": {"all_241": 241, "ECM_16": 16},
              **evaluate(panel, design, runs, expression)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result["primary"], indent=2))


if __name__ == "__main__":
    main()
