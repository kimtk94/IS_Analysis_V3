#!/usr/bin/env python3
"""Validate real GSE331114 10x matrices and derive descriptive donor-level counts.

GEO GSE331114: 4 acute ischemic stroke (AIS), 4 non-stroke brain samples.
CD45-enriched brain immune cells, NOT unbiased whole brain or an eQTL atlas.

Downloads official per-GSM processed MEX gzip files on ephemeral GitHub Actions
disk. No patient-level sequence data, cell barcodes, raw matrices or identities
are uploaded as artifacts. The reported sample-level gene panel is exploratory.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen

SOFT = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE331nnn/"
        "GSE331114/soft/GSE331114_family.soft.gz")
SAMPLES = {
    **{f"GSM973925{i}": ("AIS", i) for i in range(1, 5)},
    **{f"GSM973925{i}": ("Control", i - 4) for i in range(5, 9)},
}
CANDIDATES = (
    "ABO", "ALDH2", "CD40", "F11", "FGF5", "FURIN", "GRK5", "IL6RA",
    "KIAA0319", "MMP12", "PROCR", "SCARA5", "SWAP70", "TFPI", "TMPRSS5",
)
MARKERS = (
    "PTPRC", "CD3D", "CD8A", "CD8B", "NKG7", "LST1", "C1QA", "C1QB",
    "CD68", "MS4A7", "S100A8", "S100A9", "PECAM1", "AQP4", "RBFOX3",
    "MBP", "TYROBP", "TREM2",
)
TARGETS = set(CANDIDATES + MARKERS)
MAX_FILES = {
    "matrix": 180 * 1024 * 1024,
    "features": 3 * 1024 * 1024,
    "barcodes": 3 * 1024 * 1024,
    "soft": 18 * 1024 * 1024,
}
QC = {"min_genes": 200, "max_genes": 6500, "min_umi": 500, "max_mito_fraction": 0.20}


def safe_file_url(value, expected_gsm):
    url = value.strip().replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/", 1)
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != "ftp.ncbi.nlm.nih.gov":
        raise ValueError("Untrusted supplementary file host")
    if expected_gsm not in parsed.path or "/suppl/" not in parsed.path:
        raise ValueError("Unexpected GEO sample file path")
    return url


def download(url, dest, cap, *, source_gsm=None):
    if source_gsm:
        safe_file_url(url, source_gsm)
    elif url != SOFT:
        raise ValueError("Untrusted metadata source")
    request = Request(url, headers={"User-Agent": "IS-Analysis-V3-GSE331114-MEX-QC/1.0"})
    count = 0
    digest = hashlib.sha256()
    with urlopen(request, timeout=70) as resp:
        resolved = urlparse(resp.geturl())
        if resp.status != 200 or resolved.scheme != "https" or resolved.hostname != "ftp.ncbi.nlm.nih.gov":
            raise ValueError("Unexpected redirect or HTTP response")
        declared = resp.headers.get("Content-Length")
        if declared and declared.isdigit() and int(declared) > cap:
            raise ValueError("Source file exceeds cap")
        with open(dest, "wb") as out:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                count += len(chunk)
                if count > cap:
                    raise ValueError("File exceeds bounded download cap")
                out.write(chunk)
                digest.update(chunk)
        if declared and declared.isdigit() and count != int(declared):
            raise IOError("Incomplete source transfer")
    return {"bytes": count, "sha256_local_not_publisher_verified": digest.hexdigest()}


def parse_soft_samples(path):
    samples = {}
    current = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                current = line.split("=", 1)[1].strip()
                samples[current] = {"title": "", "files": [], "organism": None}
            elif line.startswith("^"):
                current = None
            elif current in samples:
                if line.startswith("!Sample_title = "):
                    samples[current]["title"] = line.split("=", 1)[1].strip()
                elif line.startswith("!Sample_organism_ch1 = "):
                    samples[current]["organism"] = line.split("=", 1)[1].strip()
                elif line.startswith("!Sample_supplementary_file") and " = " in line:
                    samples[current]["files"].append(line.split(" = ", 1)[1].strip())
    return samples


def select_files(entry, gsm, group, number):
    matched = {}
    expected_label = ("AIS" if group == "AIS" else "Control") + "_brain" + str(number)
    for link in entry["files"]:
        filename = urlparse(link).path.rsplit("/", 1)[-1]
        if not filename.startswith(f"{gsm}_{expected_label}_"):
            continue
        for kind, suffix in {
            "matrix": "_matrix.mtx.gz",
            "features": "_features.tsv.gz",
            "barcodes": "_barcodes.tsv.gz",
        }.items():
            if filename.endswith(suffix):
                matched[kind] = safe_file_url(link, gsm)
    if set(matched) != {"matrix", "features", "barcodes"}:
        raise ValueError(f"{gsm} missing or mislabeled 10x MEX: {sorted(matched)}")
    if entry.get("organism") != "Homo sapiens":
        raise ValueError(f"{gsm} organism mismatch")
    return matched


def load_features(path):
    features, seen_ids = [], set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("Malformed features row")
            stable_id, symbol = fields[0], fields[1]
            if stable_id in seen_ids:
                raise ValueError("Duplicate feature IDs")
            seen_ids.add(stable_id)
            features.append((stable_id, symbol))
    if len(features) < 1000:
        raise ValueError("Unusually short feature table")
    return features


def count_barcodes(path):
    n = 0
    seen = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value or value in seen:
                raise ValueError("Duplicate or empty barcode")
            seen.add(value)
            n += 1
    if n < 100:
        raise ValueError("Unusually few barcodes")
    return n


def qc_matrix(path, features, ncells):
    genes = len(features)
    gene_umi = [0] * ncells
    gene_detected = [0] * ncells
    mito_umi = [0] * ncells
    targets = {symbol: [0] * ncells for symbol in TARGETS}
    target_rows = {i for i, (_, symbol) in enumerate(features) if symbol in TARGETS}
    mito_rows = {i for i, (_, symbol) in enumerate(features) if symbol.startswith("MT-")}
    rows = 0
    nonzeros = 0
    dim = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        first = handle.readline()
        if not first.startswith("%%MatrixMarket matrix coordinate integer"):
            raise ValueError("Non-integer coordinate MatrixMarket source")
        for line in handle:
            if line.startswith("%"):
                continue
            values = line.strip().split()
            if not values:
                continue
            if dim is None:
                if len(values) != 3:
                    raise ValueError("MatrixMarket missing dimensions")
                nfeatures, nbarcodes, nonzeros = map(int, values)
                if nfeatures != genes or nbarcodes != ncells or nonzeros < 100:
                    raise ValueError("MEX dimension mismatch")
                dim = (nfeatures, nbarcodes, nonzeros)
                continue
            if len(values) != 3:
                raise ValueError("Malformed nonzero MatrixMarket row")
            gi, ci, count = map(int, values)
            if not (1 <= gi <= genes and 1 <= ci <= ncells and count > 0):
                raise ValueError("Invalid coordinate/count")
            gi -= 1
            ci -= 1
            gene_umi[ci] += count
            gene_detected[ci] += 1
            if gi in mito_rows:
                mito_umi[ci] += count
            if gi in target_rows:
                targets[features[gi][1]][ci] += count
            rows += 1
        if dim is None or rows != dim[2]:
            raise ValueError(f"MatrixMarket nnz mismatch {rows}/{nonzeros}")
    passed = [
        i for i in range(ncells)
        if QC["min_genes"] <= gene_detected[i] <= QC["max_genes"]
        and gene_umi[i] >= QC["min_umi"]
        and mito_umi[i] / max(gene_umi[i], 1) < QC["max_mito_fraction"]
    ]
    sample_total = sum(gene_umi[i] for i in passed)
    report = {
        "source_matrix": {"features": genes, "barcodes": ncells, "nnz": nonzeros},
        "qc_thresholds_exploratory": QC,
        "barcodes_qc_pass": len(passed),
        "barcodes_qc_pass_fraction": round(len(passed) / ncells, 6),
        "median_umi_unfiltered": statistics.median(gene_umi),
        "median_detected_genes_unfiltered": statistics.median(gene_detected),
        "median_mito_pct_unfiltered": round(
            statistics.median(100 * mito_umi[i] / max(gene_umi[i], 1)
                              for i in range(ncells)), 3),
        "qc_pass_library_umi": sample_total,
        "no_doublet_removal": True,
        "cell_type_labels_in_raw_matrix": False,
    }
    expression = []
    for gene in sorted(TARGETS):
        matrix_gene = targets[gene]
        total_umi = sum(matrix_gene[i] for i in passed)
        n_detected = sum(matrix_gene[i] > 0 for i in passed)
        expression.append({
            "gene_symbol": gene, "is_prespecified_IS_candidate": gene in CANDIDATES,
            "total_umi_qc_barcodes": total_umi,
            "n_expressing_qc_barcodes": n_detected,
            "fraction_expressing_qc_barcodes": round(n_detected / max(len(passed), 1), 6),
            "sample_pseudobulk_cpm": round(
                1e6 * total_umi / sample_total, 6) if sample_total > 0 else None,
        })
    if len(passed) < 100:
        report["analysis_readiness"] = "HOLD_TOO_FEW_BARCODES"
    elif sample_total == 0:
        report["analysis_readiness"] = "HOLD_EMPTY_LIBRARY"
    else:
        report["analysis_readiness"] = "PSEUDOBULK_INPUT_QC_PASS_UNANNOTATED"
    return report, expression


def run_sample(gsm, out_dir):
    group, number = SAMPLES[gsm]
    status = {
        "sample": gsm,
        "group": group, "replicate_index": number,
        "cohort": "GSE331114 human brain parenchyma CD45-enriched",
        "evidence_type": "UNANNOTATED_CELL_BARCODE_MEX",
        "disease_response_not_prestroke_susceptibility": True,
        "donor_independence_not_verified_from_GSM": True,
        "status": "SOURCE_UNAVAILABLE",
        "count_matrix_qc": "NOT_VERIFIED",
        "celltype_annotation": "NOT_RUN",
        "genotype_eQTL": "NOT_AVAILABLE_FROM_GEO",
    }
    with tempfile.TemporaryDirectory(prefix="is_gse331114_") as temp:
        work = Path(temp)
        soft_path = work / "family.soft.gz"
        status["soft_file"] = download(SOFT, soft_path, MAX_FILES["soft"])
        records = parse_soft_samples(soft_path)
        if len(records) != 8 or set(records) != set(SAMPLES):
            raise ValueError(f"Expected 8 GSM entries, found {len(records)}")
        files = select_files(records[gsm], gsm, group, number)
        paths = {}
        for kind, link in files.items():
            target = work / (kind + (".mtx.gz" if kind == "matrix" else ".tsv.gz"))
            status[kind + "_source"] = {
                "url": link,
                **download(link, target, MAX_FILES[kind], source_gsm=gsm)
            }
            paths[kind] = target
        features = load_features(paths["features"])
        ncells = count_barcodes(paths["barcodes"])
        metrics, genes = qc_matrix(paths["matrix"], features, ncells)
        status.update(metrics)
        status["status"] = "MATRIX_QC_PASS"
        status["count_matrix_qc"] = "PASS"
        status["candidate_panel_detected_in_feature_symbols"] = sum(
            any(sym == g for _, sym in features) for g in CANDIDATES
        )
        status["report_time_utc"] = datetime.now(timezone.utc).isoformat()
        sample_out = out_dir / gsm
        sample_out.mkdir(parents=True, exist_ok=True)
        (sample_out / "qc.json").write_text(
            json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        with open(sample_out / "candidate_panel.tsv", "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(genes[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(genes)
        print(f"GSE331114 {gsm} {group}: MATRIX_QC_PASS "
              f"features={metrics['source_matrix']['features']} "
              f"barcodes={ncells} qc_barcodes={metrics['barcodes_qc_pass']} "
              f"candidate_panel_found={status['candidate_panel_detected_in_feature_symbols']}/15")
    return status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", choices=sorted(SAMPLES), required=True)
    parser.add_argument("--out-dir", default="artifacts/gse331114")
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    try:
        result = run_sample(args.sample, out)
    except (OSError, ValueError, EOFError, RuntimeError, TimeoutError) as exc:
        result = {"sample": args.sample, "group": SAMPLES[args.sample][0],
                  "status": "SOURCE_UNAVAILABLE_OR_QC_ERROR",
                  "error_type": type(exc).__name__, "error_message_short": str(exc)[:250],
                  "count_matrix_qc": "HOLD", "genotype_eQTL": "NOT_AVAILABLE_FROM_GEO"}
        sample_out = out / args.sample
        sample_out.mkdir(parents=True, exist_ok=True)
        (sample_out / "qc.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result))
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
