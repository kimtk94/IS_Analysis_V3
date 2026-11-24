#!/usr/bin/env python3
"""Bounded public-data follow-up: descriptive ocular rank shifts, CKD source coverage, and missing-run inventory."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
import zipfile
from pathlib import Path

GENES = ("ACP1", "CPVL", "F12", "GSTA1", "GSTA3", "HLA-E", "INHBC", "SDCCAG8", "UMOD")
# Frozen accession mapping from the GEO GSE16643 series record.
SAMPLES = {
    729: ("93", "TA", "C"), 730: ("93", "DMSO", "C"),
    731: ("86", "TA", "B"), 732: ("93", "DMSO", "B"),
    733: ("93", "TA", "A"), 734: ("86", "DMSO", "B"),
    735: ("86", "FA", "B"), 736: ("86", "TA", "C"),
    737: ("86", "DMSO", "C"), 738: ("86", "DMSO", "A"),
    739: ("86", "TA", "A"), 740: ("93", "DMSO", "A"),
    741: ("93", "Dex", "C"), 742: ("86", "FA", "C"),
    743: ("86", "FA", "A"), 744: ("93", "Dex", "A"),
    745: ("93", "Dex", "B"), 746: ("93", "FA", "C"),
    747: ("86", "Dex", "B"), 748: ("93", "FA", "A"),
    749: ("93", "FA", "B"), 750: ("86", "Dex", "C"),
    751: ("86", "Dex", "A"),
}
SAMPLES = {"GSM417" + str(k): v for k, v in SAMPLES.items()}


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def read_table(path, begin, end):
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="strict") as f:
        lines = []
        for line in f:
            if line.strip().lower() == begin:
                break
            lines.append(line)
        else:
            raise ValueError("Missing table start")
        body = []
        for line in f:
            if line.strip().lower() == end:
                break
            body.append(line)
        else:
            raise ValueError("Missing table end")
    return lines, csv.reader(body, delimiter="\t")


def percentile(values):
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i + 1
        while j < n and values[order[i]] == values[order[j]]:
            j += 1
        r = (i + j) / (2 * n)
        for k in range(i, j):
            ranks[order[k]] = r
        i = j
    return ranks


def ocular(root):
    from ocular_frozen_241_external import read_panel
    panel = read_panel(root / "panel.zip")
    symbols = collections.defaultdict(list)
    for gid, r in panel.items():
        symbols[r["symbol"].strip().upper()].append(gid)
    allowed = {s: ids[0] for s, ids in symbols.items() if s and len(ids) == 1}
    _, table = read_table(root / "platform.gz", "!platform_table_begin", "!platform_table_end")
    h = next(table)
    required = ("ID", "CONTROL_TYPE", "GENE_SYMBOL")
    if not set(required) <= set(h) or len(set(h)) != len(h):
        raise ValueError("Unexpected GPL6480 annotation schema")
    probe_gene = {}
    for row in table:
        if len(row) != len(h):
            raise ValueError("Malformed platform row")
        d = dict(zip(h, row))
        if d["CONTROL_TYPE"].strip().upper() == "FALSE":
            p = d["ID"]
            if p in probe_gene:
                raise ValueError("Duplicate platform ID")
            probe_gene[p] = allowed.get(d["GENE_SYMBOL"].strip().upper())
    meta, table = read_table(root / "matrix.gz", "!series_matrix_table_begin", "!series_matrix_table_end")
    h = next(table)
    ids = h[1:]
    if h[0] != "ID_REF" or set(ids) != set(SAMPLES) or len(ids) != 23:
        raise ValueError("Unexpected GSE16643 matrix sample identity")
    titles = [next(csv.reader([l], delimiter="\t"))[1:] for l in meta if l.startswith("!Sample_title\t")]
    geoids = [next(csv.reader([l], delimiter="\t"))[1:] for l in meta if l.startswith("!Sample_geo_accession\t")]
    if len(titles) != 1 or len(geoids) != 1 or geoids[0] != ids or len(titles[0]) != len(ids):
        raise ValueError("Missing source sample labels")
    for accession, title in zip(ids, titles[0]):
        strain, arm, rep = SAMPLES[accession]
        expected = f"TM {strain} treated with " + ("0.1% DMSO" if arm == "DMSO" else f"{arm} at 1uM") + f", rep {rep}"
        if title.strip().casefold() != expected.casefold():
            raise ValueError("Source title differs from fixed design: " + accession)
    rows, pids, seen = [], [], set()
    missing = 0
    for row in table:
        if len(row) != len(h) or row[0] in seen:
            raise ValueError("Malformed or duplicate matrix row")
        seen.add(row[0])
        if row[0] not in probe_gene:
            continue
        try:
            a = list(map(float, row[1:]))
        except ValueError:
            missing += 1
            continue
        if not all(math.isfinite(v) for v in a):
            missing += 1
            continue
        rows.append(a)
        pids.append(row[0])
    if len(rows) < 10000:
        raise ValueError("Too few complete non-control probes")
    gene_probes = collections.defaultdict(list)
    for i, p in enumerate(pids):
        if probe_gene[p]:
            gene_probes[probe_gene[p]].append(i)
    rank_by_sample = {}
    for accession, v in zip(ids, zip(*rows)):
        ranks = percentile(v)
        rank_by_sample[accession] = {g: statistics.median(ranks[i] for i in ix) for g, ix in gene_probes.items()}
    out = {"status": "DESCRIPTIVE_ONLY", "study": "GSE16643", "platform": "GPL6480",
           "matrix_sha256": digest(root / "matrix.gz"), "platform_sha256": digest(root / "platform.gz"),
           "panel_sha256": digest(root / "panel.zip"), "complete_probes": len(rows), "excluded_missing_probes": missing,
           "biological_strains": 2, "samples": len(ids), "panels": {}}
    for key, genes in [("all241", list(panel)), ("ECM16", [g for g,r in panel.items() if r["ecm_member"]])]:
        mapped = [g for g in genes if g in gene_probes]
        if not mapped:
            raise ValueError("No exact mapped candidate symbols")
        scores = {}
        for strain in ("86", "93"):
            for arm in ("Dex", "FA", "TA"):
                a = [i for i in ids if SAMPLES[i][:2] == (strain, arm)]
                b = [i for i in ids if SAMPLES[i][:2] == (strain, "DMSO")]
                if not a or not b:
                    raise ValueError("Missing arm")
                delta = [ (statistics.mean(rank_by_sample[i][g] for i in a) -
                           statistics.mean(rank_by_sample[i][g] for i in b)) *
                           (1 if panel[g]["acute_log2fc"] > 0 else -1) for g in mapped ]
                scores[strain+"_"+arm] = {"treatment_arrays": len(a), "vehicle_arrays": len(b),
                    "mean_signed_percentile_shift": statistics.mean(delta),
                    "positive_gene_directions": sum(x > 0 for x in delta), "ties": sum(x == 0 for x in delta)}
        out["panels"][key] = {"frozen": len(genes), "mapped": len(mapped), "scores": scores}
    return out


def norm(x):
    return re.sub(r"[^a-z0-9]", "", str(x or "").lower())


def ckd(root):
    from openpyxl import load_workbook
    path = root / "ckb.xlsx"
    w = load_workbook(path, read_only=True, data_only=True)
    genes = set(GENES)
    inventory = []
    recognized = 0
    # Require an explicit gene-symbol field. No substring matching; do not infer unreported candidates.
    gene_names = {"gene", "genesymbol", "hgncsymbol", "hgnc", "genename", "proteinencodinggene", "targetgene"}
    for s in w.worksheets:
        headers = []
        chosen = None
        counts = collections.Counter()
        rowcount = 0
        preview = []
        for i, row in enumerate(s.iter_rows(values_only=True), 1):
            if i <= 12:
                nonempty = [str(x)[:100] for x in row if x is not None]
                preview.append(nonempty[:30])
            if chosen is None and i <= 20:
                ix = [j for j, value in enumerate(row) if norm(value) in gene_names]
                if len(ix) == 1:
                    chosen = ix[0]
                    headers = [str(x or "")[:160] for x in row]
                    continue
            if chosen is not None:
                rowcount += 1
                value = str(row[chosen] or "").strip().upper()
                if value in genes:
                    counts[value] += 1
        if chosen is not None:
            recognized += 1
        inventory.append({"sheet": s.title, "header_recognized": chosen is not None,
                          "header": headers, "first_rows": preview if chosen is None else [],
                          "rows_scanned_after_header": rowcount, "exact_candidate_row_counts": dict(counts)})
    w.close()
    return {"status": "SOURCE_TABLE_REVIEW_REQUIRED", "source": "10.1038/s41467-025-56935-2 Supplementary Data 4",
            "sha256": digest(path), "candidate_genes": list(GENES), "recognized_sheets": recognized,
            "sheets": inventory, "limitation": "Gene-row presence is not validated cis instrument, platform equivalence or MR."}


def recovery(root):
    from ocular_paired_panel_validation import ena_design
    from ocular_paired_panel_validation import BUNDLE_SHA
    design = ena_design(root / "ena.tsv")
    p = root / "counts.zip"
    if digest(p) != BUNDLE_SHA:
        raise ValueError("Previously pinned DEE2 bundle differs")
    with zipfile.ZipFile(p) as z:
        import io
        with z.open("SRP319815_NA/MetadataSummary.tsv") as f:
            rows = list(csv.DictReader(io.TextIOWrapper(f), delimiter="\t"))
    available = {r["SRR_accession"]: r for r in rows}
    if len(available) != len(rows):
        raise ValueError("Duplicate DEE2 run")
    with (root / "fastq.tsv").open() as f:
        raw = list(csv.DictReader(f, delimiter="\t"))
    entries = {r["run_accession"]: r for r in raw}
    all_runs = {r["run_accession"] for a in design.values() for r in a.values()}
    if set(entries) != all_runs or len(entries) != len(raw) or not set(available) <= all_runs:
        raise ValueError("ENA run manifests disagree")
    missing = []
    complete = 0
    source_bytes = 0
    for a in design.values():
        if all(r["run_accession"] in available for r in a.values()):
            complete += 1
        for arm, record in a.items():
            run = record["run_accession"]
            if run in available and (available[run]["SRX_accession"] != record["experiment_accession"] or
                                     available[run]["SRP_accession"] != "SRP319815"):
                raise ValueError("DEE2 project or experiment mismatch")
            x = entries[run]
            files = x["fastq_ftp"].split(";")
            md5 = x["fastq_md5"].split(";")
            sizes = list(map(int, x["fastq_bytes"].split(";")))
            if not len(files) == len(md5) == len(sizes) == 2 or not all(re.fullmatch("[a-fA-F0-9]{32}", m) for m in md5) or not all(n > 0 for n in sizes):
                raise ValueError("Incomplete paired-end source manifest")
            source_bytes += sum(sizes)
            if run not in available:
                missing.append({"run": run, "arm": arm, "source_group": record["experiment_title"],
                                "paired_fastq_bytes": sizes, "declared_md5": md5})
    return {"status": "HOLD_RAW_REQUANTIFICATION", "source_runs": len(all_runs), "complete_dee2_pairs": complete,
            "missing_runs": missing, "all_raw_fastq_bytes": source_bytes,
            "fastq_manifest_sha256": digest(root / "fastq.tsv"),
            "limitation": "Declared MD5 is not verification of FASTQ bytes. Requantify all 16 with one reference and pipeline."}


def selftest():
    assert percentile([4, 4, 1, 9]) == [0.5, 0.5, 0.125, 0.875]
    assert len(SAMPLES) == 23 and len(set(SAMPLES.values())) == 23
    assert norm("Gene Symbol") == "genesymbol"
    assert len(GENES) == 9
    print("Four small synthetic contract checks passed; no biological result.")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["ocular", "ckd", "recovery", "selftest"])
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--out", type=Path, default=Path("result.json"))
    a = p.parse_args()
    if a.mode == "selftest":
        selftest()
        return
    out = {"ocular": ocular, "ckd": ckd, "recovery": recovery}[a.mode](a.root)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    print(json.dumps(out, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
