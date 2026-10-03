#!/usr/bin/env python3
"""Stage 1 audit for genetic studies of resistance-training trainability.

The script intentionally separates reported claims from independently recomputed
significance tiers. It can fetch the open-access 2026 article from Europe PMC,
extract tables containing the published lead SNPs, and write an auditable set of
TSV/JSON outputs. The 2024 study remains metadata-only until approved summary
statistics are supplied by the user.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional

EUROPE_PMC_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13371584/fullTextXML"
PHDA_URL = "https://www.ncmi.cn/"

SEED_2026_RSIDS = [
    "rs10212396", "rs12519717", "rs12055037", "rs2131183", "rs75968146",
    "rs4370982", "rs74038095", "rs73966436", "rs12625907",
]

STUDIES = [
    {
        "study_id": "yang_2024_physiol_genomics",
        "year": 2024,
        "pmid": "38881426",
        "doi": "10.1152/physiolgenomics.00019.2024",
        "population": "physically inactive adults; China multi-centre cohort",
        "n_total": 440,
        "phenotype": "12-week change in rectus femoris muscle thickness (MTRF)",
        "exercise": "RT or HIIT",
        "reported_leads": "11 RT + 8 HIIT",
        "reported_threshold": "1e-5",
        "data_access": "request_required",
        "data_url": PHDA_URL,
        "role": "discovery/supporting study; full summary statistics pending approval",
    },
    {
        "study_id": "gu_2026_jcsm",
        "year": 2026,
        "pmid": "42455518",
        "doi": "10.1002/jcsm.70347",
        "population": "young Asian sedentary adults",
        "n_total": 187,
        "phenotype": "12-week change in DXA lean body mass (delta LBM)",
        "exercise": "resistance training",
        "reported_leads": "9",
        "reported_threshold": "5e-8 in abstract; manuscript/table values audited independently",
        "data_access": "open_article_tables",
        "data_url": EUROPE_PMC_XML,
        "role": "primary public lead-locus audit",
    },
]

@dataclass
class VariantRow:
    rsid: str
    p_value: Optional[float] = None
    source: str = "seed"
    raw_row: str = ""
    reported_gene: str = ""
    chr: str = ""
    pos: str = ""
    effect_allele: str = ""
    other_allele: str = ""
    beta: str = ""
    tier: str = "unresolved"
    genome_wide_5e8: bool = False
    suggestive_1e5: bool = False

def classify_p(p: Optional[float]) -> tuple[str, bool, bool]:
    if p is None or not math.isfinite(p) or p <= 0:
        return "unresolved", False, False
    if p < 5e-8:
        return "Tier_A_genome_wide", True, True
    if p < 1e-5:
        return "Tier_B_suggestive", False, True
    return "Tier_C_other", False, False

def safe_float(value: str) -> Optional[float]:
    s = (value or "").strip().replace("×", "x").replace("−", "-")
    if not s:
        return None
    s = re.sub(r"\s+", "", s)
    s = s.replace("x10^", "e").replace("x10−", "e-").replace("x10-", "e-")
    s = re.sub(r"([0-9.]+)10\^?(-?\d+)", r"\1e\2", s)
    m = re.search(r"(?<![A-Za-z])([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)", s)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None

def fetch_bytes(url: str, timeout: int = 45) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "IS-Analysis-V3-muscle-stage1/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

def flatten_text(elem: ET.Element) -> str:
    return " ".join("".join(elem.itertext()).split())

def table_matrix(table: ET.Element) -> list[list[str]]:
    out: list[list[str]] = []
    for tr in table.findall(".//tr"):
        cells = [flatten_text(x) for x in list(tr) if x.tag.split("}")[-1] in {"td", "th"}]
        if cells:
            out.append(cells)
    return out

def find_header_index(headers: list[str], patterns: Iterable[str]) -> Optional[int]:
    norm = [re.sub(r"[^a-z0-9]+", "", h.lower()) for h in headers]
    for pat in patterns:
        p = re.sub(r"[^a-z0-9]+", "", pat.lower())
        for i, h in enumerate(norm):
            if p and p in h:
                return i
    return None

def extract_2026_variants(xml_bytes: bytes) -> tuple[list[VariantRow], list[list[str]]]:
    root = ET.fromstring(xml_bytes)
    candidate_tables = []
    for table in root.findall(".//table-wrap"):
        txt = flatten_text(table)
        if any(rs in txt for rs in SEED_2026_RSIDS):
            mat = table_matrix(table)
            if mat:
                candidate_tables.append(mat)

    rows_by_rsid = {}
    raw_table = []
    for mat in candidate_tables:
        if len(mat) < 2:
            continue
        if not raw_table:
            raw_table = mat
        headers = mat[0]
        rs_idx = find_header_index(headers, ["SNP", "rsid", "variant"])
        p_idx = find_header_index(headers, ["p value", "p-value", "pvalue", "p"])
        gene_idx = find_header_index(headers, ["gene", "nearest gene", "mapped gene"])
        chr_idx = find_header_index(headers, ["chr", "chromosome"])
        pos_idx = find_header_index(headers, ["position", "pos", "bp"])
        beta_idx = find_header_index(headers, ["beta", "effect size"])
        ea_idx = find_header_index(headers, ["effect allele", "ea", "a1"])
        oa_idx = find_header_index(headers, ["other allele", "nea", "a2"])

        for cells in mat[1:]:
            joined = " | ".join(cells)
            rsids = re.findall(r"rs\d+", joined)
            for rsid in rsids:
                if rsid not in SEED_2026_RSIDS:
                    continue
                def get(idx):
                    return cells[idx] if idx is not None and idx < len(cells) else ""
                p = safe_float(get(p_idx)) if p_idx is not None else None
                tier, gw, sugg = classify_p(p)
                rows_by_rsid[rsid] = VariantRow(
                    rsid=rsid, p_value=p, source="PMC13371584_table",
                    raw_row=joined, reported_gene=get(gene_idx), chr=get(chr_idx),
                    pos=get(pos_idx), effect_allele=get(ea_idx), other_allele=get(oa_idx),
                    beta=get(beta_idx), tier=tier, genome_wide_5e8=gw, suggestive_1e5=sugg,
                )

    rows = []
    for rsid in SEED_2026_RSIDS:
        row = rows_by_rsid.get(rsid, VariantRow(rsid=rsid))
        if row.p_value is None:
            row.tier, row.genome_wide_5e8, row.suggestive_1e5 = classify_p(None)
        rows.append(row)
    return rows, raw_table

def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def run(outdir: Path, offline: bool = False, xml_input: Optional[Path] = None) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    rawdir = outdir / "raw"
    rawdir.mkdir(exist_ok=True)

    write_tsv(outdir / "MUSCLE_STAGE1_STUDY_REGISTRY.tsv", STUDIES, list(STUDIES[0].keys()))

    xml_bytes = None
    fetch_status = "offline"
    if xml_input:
        xml_bytes = xml_input.read_bytes()
        fetch_status = f"local:{xml_input}"
    elif not offline:
        try:
            xml_bytes = fetch_bytes(EUROPE_PMC_XML)
            (rawdir / "PMC13371584.xml").write_bytes(xml_bytes)
            fetch_status = "downloaded"
        except Exception as e:
            fetch_status = f"download_failed:{type(e).__name__}:{e}"

    if xml_bytes:
        try:
            variants, raw_table = extract_2026_variants(xml_bytes)
        except Exception as e:
            variants = [VariantRow(rsid=x) for x in SEED_2026_RSIDS]
            raw_table = []
            fetch_status += f";parse_failed:{type(e).__name__}:{e}"
    else:
        variants = [VariantRow(rsid=x) for x in SEED_2026_RSIDS]
        raw_table = []

    var_dicts = [asdict(v) for v in variants]
    write_tsv(outdir / "MUSCLE_STAGE1_2026_VARIANTS.tsv", var_dicts, list(asdict(variants[0]).keys()))

    if raw_table:
        maxcols = max(len(x) for x in raw_table)
        padded = [x + [""] * (maxcols - len(x)) for x in raw_table]
        with (outdir / "MUSCLE_STAGE1_2026_RAW_TABLE.tsv").open("w", encoding="utf-8", newline="") as f:
            csv.writer(f, delimiter="\t").writerows(padded)

    numeric = [v for v in variants if v.p_value is not None]
    gw = [v for v in numeric if v.genome_wide_5e8]
    suggestive = [v for v in numeric if v.suggestive_1e5]
    unresolved = [v for v in variants if v.p_value is None]

    flags = []
    if numeric and len(gw) < len(numeric):
        flags.append({
            "severity": "HIGH",
            "flag": "reported_vs_recomputed_significance",
            "detail": f"Article reports genome-wide significance at 5e-8, but recomputation finds {len(gw)}/{len(numeric)} parsed lead rows below 5e-8.",
        })
    if unresolved:
        flags.append({
            "severity": "MEDIUM",
            "flag": "unresolved_numeric_p_values",
            "detail": f"Could not parse numeric P values for {len(unresolved)} of 9 seed rsIDs; inspect RAW_TABLE.tsv/XML before inference.",
        })
    flags.append({
        "severity": "MEDIUM",
        "flag": "possible_cohort_overlap_2024_2026",
        "detail": "Overlapping authors/grant/training design warrant participant-level overlap verification before treating studies as independent replication.",
    })
    flags.append({
        "severity": "INFO",
        "flag": "2024_summary_statistics_access",
        "detail": "2024 full dataset is described as available from Population Health Data Archive on reasonable/approved request; do not fabricate unavailable summary statistics.",
    })
    write_tsv(outdir / "MUSCLE_STAGE1_QC_FLAGS.tsv", flags, ["severity", "flag", "detail"])

    summary = {
        "stage": "MUSCLE_STAGE1_RT_GWAS_AUDIT",
        "status": "PASS_WITH_QC" if flags else "PASS",
        "europe_pmc_fetch_status": fetch_status,
        "studies_registered": len(STUDIES),
        "seed_2026_rsids": len(SEED_2026_RSIDS),
        "parsed_numeric_p_values": len(numeric),
        "recomputed_genome_wide_p_lt_5e8": len(gw),
        "recomputed_suggestive_p_lt_1e5": len(suggestive),
        "unresolved_p_values": len(unresolved),
        "next_stage_gate": "Proceed to variant-to-gene mapping only after numeric/allele/build fields are verified; coloc/MR require full summary statistics.",
    }
    (outdir / "MUSCLE_STAGE1_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    readme = f"""# MUSCLE Stage 1 RT-GWAS Audit

Status: **{summary['status']}**

- Studies registered: {summary['studies_registered']}
- 2026 lead rsIDs: {summary['seed_2026_rsids']}
- Numeric P values parsed: {summary['parsed_numeric_p_values']}
- Recomputed P < 5e-8: {summary['recomputed_genome_wide_p_lt_5e8']}
- Recomputed P < 1e-5: {summary['recomputed_suggestive_p_lt_1e5']}
- Unresolved P values: {summary['unresolved_p_values']}
- Europe PMC: `{fetch_status}`

## Gate
{summary['next_stage_gate']}

Inspect `MUSCLE_STAGE1_QC_FLAGS.tsv` before promoting any locus.
"""
    (outdir / "MUSCLE_STAGE1_README.md").write_text(readme, encoding="utf-8")
    return summary

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--offline", action="store_true")
    p.add_argument("--xml-input", type=Path)
    return p

def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    summary = run(args.outdir, args.offline, args.xml_input)
    print(json.dumps(summary, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
