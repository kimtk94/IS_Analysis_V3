#!/usr/bin/env python3
"""Audit public GEO human ischemic stroke single-cell series metadata.

Purpose: choose a disease-context substitute for ROSMAP expression annotation,
NOT an eQTL substitute. Downloads only GEO SOFT metadata (not cell count
matrices, FASTQ, identifiable donor material or restricted data).

Checks GSE331114 (brain CD45+ enriched), GSE267023 (PBMC phases),
GSE285659 (LVO intravascular PBMC). Writes fail-closed per-study evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import io
import json
import os
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request

SERIES = {
    "GSE331114": {
        "role": "human_acute_IS_brain_immune_cell_context",
        "species": "Homo sapiens",
        "expected_sample_count": 8,
        "expected_groups": {"AIS": 4, "Control": 4},
        "expected_tissue_context": "brain_CD45_enriched",
        "genotype_mapped_eqtl": False,
    },
    "GSE267023": {
        "role": "human_acute_chronic_IS_PBMC_context",
        "species": "Homo sapiens",
        "expected_sample_count": 16,
        "expected_groups": None,
        "expected_tissue_context": "PBMC",
        "genotype_mapped_eqtl": False,
    },
    "GSE285659": {
        "role": "human_LVO_arterial_compartment_PBMC",
        "species": "Homo sapiens",
        "expected_sample_count": None,
        "expected_groups": None,
        "expected_tissue_context": "stroke_intravascular_PBMC",
        "genotype_mapped_eqtl": False,
    },
}
MAX_COMPRESSED = 16 * 1024 * 1024
MAX_EXPANDED = 36 * 1024 * 1024


def geo_soft_url(gse):
    if re.fullmatch(r"GSE[0-9]{4,8}", gse) is None:
        raise ValueError("Invalid GSE accession")
    bucket = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}/soft/{gse}_family.soft.gz"


def download_metadata(url):
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != "ftp.ncbi.nlm.nih.gov":
        raise ValueError("Non-NCBI metadata URL")
    req = urllib.request.Request(url, headers={
        "User-Agent": "IS-Analysis-V3-StrokeGEO-AccessAudit/1.0"
    })
    with urllib.request.urlopen(req, timeout=35) as resp:
        final = urllib.parse.urlparse(resp.geturl())
        if final.scheme != "https" or final.hostname != "ftp.ncbi.nlm.nih.gov":
            raise ValueError("Unexpected metadata redirect origin")
        chunks = []
        total = 0
        while True:
            part = resp.read(256 * 1024)
            if not part:
                break
            total += len(part)
            if total > MAX_COMPRESSED:
                raise ValueError("Compressed metadata size cap exceeded")
            chunks.append(part)
    raw = gzip.decompress(b"".join(chunks))
    if len(raw) > MAX_EXPANDED:
        raise ValueError("Expanded metadata size cap exceeded")
    return raw.decode("utf-8", errors="replace"), total


def parse_soft(text):
    """Read per-sample metadata only, preserving sample as unit of analysis."""
    series = {}
    samples = {}
    current = None
    for line in io.StringIO(text):
        if line.startswith("^SERIES = "):
            series["accession"] = line.split(" = ", 1)[1].strip()
            current = "series"
        elif line.startswith("^SAMPLE = "):
            sid = line.split(" = ", 1)[1].strip()
            if sid in samples:
                raise ValueError("Duplicated GSM entry")
            samples[sid] = {"accession": sid, "characteristics": [], "files": [], "title": ""}
            current = sid
        elif line.startswith("^"):
            current = None
        elif line.startswith("!Series_") and current == "series":
            if line.startswith("!Series_title = "):
                series["title"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Series_summary = "):
                series.setdefault("summaries", []).append(line.split(" = ", 1)[1].strip()[:1000])
            elif line.startswith("!Series_overall_design = "):
                series.setdefault("design", []).append(line.split(" = ", 1)[1].strip()[:1000])
        elif current in samples:
            sample = samples[current]
            if line.startswith("!Sample_title = "):
                sample["title"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                sample["source_name"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_organism_ch1 = "):
                sample["organism"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                sample["characteristics"].append(line.split(" = ", 1)[1].strip())
            elif line.startswith("!Sample_supplementary_file = "):
                sample["files"].append(line.split(" = ", 1)[1].strip())
    return series, samples


def derive_group(sample):
    for entry in sample["characteristics"]:
        match = re.match(r"^group\s*:\s*(.*?)\s*$", entry, flags=re.IGNORECASE)
        if match:
            name = match.group(1).lower()
            if name in {"ais", "ischemic stroke", "stroke"}:
                return "AIS"
            if name in {"control", "non-stroke", "healthy"}:
                return "Control"
            return "OTHER"
    title = sample["title"].lower()
    if "ais_brain" in title:
        return "AIS"
    if "control_brain" in title:
        return "Control"
    return "UNKNOWN"


def audit_one(accession, expected):
    url = geo_soft_url(accession)
    outcome = {
        "accession": accession,
        "role": expected["role"],
        "metadata_url": url,
        "metadata_status": "SOURCE_UNAVAILABLE",
        "expected_tissue": expected["expected_tissue_context"],
        "genotype_mapped_eqtl_in_submission": False,
        "full_variant_level_eqtl": "NOT_PROVIDED_BY_THIS_SERIES",
        "actual_expression_data_downloaded": False,
        "individual_sample_count_validated": False,
        "patient_donor_level_independence": "NOT_VERIFIED",
        "GWAS_locus_colocalization": "NOT_APPLICABLE_WITHOUT_GENOTYPES",
    }
    try:
        text, size = download_metadata(url)
        info, sample_records = parse_soft(text)
        outcome["compressed_soft_bytes"] = size
        outcome["title"] = info.get("title")
        outcome["sample_count"] = len(sample_records)
        outcome["sample_ids"] = sorted(sample_records)
        organism = Counter(x.get("organism", "UNKNOWN") for x in sample_records.values())
        outcome["organism_counts"] = dict(organism)
        group_counts = Counter(derive_group(x) for x in sample_records.values())
        outcome["sample_groups"] = dict(group_counts)
        outcome["supplementary_file_count"] = sum(len(x.get("files", [])) for x in sample_records.values())
        outcome["mtx_file_count_in_soft"] = sum(
            f.endswith((".mtx.gz", ".mtx")) for x in sample_records.values() for f in x.get("files", [])
        )
        outcome["metadata_status"] = (
            "PUBLIC_SOFT_METADATA_PASS" if info.get("accession") == accession and
            len(sample_records) > 0 and set(organism) == {expected["species"]}
            else "METADATA_REVIEW"
        )
        if expected["expected_sample_count"] is not None:
            outcome["individual_sample_count_validated"] = (
                len(sample_records) == expected["expected_sample_count"]
            )
            if not outcome["individual_sample_count_validated"]:
                outcome["metadata_status"] = "COUNT_MISMATCH_REVIEW"
        if expected["expected_groups"] is not None:
            outcome["group_counts_match_expected"] = (
                dict(group_counts) == expected["expected_groups"]
            )
            if not outcome["group_counts_match_expected"]:
                outcome["metadata_status"] = "GROUP_MISMATCH_REVIEW"
        outcome["selected_sample_metadata"] = [{
            "id": x["accession"],
            "title": x["title"],
            "source_name": x.get("source_name"),
            "characteristics": x["characteristics"],
            "file_extensions": sorted(set(
                (".mtx.gz" if f.endswith(".mtx.gz") else
                 ".tsv.gz" if f.endswith(".tsv.gz") else
                 ".other") for f in x.get("files", [])
            )),
        } for x in sample_records.values()]
    except (urllib.error.URLError, OSError, ValueError, TimeoutError, EOFError) as exc:
        outcome["error_type"] = type(exc).__name__
        outcome["error_message_short"] = str(exc)[:180]
    return outcome


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="artifacts/is_stroke_alternative_geo")
    args = p.parse_args()
    directory = Path(args.out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    records = [audit_one(x, y) for x, y in SERIES.items()]
    payload = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "OFFICIAL_GEO_SERIES_METADATA_ONLY",
        "do_not_call_eqtl_or_genetic_validation": True,
        "results": records,
    }
    (directory / "stroke_alternative_geo.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Human ischemic stroke alternatives to ROSMAP",
        "",
        "**Clinical disease expression context, not eQTL replacement.**",
        "",
        "| GEO | Metadata | GSM | AIS/control | Matrix file count |",
        "|---|---|---:|---|---:|",
    ]
    for row in records:
        lines.append(
            f"| {row['accession']} | {row['metadata_status']} | "
            f"{row.get('sample_count', 'unknown')} | "
            f"{row.get('sample_groups', {})} | "
            f"{row.get('mtx_file_count_in_soft', 'unknown')} |"
        )
    lines += [
        "",
        "GSE331114 was designed as human stroke brain CD45+-enriched immune cell",
        "scRNA-seq, not an unbiased entire human cortex single-nucleus atlas.",
        "The study's series description mentions skull marrow (n=7+7) as well as",
        "brain (n=4+4); **only actual GEO GSM counts in THIS series count**.",
        "All three GEO records lack verified paired donor genotypes/eQTL summary files.",
        "Post-stroke response states do not equal pre-stroke inherited susceptibility.",
    ]
    result = "\n".join(lines) + "\n"
    (directory / "summary.md").write_text(result, encoding="utf-8")
    print(result)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
