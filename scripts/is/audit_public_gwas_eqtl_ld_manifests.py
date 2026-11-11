#!/usr/bin/env python3
"""Audit externally published manifests for GIGASTROKE, MAEEA eQTL and 1KG reference.

Metadata-only, bounded HTTP. Downloading a Zenodo file or performing analyses
is explicitly NOT performed. Missing or ambiguous fields remain HOLD.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import ssl
import urllib.error
import urllib.parse
import urllib.request

ZENODO = {
    "MAEEA_EAS_CIS_EQTL": "https://zenodo.org/api/records/21296030",
    "ONE_KG_GRCH37_EUR_EAS_REFERENCE": "https://zenodo.org/api/records/7768714",
}
FTP_ROOT = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90104001-GCST90105000")
GWAS = {
    "EUR_IS": "GCST90104540",
    "EAS_IS": "GCST90104545",
}
CAP = 2_000_000


def safe_fetch_json_or_text(url, cap=CAP):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "IS-Analysis-V3-ResearchSourceManifestAudit/1.0",
                 "Accept": "application/json,text/plain,*/*"},
    )
    with urllib.request.urlopen(req, timeout=25) as resp:
        real_url = urllib.parse.urlparse(resp.geturl())
        intended = urllib.parse.urlparse(url)
        if real_url.scheme != "https" or real_url.hostname != intended.hostname:
            raise ValueError("Source redirected outside trusted origin")
        length = resp.headers.get("Content-Length")
        if length and int(length) > cap:
            raise ValueError("Metadata exceeds size cap")
        raw = resp.read(cap + 1)
        if len(raw) > cap:
            raise ValueError("Metadata response exceeds cap")
        return raw.decode("utf-8", errors="replace"), resp.status


def zenodo_metadata(label, url):
    rec = {
        "dataset": label,
        "record_url": url.replace("/api/records/", "/records/"),
        "api_url": url,
        "access_status": "SOURCE_UNAVAILABLE",
        "files_checked_for_content": False,
        "full_variant_level_eQTL_QC": "HOLD" if "EQTL" in label else "NOT_APPLICABLE",
        "reference_LD_matched_and_QC": "HOLD" if "ONE_KG" in label else "NOT_APPLICABLE",
    }
    try:
        raw, code = safe_fetch_json_or_text(url)
        payload = json.loads(raw)
        rec["access_status"] = "PUBLIC_METADATA_ACCESSIBLE" if code == 200 else "HTTP_REVIEW"
        rec["record_id"] = payload.get("id")
        meta = payload.get("metadata") or {}
        rec["title"] = meta.get("title")
        rec["publication_date"] = meta.get("publication_date")
        rec["resource_type"] = meta.get("resource_type")
        rights = meta.get("rights") or {}
        rec["rights"] = rights if isinstance(rights, str) else rights.get("id")
        rec["versions"] = payload.get("revision")
        files = payload.get("files", [])
        result_files = []
        for f in files:
            links = f.get("links") or {}
            link = links.get("self") if isinstance(links, dict) else None
            if link:
                parsed = urllib.parse.urlparse(link)
                if parsed.scheme != "https" or parsed.hostname not in {"zenodo.org", "www.zenodo.org"}:
                    link = None
            checksum = f.get("checksum", "")
            result_files.append({
                "key": f.get("key"),
                "size_bytes": f.get("size"),
                "checksum_publisher": checksum,
                "file_url": link,
                "content_status": "NOT_DOWNLOADED",
            })
        rec["file_count"] = len(result_files)
        rec["files"] = result_files
        rec["full_cis_eqtl_availability"] = (
            "FILE_MANIFEST_ONLY" if "EQTL" in label else "NOT_APPLICABLE"
        )
    except (urllib.error.URLError, TimeoutError, OSError, ValueError,
            json.JSONDecodeError) as exc:
        rec["error_type"] = type(exc).__name__
        rec["error_short"] = str(exc)[:200]
    return rec


def scan_yaml_text(text):
    fields = {}
    # Scalar values only, no full YAML parser required for chosen fields.
    for key in ("gwas_id", "data_file_name", "data_file_md5sum",
                "genome_assembly", "case_count", "control_count",
                "sample_size", "sample_ancestry", "trait_description", "file_type"):
        match = re.search(r"^\s*" + re.escape(key) + r":\s*(.*?)\s*$", text,
                          flags=re.MULTILINE | re.IGNORECASE)
        if match:
            val = match.group(1).strip().strip("'\"")
            fields[key] = val[:250]
    return fields


def gwas_yaml_manifest(label, accession):
    file = f"{accession}_buildGRCh37.tsv.gz"
    url = f"{FTP_ROOT}/{accession}/{file}-meta.yaml"
    rec = {
        "dataset": label,
        "accession": accession,
        "metadata_url": url,
        "status": "SOURCE_UNAVAILABLE",
        "original_file_checksum_match": "NOT_TESTED",
    }
    try:
        raw, code = safe_fetch_json_or_text(url, cap=500_000)
        fields = scan_yaml_text(raw)
        rec["metadata_fields"] = fields
        if fields.get("gwas_id") == accession and fields.get("genome_assembly") in ("GRCh37", "hg19"):
            rec["status"] = "META_ACCESSION_BUILD_MATCH"
        elif fields:
            rec["status"] = "META_REVIEW"
        else:
            rec["status"] = "META_EMPTY_OR_UNSUPPORTED"
        rec["publisher_md5_reported"] = bool(re.fullmatch(
            r"[0-9a-fA-F]{32}", fields.get("data_file_md5sum", "")
        ))
    except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
        rec["error_type"] = type(exc).__name__
        rec["error_short"] = str(exc)[:200]
    return rec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="artifacts/is_eas_public_source_manifests")
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    zenodo = [zenodo_metadata(name, url) for name, url in ZENODO.items()]
    gwases = [gwas_yaml_manifest(name, acc) for name, acc in GWAS.items()]
    report = {
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "metadata_only": True,
        "eqtl_statistical_feasibility": "HOLD",
        "matched_ld_quality": "HOLD",
        "gwas_full_file_independent_checksum": "NOT_TESTED",
        "zenodo": zenodo, "gwas_metadata": gwases,
    }
    (out / "public_source_manifest_audit.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = ["# Public eQTL / LD / GWAS source metadata discovery", "",
             "Metadata discovery is not verification of full cis-eQTL SNP availability,",
             "effect/SE format, variant coverage, biological significance, or ancestry-matched LD.",
             "",
             "| Source | Metadata status | Files | Actual analysis |",
             "|---|---|---:|---|"]
    for r in zenodo:
        lines.append(f"| {r['dataset']} | {r['access_status']} | "
                     f"{r.get('file_count', 'unknown')} | HOLD |")
    for r in gwases:
        lines.append(f"| {r['dataset']} YAML | {r['status']} | - | HOLD checksum |")
    lines += ["", "## File manifest (source metadata only)", ""]
    for r in zenodo:
        lines.append("### " + r["dataset"])
        for f in r.get("files", []):
            lines.append(f"- {f.get('key')}: {f.get('size_bytes')} bytes; "
                         f"source checksum `{f.get('checksum_publisher')}`")
    for r in gwases:
        lines.append("### " + r["dataset"] + " YAML")
        lines.append("`" + json.dumps(r.get("metadata_fields", {}), ensure_ascii=False) + "`")
    summary = "\n".join(lines) + "\n"
    (out / "public_source_manifest_audit.md").write_text(summary, encoding="utf-8")
    print(summary)
    target = os.getenv("GITHUB_STEP_SUMMARY")
    if target:
        with open(target, "a", encoding="utf-8") as stream:
            stream.write(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
