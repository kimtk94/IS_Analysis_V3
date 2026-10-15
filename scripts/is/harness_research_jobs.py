#!/usr/bin/env python3
"""Restricted, evidence-first IS research audit jobs. Python stdlib only.

Reads existing server scientific outputs; never changes them. All generated JSON
belongs solely under /srv/is-analysis/results/is/github_actions_runs/evidence/.
"""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
import urllib.request

ROOT = Path("/srv/is-analysis/results/is")
AUDITS = ROOT / "audits"
OUTPUT = ROOT / "github_actions_runs" / "evidence"
FROZEN_COMMIT = "fdb1983c436d8f438a606c14c3322bb157349d77"
READINESS_URL = (
    "https://raw.githubusercontent.com/kimtk94/IS_Analysis_V3/"
    + FROZEN_COMMIT
    + "/docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv"
)
EAS5_VARIANTS = AUDITS / "is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1" / "IS_EAS5_OFFICIAL_QTL_ALT_GWAS_HARMONIZED_VARIANTS.tsv"
EAS5_READINESS = AUDITS / "is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1" / "IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv"
EAS5_OVERLAP = AUDITS / "is_eas5_qtl_gwas_snp_ref37_20261011_v2" / "IS_EAS5_GENE_TISSUE_GWAS_SNP_OVERLAP_QC.tsv"

AUDIT_SCHEMA_REV = "20261016-p2-p4-qc-v1"
EXPECTED_READY = {
    "NOT_IN_EAS_30_GROUPS": 1356,
    "PUBLIC_NOMINAL_QTL_REGIONAL_PILOT_NOT_EXECUTED": 629,
    "RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC": 149,
    "NO_NOMINAL_ROWS_IN_TARGETED_SLICES_NOT_NEGATIVE": 91,
}
REQUIRED_VARIANT_COLUMNS = (
    "region_id", "QTL_dataset", "gene_id",
    "QTL_beta_mapped_GRCh37_ALT_OFFICIAL_DEFINITION", "QTL_p",
    "BBJ_exact_SNP", "EAS_AIS_exact_SNP",
    "both_GWAS_sources_nonpalindromic", "status",
    "source_cohort_ancestry_verified", "matched_QTL_LD_verified",
    "coloc_recomputed"
)
REQUIRED_READINESS_COLUMNS = (
    "region_id", "QTL_dataset", "gene_id", "scientific_status",
    "QTL_source_ancestry_verified", "molecular_QTL_matched_LD_verified",
    "dual_GWAS_and_QTL_coloc_results",
)
REQUIRED_OVERLAP_COLUMNS = (
    "region_id", "QTL_dataset", "gene_id",
    "both_GWAS_exact_nonpalindromic_rows", "eligible_for_coloc", "blocked_on"
)



def input_metadata_signature() -> str:
    """Fast change hint; metadata is not content proof. Hourly full hashing remains required."""
    inputs = [EAS5_VARIANTS, EAS5_READINESS, EAS5_OVERLAP]
    items = []
    for path in inputs:
        stat = path.stat()
        items.append([str(path), stat.st_size, stat.st_mtime_ns])
    raw = json.dumps({
        "frozen_commit": FROZEN_COMMIT,
        "schema_rev": AUDIT_SCHEMA_REV,
        "files": items,
    }, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def valid_numeric(value, *, pvalue=False):
    """Return (present, valid) without interpreting missing statistics as negative evidence."""
    if str(value).strip() in ("", "NA", "NaN", "nan", "None", "."):
        return False, False
    try:
        numeric = float(value)
        return True, math.isfinite(numeric) and (not pvalue or 0 <= numeric <= 1)
    except (ValueError, TypeError):
        return True, False


def detail_from_variant_groups(group_counts, qc_flags, total):
    """P2: descriptive eligibility only; never asserts full allele validation."""
    return {
        "schema": "is-p2-eas5-qc-v1",
        "generated_at": now(),
        "rows": total,
        "group_count": len(group_counts),
        "groups": [
            {"region_id": key[0], "QTL_dataset": key[1], "gene_id": key[2], **dict(value)}
            for key, value in sorted(group_counts.items())
        ],
        "qc": dict(sorted(qc_flags.items())),
        "statement": "Existing harmonized-table descriptive QC; not a new FASTA/EAF audit or independent coloc",
        "integrity": "PASS",
    }

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_table(path: Path, required: tuple[str, ...]):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Required audit input unavailable: {path.name}")
    with path.open("r", encoding="utf-8-sig", errors="strict", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        header = reader.fieldnames or []
        if len(header) != len(set(header)):
            raise ValueError(f"Duplicate columns: {path.name}")
        missing = sorted(set(required).difference(header))
        if missing:
            raise ValueError(f"Missing {path.name} columns: {missing}")
        count = 0
        yield header, None
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"Malformed TSV row {count + 2}: {path.name}")
            count += 1
            yield None, row
        if not count:
            raise ValueError(f"No data rows: {path.name}")


def table_summary(path: Path, columns: tuple[str, ...]) -> dict:
    rows = 0
    regions = set()
    for header, row in safe_table(path, columns):
        if row is None:
            continue
        rows += 1
        regions.add(row["region_id"])
    return {
        "rows": rows,
        "regions": len(regions),
        "sha256": sha256_path(path),
        "bytes": path.stat().st_size,
    }


def audit_eas5() -> dict:
    rows = 0
    region_counts = collections.Counter()
    status_counts = collections.Counter()
    group_counts = collections.defaultdict(collections.Counter)
    qc_flags = collections.Counter()
    dual = nonpal = ld = coloc = ancestry = 0
    for header, row in safe_table(EAS5_VARIANTS, REQUIRED_VARIANT_COLUMNS):
        if row is None:
            continue
        rows += 1
        region_counts[row["region_id"]] += 1
        status_counts[row["status"] or "MISSING_STATUS"] += 1
        group = group_counts[(row["region_id"], row["QTL_dataset"], row["gene_id"])]
        group["rows"] += 1
        for field in ("QTL_p", "BBJ_p", "EAS_AIS_p"):
            present, valid = valid_numeric(row.get(field, ""), pvalue=True)
            if present and not valid:
                qc_flags[field + "_invalid_nonmissing"] += 1
                group[field + "_invalid_nonmissing"] += 1
            if not present:
                qc_flags[field + "_missing"] += 1
        for field in ("QTL_beta_mapped_GRCh37_ALT_OFFICIAL_DEFINITION",
                      "BBJ_beta_mapped_ALT", "EAS_AIS_beta_mapped_ALT"):
            present, valid = valid_numeric(row.get(field, ""))
            if present and not valid:
                qc_flags[field + "_invalid_nonmissing"] += 1
        if row["BBJ_exact_SNP"] and row["EAS_AIS_exact_SNP"]:
            group["dual_gwas_exact_snp_rows"] += 1
        if str(row["both_GWAS_sources_nonpalindromic"]).strip().lower() in {"true", "1", "yes"}:
            group["both_nonpalindromic_rows"] += 1
        if row["BBJ_exact_SNP"] and row["EAS_AIS_exact_SNP"]:
            dual += 1
        if str(row["both_GWAS_sources_nonpalindromic"]).strip().lower() in {"true", "1", "yes"}:
            nonpal += 1
        if str(row["matched_QTL_LD_verified"]).strip().lower() in {"true", "1", "yes"}:
            ld += 1
        if str(row["coloc_recomputed"]).strip().lower() in {"true", "1", "yes"}:
            coloc += 1
        if str(row["source_cohort_ancestry_verified"]).strip().lower() in {"true", "1", "yes"}:
            ancestry += 1
    if len(region_counts) != 5:
        raise ValueError(f"Expected 5 EAS regions, found {len(region_counts)}")
    if not (0 <= nonpal <= rows and 0 <= dual <= rows):
        raise AssertionError("Invalid counts")
    result = {
        "schema": "is-eas5-audit-v1",
        "generated_at": now(), "integrity": "PASS",
        "source": "existing server TSV, no recomputation",
        "rows": rows,
        "regions": len(region_counts),
        "region_rows": dict(sorted(region_counts.items())),
        "status_rows": dict(status_counts.most_common(12)),
        "dual_gwas_exact_snp_rows": dual,
        "both_nonpalindromic_rows": nonpal,
        "ancestry_verified_rows": ancestry,
        "matched_qtl_ld_verified_rows": ld,
        "coloc_recomputed_rows": coloc,
        "variants_file_sha256": sha256_path(EAS5_VARIANTS),
        "readiness_table": table_summary(EAS5_READINESS, REQUIRED_READINESS_COLUMNS),
        "overlap_qc_table": table_summary(EAS5_OVERLAP, REQUIRED_OVERLAP_COLUMNS),
        "p2_detail": detail_from_variant_groups(group_counts, qc_flags, rows),
        "causal_gene_conclusions": "NOT_ESTABLISHED_BY_THIS_AUDIT",
    }
    return result


def download_frozen_readiness() -> bytes:
    req = urllib.request.Request(READINESS_URL, headers={"User-Agent": "is-harness-evidence/1"})
    with urllib.request.urlopen(req, timeout=25) as response:
        data = response.read(2_000_001)
    if len(data) > 2_000_000 or not data:
        raise ValueError("Invalid frozen public readiness snapshot size")
    return data


def audit_candidates(data: bytes | None = None) -> dict:
    data = download_frozen_readiness() if data is None else data
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter="\t")
    required = ("gene_id_stable", "new_nominal_QTL_gene_source_status",
                "ischemic_stroke_causal_gene", "coloc_recomputed",
                "GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified")
    if not reader.fieldnames or set(required).difference(reader.fieldnames):
        raise ValueError("Readiness TSV missing required columns")
    statuses = collections.Counter()
    unique = set()
    causal_status = collections.Counter()
    flags = collections.Counter()
    rows = 0
    for row in reader:
        if None in row or None in row.values():
            raise ValueError("Malformed 2,225 candidate TSV")
        gene = row["gene_id_stable"].strip()
        if not gene.startswith("ENSG") or gene in unique:
            raise ValueError("Duplicate or invalid Ensembl gene identifier")
        unique.add(gene)
        rows += 1
        statuses[row["new_nominal_QTL_gene_source_status"]] += 1
        causal_status[row["ischemic_stroke_causal_gene"]] += 1
        for key in ("coloc_recomputed", "GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified"):
            if row[key].strip().lower() in {"true", "1", "yes"}:
                flags[key] += 1
    if rows != 2225 or dict(statuses) != EXPECTED_READY:
        raise ValueError(f"Snapshot integrity failed (rows={rows}, statuses={dict(statuses)})")
    return {
        "schema": "is-candidate-audit-v1",
        "generated_at": now(), "integrity": "PASS",
        "source_commit": FROZEN_COMMIT,
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "unique_ensembl_genes": rows,
        "readiness_status_counts": dict(sorted(statuses.items())),
        "flags_true": dict(sorted(flags.items())),
        "causal_gene_status_counts": dict(sorted(causal_status.items())),
        "causal_gene_conclusions": "NOT_ESTABLISHED_BY_THIS_AUDIT",
    }


def atomic_json(path: Path, obj: dict):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=OUTPUT, prefix=".tmp_", suffix=".json",
                                 encoding="utf-8", delete=False) as tmp:
        json.dump(obj, tmp, ensure_ascii=False, indent=2, sort_keys=True)
        tmp.write("\n")
        tmp_path = Path(tmp.name)
    try:
        os.replace(tmp_path, path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def load_validated(path: Path, schema: str) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema") != schema or value.get("integrity") != "PASS":
        raise ValueError(f"Invalid persisted audit: {path.name}")
    return value


def manifest(eas5: dict, cand: dict) -> dict:
    if eas5.get("integrity") != "PASS" or cand.get("integrity") != "PASS":
        raise ValueError("Unvalidated source audits")
    statuses = cand["readiness_status_counts"]
    bars = [
        {"label": "Nominal QTL rows present", "value": statuses["RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC"]},
        {"label": "Pilot not executed", "value": statuses["PUBLIC_NOMINAL_QTL_REGIONAL_PILOT_NOT_EXECUTED"]},
        {"label": "No rows in targeted slices", "value": statuses["NO_NOMINAL_ROWS_IN_TARGETED_SLICES_NOT_NEGATIVE"]},
        {"label": "Outside EAS 30 groups", "value": statuses["NOT_IN_EAS_30_GROUPS"]},
    ]
    return {
        "schema": "is-research-evidence-v1",
        "generated_at": now(), "integrity": "PASS",
        "figures": [
            {"id": "candidate-readiness", "title": "2,225 candidate genes: public QTL readiness",
             "denominator": cand["unique_ensembl_genes"], "bars": bars},
            {"id": "eas5-evidence", "title": "EAS5 matched QTL/GWAS variants: QC coverage",
             "denominator": eas5["rows"], "bars": [
                {"label": "BBJ + EAS exact SNP overlap", "value": eas5["dual_gwas_exact_snp_rows"]},
                {"label": "Both GWAS nonpalindromic", "value": eas5["both_nonpalindromic_rows"]},
                {"label": "Matched QTL LD verified", "value": eas5["matched_qtl_ld_verified_rows"]},
                {"label": "Colocalization recomputed", "value": eas5["coloc_recomputed_rows"]},
             ]},
        ],
        "metrics": {
            "candidate_genes": cand["unique_ensembl_genes"],
            "eas_regions": eas5["regions"],
            "eas_variant_rows": eas5["rows"],
        },
        "provenance": {
            "readiness_source_sha256": cand["source_sha256"],
            "eas5_variants_sha256": eas5["variants_file_sha256"],
        },
        "scientific_limitations": [
            "Counts are data-readiness and QC measures, not causal gene evidence.",
            "Nominal QTL presence does not establish allele harmonization or colocalization.",
            "EAS5 variant counts may represent repeated gene/tissue measurements, not unique loci.",
            "No MR, fine-mapping or coloc was recomputed by this report.",
        ],
    }


def dispatch(op: str) -> dict:
    if op == "eas5-audit":
        e = audit_eas5()
        atomic_json(OUTPUT / "eas5_audit.json", e)
        return {"operation": op, "integrity": "PASS", "rows": e["rows"],
                "regions": e["regions"], "ld_verified_rows": e["matched_qtl_ld_verified_rows"],
                "coloc_recomputed_rows": e["coloc_recomputed_rows"]}
    if op == "candidate-2225-audit":
        c = audit_candidates()
        atomic_json(OUTPUT / "candidate_2225_audit.json", c)
        return {"operation": op, "integrity": "PASS",
                "unique_genes": c["unique_ensembl_genes"],
                "status_counts": c["readiness_status_counts"]}
    if op == "figure-manifest":
        e = load_validated(OUTPUT / "eas5_audit.json", "is-eas5-audit-v1")
        c = load_validated(OUTPUT / "candidate_2225_audit.json", "is-candidate-audit-v1")
        m = manifest(e, c)
        atomic_json(OUTPUT / "research_manifest.json", m)
        return m
    if op == "research-refresh":
        e = audit_eas5()
        c = audit_candidates()
        m = manifest(e, c)
        atomic_json(OUTPUT / "eas5_audit.json", e)
        atomic_json(OUTPUT / "candidate_2225_audit.json", c)
        atomic_json(OUTPUT / "research_manifest.json", m)
        return m
    raise ValueError("Unsupported research operation")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--operation", required=True,
                        choices=["eas5-audit", "candidate-2225-audit",
                                 "figure-manifest", "research-refresh"])
    args = parser.parse_args()
    result = dispatch(args.operation)
    # Deliberately publish only a tightly scoped aggregate summary, not TSV rows.
    print("IS_EVIDENCE_JSON_V1=" + json.dumps(result, ensure_ascii=False, separators=(",", ":"), sort_keys=True))


if __name__ == "__main__":
    main()
