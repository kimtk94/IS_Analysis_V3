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
import importlib.util
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

AUDIT_SCHEMA_REV = "20261018-p3a-official-gwas-provenance-v8"
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
    inputs = [EAS5_VARIANTS, EAS5_READINESS, EAS5_OVERLAP,
        AUDITS / "is_eas5_98_egene_fullcis_gwas_overlap_20261011_v3_fastaanchored" /
        "IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz"]
    optional_inputs = [
        OUTPUT / "science_inputs" / "P1_APPROVED_QTL_LD.tsv",
        OUTPUT / "science_inputs" / "P2_APPROVED_GRCH37_REFERENCE.json",
        OUTPUT / "science_inputs" / "P3_APPROVED_COLOC_INPUTS.tsv",
        OUTPUT / "science_inputs" / "P3A_APPROVED_GWAS_SOURCES.tsv",
    ]
    items = []
    for path in inputs:
        if path.is_file():
            stat = path.stat()
            items.append([str(path), stat.st_size, stat.st_mtime_ns])
        else:
            items.append([str(path), "MISSING"])
    for path in optional_inputs:
        if path.is_file():
            stat = path.stat()
            items.append([str(path), stat.st_size, stat.st_mtime_ns])
        else:
            items.append([str(path), "MISSING"])
    # Changes to the science implementation must force an immediate new audit.
    stage_script = Path(__file__).resolve().with_name("harness_science_pipeline.py")
    items.append([str(stage_script), sha256_path(stage_script) if stage_script.is_file() else "MISSING"])
    archived_script = Path(__file__).resolve().with_name("harness_archived_evidence.py")
    items.append([str(archived_script),sha256_path(archived_script) if archived_script.is_file() else "MISSING"])
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


def source_exact_match(value):
    """Interpret source booleans correctly; unknown encodings fail closed."""
    raw=str(value or "").strip()
    low=raw.lower()
    if low in {"true","1","yes"}: return True
    if low in {"false","0","no","","na","none"}: return False
    import re
    if re.fullmatch(r"(?:chr)?(?:[0-9]{1,2}|X|Y|MT):[0-9]+:[ACGT]+:[ACGT]+",raw,re.I):
        return True
    raise ValueError("Unknown source exact-SNP field encoding")


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
        if source_exact_match(row["BBJ_exact_SNP"]) and source_exact_match(row["EAS_AIS_exact_SNP"]):
            group["dual_gwas_exact_snp_rows"] += 1
        if str(row["both_GWAS_sources_nonpalindromic"]).strip().lower() in {"true", "1", "yes"}:
            group["both_nonpalindromic_rows"] += 1
        if source_exact_match(row["BBJ_exact_SNP"]) and source_exact_match(row["EAS_AIS_exact_SNP"]):
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


def audit_p1_ld_gates() -> dict:
    """Review EAS5 molecular-QTL ancestry/LD readiness; no external LD is inferred."""
    statuses = collections.Counter()
    per_region = collections.defaultdict(collections.Counter)
    representative = []
    count = 0
    for _, row in safe_table(EAS5_READINESS, REQUIRED_READINESS_COLUMNS):
        if row is None:
            continue
        count += 1
        region = row["region_id"]
        key = per_region[region]
        key["gene_tissue_rows"] += 1
        ancestry = str(row["QTL_source_ancestry_verified"]).strip().lower() in {"true","1","yes"}
        ld = str(row["molecular_QTL_matched_LD_verified"]).strip().lower() in {"true","1","yes"}
        direction = str(row.get("can_now_assert_OFFICIAL_ALT_EFFECT_DIRECTION", "")).strip().lower() in {"true","1","yes"}
        if ancestry:key["source_ancestry_verified"] += 1
        if ld:key["molecular_qtl_ld_verified"] += 1
        if direction:key["official_alt_direction_assertable"] += 1
        statuses[row["scientific_status"] or "MISSING_STATUS"] += 1
        if len(representative)<20 and not ancestry and not ld:
            representative.append({
                "region_id":region,"QTL_dataset":row["QTL_dataset"],
                "gene_id":row["gene_id"],"status":row["scientific_status"]
            })
    if not count:
        raise ValueError("LD readiness table empty")
    return {
        "schema":"is-p1-ld-feasibility-v1", "generated_at":now(), "integrity":"PASS",
        "rows":count,"scientific_status_counts":dict(statuses),
        "by_region":dict(sorted((k,dict(v)) for k,v in per_region.items())),
        "illustrative_unverified":representative,
        "interpretation":"Review-only cohort/source LD gates; no LD reference, matrix or fine-mapping generated.",
    }


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
    p4_candidates = []
    rows = 0
    for row in reader:
        if None in row or None in row.values():
            raise ValueError("Malformed 2,225 candidate TSV")
        gene = row["gene_id_stable"].strip()
        if not gene.startswith("ENSG") or gene in unique:
            raise ValueError("Duplicate or invalid Ensembl gene identifier")
        unique.add(gene)
        rows += 1
        status = row["new_nominal_QTL_gene_source_status"]
        statuses[status] += 1
        causal_status[row["ischemic_stroke_causal_gene"]] += 1
        if status == "RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC":
            raw_count = row.get("n_new_nominal_QTL_gene_variant_rows", "")
            try:
                variant_rows = int(raw_count)
                if variant_rows < 0: raise ValueError("negative rows")
            except (ValueError, TypeError):
                raise ValueError("Invalid nominal QTL candidate row count")
            reported_p = row.get("min_source_nominal_QTL_p", "")
            p_present,p_valid = valid_numeric(reported_p,pvalue=True)
            if p_present and not p_valid:
                raise ValueError("Invalid nominal QTL candidate minimum p")
            p4_candidates.append({
                "gene_id": gene,
                "gene_symbol": row.get("canonical_gene_symbol",""),
                "region_id":row.get("region_id",""),
                "source_nominal_QTL_ids":row.get("source_nominal_QTL_ids",""),
                "nominal_variant_rows":variant_rows,
                "nominal_min_p":float(reported_p) if p_present else None,
                "alleles_coharmonized": row["GWAS_and_QTL_alleles_coharmonized"],
                "matched_qtl_ld_verified":row["matched_QTL_LD_verified"],
                "coloc_recomputed":row["coloc_recomputed"],
                "interpretation":"nominal availability only; not a validated stroke target",
            })
        for key in ("coloc_recomputed", "GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified"):
            if row[key].strip().lower() in {"true", "1", "yes"}:
                flags[key] += 1
    if rows != 2225 or dict(statuses) != EXPECTED_READY:
        raise ValueError(f"Snapshot integrity failed (rows={rows}, statuses={dict(statuses)})")
    p4_candidates.sort(key=lambda r: (
        r["nominal_min_p"] is None, r["nominal_min_p"] if r["nominal_min_p"] is not None else float("inf"),
        -r["nominal_variant_rows"], r["gene_id"]
    ))
    p4_detail = {
        "schema":"is-p4-candidate-triage-v1","generated_at":now(),"integrity":"PASS",
        "screened_gene_count":rows,"nominal_qtl_candidate_count":len(p4_candidates),
        "candidates":p4_candidates,
        "classification":"Descriptive triage sorted by nominal p; no FDR/MR/coloc causal ranking",
    }
    return {
        "schema": "is-candidate-audit-v1",
        "generated_at": now(), "integrity": "PASS",
        "source_commit": FROZEN_COMMIT,
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "unique_ensembl_genes": rows,
        "readiness_status_counts": dict(sorted(statuses.items())),
        "flags_true": dict(sorted(flags.items())),
        "causal_gene_status_counts": dict(sorted(causal_status.items())),
        "p4_detail":p4_detail,
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


def manifest(eas5: dict, cand: dict, p1: dict | None = None, input_signature: str | None = None) -> dict:
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
        "input_metadata_signature": input_signature if input_signature is not None else input_metadata_signature(),
        "research_qc": {
            "p1_ld": {
                "readiness_rows": p1["rows"] if p1 else None,
                "matched_ld_verified_rows": sum(v.get("molecular_qtl_ld_verified",0) for v in p1["by_region"].values()) if p1 else None,
                "ancestry_verified_rows": sum(v.get("source_ancestry_verified",0) for v in p1["by_region"].values()) if p1 else None,
                "scope": "QTL cohort ancestry/LD evidence flags, no new LD reference computed",
            },
            "p2_qc": {
                "region_gene_dataset_groups":eas5["p2_detail"]["group_count"],
                "invalid_qtl_p":eas5["p2_detail"]["qc"].get("QTL_p_invalid_nonmissing",0),
                "invalid_bbj_p":eas5["p2_detail"]["qc"].get("BBJ_p_invalid_nonmissing",0),
                "invalid_eas_p":eas5["p2_detail"]["qc"].get("EAS_AIS_p_invalid_nonmissing",0),
                "scope":"Descriptive existing summary rows, not a new allele FASTA validation",
            },
            "p4_candidates": {
                "nominal_qtl_candidate_count":cand["p4_detail"]["nominal_qtl_candidate_count"],
                "top_gene_symbols":[r["gene_symbol"] for r in cand["p4_detail"]["candidates"][:10]],
                "scope":"Nominal p ordering only, not a validated gene discovery or causal ranking",
            },
        },
        "scientific_limitations": [
            "Counts are data-readiness and QC measures, not causal gene evidence.",
            "Nominal QTL presence does not establish allele harmonization or colocalization.",
            "EAS5 variant counts may represent repeated gene/tissue measurements, not unique loci.",
            "No MR, fine-mapping or coloc was recomputed by this report.",
        ],
    }


def perform_archived_fasta_abf_audit():
    script=Path(__file__).resolve().with_name("harness_archived_evidence.py")
    if not script.is_file():raise RuntimeError("Archived FASTA audit module not installed")
    spec=importlib.util.spec_from_file_location("is_archived_fasta_audit",script)
    if not spec or not spec.loader:raise RuntimeError("Archived audit module not loadable")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result=module.audit(AUDITS)
    atomic_json(OUTPUT/"science"/"archived_fasta_abf_audit.json",result)
    return result


def perform_staged_science(eas5, candidates):
    script=Path(__file__).resolve().with_name("harness_science_pipeline.py")
    if not script.is_file():
        raise RuntimeError("Staged science module missing: deployment incomplete")
    spec=importlib.util.spec_from_file_location("is_harness_science_pipeline",script)
    if not spec or not spec.loader:
        raise RuntimeError("Staged science module cannot load")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.advance(EAS5_VARIANTS,EAS5_READINESS,candidates["p4_detail"],output=OUTPUT)


def p1_completed_snapshot() -> dict:
    """Merge independently audited P1 evidence into the existing 5-minute manifest."""
    folder = AUDITS / "is_p1_full_universe_gh_39124908298_attempt_1"
    provenance_file = folder / "P1_GITHUB_RUN_PROVENANCE.json"
    if not provenance_file.is_file():
        return {"status":"P1_RUN_MISSING","independent_replication":"NOT_ESTABLISHED"}
    provenance=json.loads(provenance_file.read_text())
    if provenance.get("status")!="P1_COMPLETE_NONCAUSAL":
        raise ValueError("P1 status unexpected")
    outputs=provenance.get("output_hashes",{})
    required={"IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv",
              "IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv",
              "IS_2225_WORKFLOW_PRIORITY_MANIFEST.json"}
    if not required.issubset(outputs):
        raise ValueError("P1 missing output hash")
    for name in sorted(required):
        fp=folder/name
        if not fp.is_file() or sha256_path(fp)!=outputs[name]:
            raise ValueError("P1 output content changed")
    summary=json.loads((folder/"IS_2225_WORKFLOW_PRIORITY_MANIFEST.json").read_text())
    if summary.get("total_source_unique_genes")!=2225 or summary.get("total_source_gene_region_pairs")!=2425:
        raise ValueError("P1 2225/2425 drift")
    return {"status":"P1_COMPLETE_NONCAUSAL","source_run":39124908298,
            "genes":2225,"gene_region_pairs":2425,"tier_counts":summary.get("tier_counts"),
            "top40_count":len(summary.get("TOP40_non_G0022_workflow_genes",[])),
            "independent_replication":"NOT_ESTABLISHED",
            "gwas_source_overlap":"BBJ_PRESENT_IN_GIGASTROKE_EAS_DISCOVERY",
            "causal_gene_status":"NOT_ESTABLISHED",
            "manifest_sha256":sha256_path(folder/"IS_2225_WORKFLOW_PRIORITY_MANIFEST.json")}


def dispatch(op: str) -> dict:
    if op == "eas5-audit":
        e = audit_eas5()
        p1 = audit_p1_ld_gates()
        atomic_json(OUTPUT / "eas5_audit.json", {k:v for k,v in e.items() if k != "p2_detail"})
        atomic_json(OUTPUT / "p2_region_qc.json", e["p2_detail"])
        atomic_json(OUTPUT / "p1_ld_gate.json", p1)
        return {"operation": op, "integrity": "PASS", "rows": e["rows"],
                "regions": e["regions"], "ld_verified_rows": e["matched_qtl_ld_verified_rows"],
                "coloc_recomputed_rows": e["coloc_recomputed_rows"],
                "p2_groups":e["p2_detail"]["group_count"], "p1_review_rows":p1["rows"]}
    if op == "candidate-2225-audit":
        c = audit_candidates()
        atomic_json(OUTPUT / "candidate_2225_audit.json", {k:v for k,v in c.items() if k != "p4_detail"})
        atomic_json(OUTPUT / "p4_nominal_candidate_triage.json", c["p4_detail"])
        return {"operation": op, "integrity": "PASS",
                "unique_genes": c["unique_ensembl_genes"],
                "nominal_triage_genes":c["p4_detail"]["nominal_qtl_candidate_count"],
                "status_counts": c["readiness_status_counts"]}
    if op == "figure-manifest":
        e = load_validated(OUTPUT / "eas5_audit.json", "is-eas5-audit-v1")
        c = load_validated(OUTPUT / "candidate_2225_audit.json", "is-candidate-audit-v1")
        e["p2_detail"] = load_validated(OUTPUT / "p2_region_qc.json", "is-p2-eas5-qc-v1")
        c["p4_detail"] = load_validated(OUTPUT / "p4_nominal_candidate_triage.json", "is-p4-candidate-triage-v1")
        p1 = load_validated(OUTPUT / "p1_ld_gate.json", "is-p1-ld-feasibility-v1")
        m = manifest(e, c, p1)
        atomic_json(OUTPUT / "research_manifest.json", m)
        return m
    if op == "research-refresh":
        e = audit_eas5()
        c = audit_candidates()
        p1 = audit_p1_ld_gates()
        science = perform_staged_science(e,c)
        m = manifest(e, c, p1)
        archived_module=Path(__file__).resolve().with_name("harness_archived_evidence.py")
        if archived_module.is_file():
            archive=perform_archived_fasta_abf_audit()
            m["archived_fasta_abf"]={
                "schema":archive["schema"],
                "archived":archive["archived"],
                "overlap_report":archive["overlap_report"],
                "readiness":archive["readiness"],
                "abf_screen":archive["abf_screen"],
                "abf_gate_summary":archive["abf_gate_summary"],
                "p3a_source_header_audit":archive["p3a_source_header_audit"],
                "p3a_gwas_source_identity":archive["p3a_gwas_source_identity"],
                "scientific_status":"ARCHIVED_SOURCE_FLAGS_AND_ABF_READINESS_ONLY",
            }
        else:
            m["archived_fasta_abf"]={"scientific_status":"MODULE_NOT_INSTALLED_NOT_AUDITED"}
        m["p1_full_universe"] = p1_completed_snapshot()
        m["science_pipeline"] = {
            "schema":science["schema"],
            "p1":science["p1"],
            "p2":science["p2"],
            "p3":science["p3"],
            "p4":science["p4"],
            "p2_reference_state":science["p2_reference_state"],
            "scientific_result":science["scientific_result"],
        }
        atomic_json(OUTPUT / "eas5_audit.json", {k:v for k,v in e.items() if k != "p2_detail"})
        atomic_json(OUTPUT / "p2_region_qc.json", e["p2_detail"])
        atomic_json(OUTPUT / "candidate_2225_audit.json", {k:v for k,v in c.items() if k != "p4_detail"})
        atomic_json(OUTPUT / "p4_nominal_candidate_triage.json", c["p4_detail"])
        atomic_json(OUTPUT / "p1_ld_gate.json", p1)
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
