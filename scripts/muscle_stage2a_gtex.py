#!/usr/bin/env python3
"""MUSCLE Stage 2A v2: normalize RT GWAS variants and query GTEx skeletal-muscle QTLs.

Key rule:
GTEx /dataset/variant accepts rsIDs (snpId), but the static association endpoints
are queried with GTEx variant IDs (variantId). Therefore each rsID is first
resolved to one or more GTEx variant IDs, and those IDs are used for eQTL/sQTL
queries.

A variant not present in GTEx is not an API error. It is recorded as
SKIPPED_NO_GTEX_VARIANT for the QTL steps.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

VERSION = "2A.2"
ENSEMBL_BASE = "https://rest.ensembl.org"
GTEX_BASE = "https://gtexportal.org/api/v2"
DEFAULT_DATASET = "gtex_v10"
DEFAULT_TISSUE = "Muscle_Skeletal"

REQUIRED_STAGE1_COLUMNS = {
    "study_id", "rsid", "chr", "pos", "gene", "p", "recalculated_tier",
    "functional_evidence_points", "stage1_priority",
}

VARIANT_FIELDS = [
    "study_id", "rsid", "source_gene", "source_chr", "source_pos",
    "source_p", "stage1_tier", "stage1_points", "stage1_priority",
    "ensembl_grch38_chr", "ensembl_grch38_pos", "ensembl_allele_string",
    "source_position_matches_grch38",
    "gtex_variant_ids", "gtex_b37_variant_ids",
    "source_position_matches_gtex_b37",
    "gtex_variant_records", "gtex_eqtl_count", "gtex_sqtl_count",
    "gtex_eqtl_genes", "gtex_sqtl_genes",
    "ensembl_status", "gtex_variant_status", "gtex_eqtl_status",
    "gtex_sqtl_status",
]

EQTL_FIELDS = [
    "study_id", "rsid", "source_gene", "stage1_priority",
    "gencode_id", "gene_symbol", "variant_id", "snp_id",
    "tissue", "dataset_id", "p_value", "q_value", "nes",
    "slope", "slope_se", "maf",
]

SQTL_FIELDS = [
    "study_id", "rsid", "source_gene", "stage1_priority",
    "gencode_id", "gene_symbol", "variant_id", "snp_id",
    "phenotype_id", "tissue", "dataset_id", "p_value", "q_value",
    "nes", "slope", "slope_se", "maf",
]

CANDIDATE_FIELDS = [
    "gene", "source_gene_flag", "gtex_eqtl_flag", "gtex_sqtl_flag",
    "supporting_rsids", "supporting_studies", "n_supporting_variants",
    "best_stage1_points", "best_stage1_priority", "best_source_p",
    "stage2a_evidence_score",
]


def _missing(value: Any) -> bool:
    return value is None or str(value).strip() in {"", "NA", "N/A", ".", "None", "null"}


def _float(value: Any, default: Optional[float] = None) -> Optional[float]:
    if _missing(value):
        return default
    try:
        x = float(value)
    except (TypeError, ValueError):
        return default
    return x if math.isfinite(x) else default


def _first(record: Dict[str, Any], *names: str, default: Any = "") -> Any:
    for name in names:
        if name in record and record[name] is not None:
            return record[name]
    return default


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_STAGE1_COLUMNS - fields)
        if missing:
            raise ValueError(f"Missing Stage 1 columns: {', '.join(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def write_tsv(path: Path, rows: Iterable[Dict[str, Any]], fields: Sequence[str]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")


def http_get_json(
    base_url: str,
    params: Optional[Dict[str, Any]] = None,
    timeout: int = 30,
    retries: int = 3,
    pause: float = 1.0,
) -> Any:
    url = base_url
    if params:
        clean = {k: v for k, v in params.items() if v not in (None, "")}
        url = f"{base_url}?{urlencode(clean, doseq=True)}"

    headers = {
        "Accept": "application/json",
        "User-Agent": "IS_Analysis_V3-MUSCLE-Stage2A/2.0",
    }

    last: Optional[Exception] = None
    last_body = ""
    for attempt in range(1, retries + 1):
        try:
            with urlopen(Request(url, headers=headers), timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            last = exc
            try:
                last_body = exc.read().decode("utf-8", errors="replace")
            except Exception:
                last_body = ""
            if attempt < retries and exc.code >= 500:
                time.sleep(pause * attempt)
                continue
            break
        except (URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            if attempt < retries:
                time.sleep(pause * attempt)

    raise RuntimeError(
        f"GET failed: {url}: {last}; body={last_body[:1000]}"
    )


def extract_records(payload: Any, preferred_keys: Sequence[str] = ()) -> List[Dict[str, Any]]:
    if isinstance(payload, list):
        return [x for x in payload if isinstance(x, dict)]
    if not isinstance(payload, dict):
        return []

    for key in ("data", *preferred_keys):
        value = payload.get(key)
        if isinstance(value, list):
            return [x for x in value if isinstance(x, dict)]
        if isinstance(value, dict):
            for nested in ("data", "items"):
                nv = value.get(nested)
                if isinstance(nv, list):
                    return [x for x in nv if isinstance(x, dict)]

    for value in payload.values():
        if isinstance(value, list) and all(isinstance(x, dict) for x in value):
            return value
    return []


def select_grch38_mapping(payload: Dict[str, Any]) -> Dict[str, Any]:
    mappings = payload.get("mappings", []) if isinstance(payload, dict) else []
    valid_chroms = {str(i) for i in range(1, 23)} | {"X", "Y"}
    candidates = []
    for m in mappings:
        if not isinstance(m, dict):
            continue
        if m.get("assembly_name") != "GRCh38":
            continue
        if m.get("coord_system") != "chromosome":
            continue
        chrom = str(m.get("seq_region_name", "")).replace("chr", "")
        if chrom not in valid_chroms:
            continue
        candidates.append(m)
    if not candidates:
        return {}
    candidates.sort(key=lambda x: (str(x.get("seq_region_name", "")), int(x.get("start", 0))))
    return candidates[0]


def query_ensembl(rsid: str) -> Tuple[Dict[str, Any], str]:
    try:
        payload = http_get_json(
            f"{ENSEMBL_BASE}/variation/human/{rsid}",
            params={"content-type": "application/json"},
        )
        if not isinstance(payload, dict):
            return {}, "BAD_RESPONSE"
        return payload, "OK"
    except Exception as exc:
        return {"error": str(exc)}, "ERROR"


def query_gtex_variant(rsid: str, dataset: str) -> Tuple[Any, str]:
    try:
        payload = http_get_json(
            f"{GTEX_BASE}/dataset/variant",
            params={
                "snpId": rsid,
                "datasetId": dataset,
                "itemsPerPage": 2000,
            },
        )
        return payload, "OK"
    except Exception as exc:
        return {"error": str(exc)}, "ERROR"


def query_gtex_qtl(
    variant_id: str,
    dataset: str,
    tissue: str,
    qtl: str,
) -> Tuple[Any, str]:
    endpoint = "singleTissueEqtl" if qtl == "eqtl" else "singleTissueSqtl"
    try:
        payload = http_get_json(
            f"{GTEX_BASE}/association/{endpoint}",
            params={
                "variantId": variant_id,
                "tissueSiteDetailId": tissue,
                "datasetId": dataset,
                "itemsPerPage": 10000,
            },
        )
        return payload, "OK"
    except Exception as exc:
        return {"error": str(exc)}, "ERROR"


def normalize_qtl_record(
    stage1: Dict[str, str], rec: Dict[str, Any], qtl: str
) -> Dict[str, Any]:
    out = {
        "study_id": stage1["study_id"],
        "rsid": stage1["rsid"],
        "source_gene": stage1["gene"],
        "stage1_priority": stage1["stage1_priority"],
        "gencode_id": _first(rec, "gencodeId", "gencode_id"),
        "gene_symbol": _first(rec, "geneSymbol", "gene_symbol"),
        "variant_id": _first(rec, "variantId", "variant_id"),
        "snp_id": _first(rec, "snpId", "snp_id", default=stage1["rsid"]),
        "tissue": _first(rec, "tissueSiteDetailId", "tissueSiteDetail", default=DEFAULT_TISSUE),
        "dataset_id": _first(rec, "datasetId", "dataset_id"),
        "p_value": _first(rec, "pValue", "pval", "p_value"),
        "q_value": _first(rec, "qValue", "qval", "q_value"),
        "nes": _first(rec, "nes"),
        "slope": _first(rec, "slope"),
        "slope_se": _first(rec, "slopeStandardError", "slopeSe", "slope_se"),
        "maf": _first(rec, "maf"),
    }
    if qtl == "sqtl":
        out["phenotype_id"] = _first(
            rec, "phenotypeId", "phenotype_id", "intronId", "junctionId"
        )
    return out


def split_source_genes(value: str) -> List[str]:
    genes = []
    for token in str(value or "").replace(",", ";").split(";"):
        token = token.strip()
        if token and token.upper() not in {"NA", "N/A", "."}:
            genes.append(token)
    return genes


def parse_gtex_b37_variant_id(value: str) -> Optional[Tuple[str, int]]:
    if not value:
        return None
    m = re.match(r"^(?:chr)?([^_]+)_(\d+)_.*_b37$", str(value))
    if not m:
        return None
    return m.group(1).replace("chr", ""), int(m.group(2))


def source_matches_any_b37(source_chr: str, source_pos: int, b37_ids: Sequence[str]) -> str:
    for vid in b37_ids:
        parsed = parse_gtex_b37_variant_id(vid)
        if not parsed:
            continue
        chrom, pos = parsed
        if chrom.upper() == source_chr.upper() and pos == source_pos:
            return "1"
    return "0" if b37_ids else ""


def build_candidate_genes(
    stage1_rows: List[Dict[str, str]],
    eqtl_rows: List[Dict[str, Any]],
    sqtl_rows: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    by_rsid = {r["rsid"]: r for r in stage1_rows}
    evidence: Dict[str, Dict[str, Any]] = defaultdict(
        lambda: {"source": set(), "eqtl": set(), "sqtl": set(), "rsids": set(), "studies": set()}
    )

    for row in stage1_rows:
        for gene in split_source_genes(row["gene"]):
            e = evidence[gene]
            e["source"].add(row["rsid"])
            e["rsids"].add(row["rsid"])
            e["studies"].add(row["study_id"])

    for qtl_rows, key in ((eqtl_rows, "eqtl"), (sqtl_rows, "sqtl")):
        for row in qtl_rows:
            gene = str(row.get("gene_symbol") or "").strip()
            if not gene:
                gene = str(row.get("gencode_id") or "").split(".")[0].strip()
            if not gene:
                continue
            e = evidence[gene]
            e[key].add(row["rsid"])
            e["rsids"].add(row["rsid"])
            e["studies"].add(row["study_id"])

    priority_order = {
        "HIGH_FUNCTIONAL_PRIORITY": 3,
        "MEDIUM_FUNCTIONAL_PRIORITY": 2,
        "EXPLORATORY": 1,
    }
    out = []
    for gene, e in evidence.items():
        rsids = sorted(e["rsids"])
        supporting = [by_rsid[r] for r in rsids if r in by_rsid]
        best_points = max((int(r["functional_evidence_points"]) for r in supporting), default=0)
        best_priority = max(
            (r["stage1_priority"] for r in supporting),
            key=lambda x: priority_order.get(x, 0),
            default="",
        )
        pvals = [_float(r["p"]) for r in supporting]
        pvals = [p for p in pvals if p is not None]

        source_flag = bool(e["source"])
        eqtl_flag = bool(e["eqtl"])
        sqtl_flag = bool(e["sqtl"])
        score = best_points + (1 if source_flag else 0) + (3 if eqtl_flag else 0) + (2 if sqtl_flag else 0)

        out.append({
            "gene": gene,
            "source_gene_flag": int(source_flag),
            "gtex_eqtl_flag": int(eqtl_flag),
            "gtex_sqtl_flag": int(sqtl_flag),
            "supporting_rsids": ";".join(rsids),
            "supporting_studies": ";".join(sorted(e["studies"])),
            "n_supporting_variants": len(rsids),
            "best_stage1_points": best_points,
            "best_stage1_priority": best_priority,
            "best_source_p": min(pvals) if pvals else "",
            "stage2a_evidence_score": score,
        })

    out.sort(
        key=lambda r: (
            -int(r["stage2a_evidence_score"]),
            -int(r["gtex_eqtl_flag"]),
            -int(r["gtex_sqtl_flag"]),
            _float(r["best_source_p"], 1.0) or 1.0,
            r["gene"],
        )
    )
    return out


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def run(
    input_path: Path,
    outdir: Path,
    dataset: str = DEFAULT_DATASET,
    tissue: str = DEFAULT_TISSUE,
    pause: float = 0.15,
) -> Dict[str, Any]:
    stage1_rows = read_tsv(input_path)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_root = outdir / "raw"
    variant_rows: List[Dict[str, Any]] = []
    eqtl_rows: List[Dict[str, Any]] = []
    sqtl_rows: List[Dict[str, Any]] = []

    for idx, row in enumerate(stage1_rows, start=1):
        rsid = row["rsid"]
        print(f"[{idx}/{len(stage1_rows)}] {rsid} {row['gene']}")

        ens_payload, ens_status = query_ensembl(rsid)
        write_json(raw_root / "ensembl" / f"{rsid}.json", ens_payload)
        mapping = select_grch38_mapping(ens_payload if isinstance(ens_payload, dict) else {})

        gv_payload, gv_status = query_gtex_variant(rsid, dataset)
        write_json(raw_root / "gtex_variant" / f"{rsid}.json", gv_payload)
        gv_records = extract_records(gv_payload, ("variant", "variants"))

        variant_ids = sorted({
            str(_first(x, "variantId", "variant_id")).strip()
            for x in gv_records
            if str(_first(x, "variantId", "variant_id")).strip()
        })
        b37_ids = sorted({
            str(_first(x, "b37VariantId", "b37_variant_id")).strip()
            for x in gv_records
            if str(_first(x, "b37VariantId", "b37_variant_id")).strip()
        })

        normalized_eq: List[Dict[str, Any]] = []
        normalized_sq: List[Dict[str, Any]] = []
        eq_statuses: List[str] = []
        sq_statuses: List[str] = []

        if not variant_ids:
            eq_status = "SKIPPED_NO_GTEX_VARIANT"
            sq_status = "SKIPPED_NO_GTEX_VARIANT"
            write_json(
                raw_root / "gtex_eqtl" / f"{rsid}.json",
                {"status": eq_status, "rsid": rsid, "variantIds": []},
            )
            write_json(
                raw_root / "gtex_sqtl" / f"{rsid}.json",
                {"status": sq_status, "rsid": rsid, "variantIds": []},
            )
        else:
            eq_payloads = []
            sq_payloads = []
            for variant_id in variant_ids:
                eq_payload, one_eq_status = query_gtex_qtl(
                    variant_id, dataset, tissue, "eqtl"
                )
                sq_payload, one_sq_status = query_gtex_qtl(
                    variant_id, dataset, tissue, "sqtl"
                )
                eq_statuses.append(one_eq_status)
                sq_statuses.append(one_sq_status)
                eq_payloads.append({
                    "variantId": variant_id,
                    "status": one_eq_status,
                    "payload": eq_payload,
                })
                sq_payloads.append({
                    "variantId": variant_id,
                    "status": one_sq_status,
                    "payload": sq_payload,
                })

                eq_records = extract_records(eq_payload, ("singleTissueEqtl",))
                sq_records = extract_records(sq_payload, ("singleTissueSqtl",))
                normalized_eq.extend(
                    normalize_qtl_record(row, x, "eqtl") for x in eq_records
                )
                normalized_sq.extend(
                    normalize_qtl_record(row, x, "sqtl") for x in sq_records
                )

            eq_status = "ERROR" if "ERROR" in eq_statuses else "OK"
            sq_status = "ERROR" if "ERROR" in sq_statuses else "OK"
            write_json(
                raw_root / "gtex_eqtl" / f"{rsid}.json",
                {"rsid": rsid, "queries": eq_payloads},
            )
            write_json(
                raw_root / "gtex_sqtl" / f"{rsid}.json",
                {"rsid": rsid, "queries": sq_payloads},
            )

        eqtl_rows.extend(normalized_eq)
        sqtl_rows.extend(normalized_sq)

        source_chr = str(row["chr"]).replace("chr", "")
        source_pos = int(float(row["pos"]))
        grch38_chr = str(mapping.get("seq_region_name", "")).replace("chr", "")
        grch38_pos = mapping.get("start", "")
        pos_match38 = (
            "1"
            if grch38_chr and grch38_pos
            and source_chr.upper() == grch38_chr.upper()
            and source_pos == int(grch38_pos)
            else "0" if grch38_chr and grch38_pos else ""
        )

        eq_genes = sorted({
            str(x.get("gene_symbol") or x.get("gencode_id") or "").strip()
            for x in normalized_eq
            if str(x.get("gene_symbol") or x.get("gencode_id") or "").strip()
        })
        sq_genes = sorted({
            str(x.get("gene_symbol") or x.get("gencode_id") or "").strip()
            for x in normalized_sq
            if str(x.get("gene_symbol") or x.get("gencode_id") or "").strip()
        })

        variant_rows.append({
            "study_id": row["study_id"],
            "rsid": rsid,
            "source_gene": row["gene"],
            "source_chr": source_chr,
            "source_pos": source_pos,
            "source_p": row["p"],
            "stage1_tier": row["recalculated_tier"],
            "stage1_points": row["functional_evidence_points"],
            "stage1_priority": row["stage1_priority"],
            "ensembl_grch38_chr": grch38_chr,
            "ensembl_grch38_pos": grch38_pos,
            "ensembl_allele_string": mapping.get("allele_string", ""),
            "source_position_matches_grch38": pos_match38,
            "gtex_variant_ids": ";".join(variant_ids),
            "gtex_b37_variant_ids": ";".join(b37_ids),
            "source_position_matches_gtex_b37": source_matches_any_b37(
                source_chr, source_pos, b37_ids
            ),
            "gtex_variant_records": len(gv_records),
            "gtex_eqtl_count": len(normalized_eq),
            "gtex_sqtl_count": len(normalized_sq),
            "gtex_eqtl_genes": ";".join(eq_genes),
            "gtex_sqtl_genes": ";".join(sq_genes),
            "ensembl_status": ens_status,
            "gtex_variant_status": gv_status,
            "gtex_eqtl_status": eq_status,
            "gtex_sqtl_status": sq_status,
        })

        if pause > 0:
            time.sleep(pause)

    candidate_rows = build_candidate_genes(stage1_rows, eqtl_rows, sqtl_rows)

    write_tsv(outdir / "MUSCLE_STAGE2A_VARIANT_NORMALIZATION.tsv", variant_rows, VARIANT_FIELDS)
    write_tsv(outdir / "MUSCLE_STAGE2A_EQTL.tsv", eqtl_rows, EQTL_FIELDS)
    write_tsv(outdir / "MUSCLE_STAGE2A_SQTL.tsv", sqtl_rows, SQTL_FIELDS)
    write_tsv(outdir / "MUSCLE_STAGE2A_CANDIDATE_GENES.tsv", candidate_rows, CANDIDATE_FIELDS)

    api_error_variants = sorted({
        r["rsid"]
        for r in variant_rows
        if "ERROR" in {
            r["ensembl_status"], r["gtex_variant_status"],
            r["gtex_eqtl_status"], r["gtex_sqtl_status"],
        }
    })
    unresolved_gtex = sorted({
        r["rsid"] for r in variant_rows if int(r["gtex_variant_records"]) == 0
    })

    summary = {
        "version": VERSION,
        "dataset": dataset,
        "tissue": tissue,
        "n_input_variants": len(stage1_rows),
        "n_ensembl_grch38_resolved": sum(bool(r["ensembl_grch38_chr"]) for r in variant_rows),
        "n_source_position_matches_grch38": sum(
            r["source_position_matches_grch38"] == "1" for r in variant_rows
        ),
        "n_source_position_matches_gtex_b37": sum(
            r["source_position_matches_gtex_b37"] == "1" for r in variant_rows
        ),
        "n_gtex_variant_resolved": sum(int(r["gtex_variant_records"]) > 0 for r in variant_rows),
        "n_gtex_variant_unresolved": len(unresolved_gtex),
        "gtex_variant_unresolved_rsids": unresolved_gtex,
        "n_variants_with_significant_eqtl": sum(int(r["gtex_eqtl_count"]) > 0 for r in variant_rows),
        "n_variants_with_significant_sqtl": sum(int(r["gtex_sqtl_count"]) > 0 for r in variant_rows),
        "n_eqtl_associations": len(eqtl_rows),
        "n_sqtl_associations": len(sqtl_rows),
        "n_candidate_genes": len(candidate_rows),
        "top_candidate_genes": [
            {
                "gene": r["gene"],
                "score": r["stage2a_evidence_score"],
                "rsids": r["supporting_rsids"],
                "eqtl": r["gtex_eqtl_flag"],
                "sqtl": r["gtex_sqtl_flag"],
            }
            for r in candidate_rows[:20]
        ],
        "n_api_error_variants": len(api_error_variants),
        "api_error_variants": api_error_variants,
        "interpretation": (
            "GTEx associations are queried by GTEx variantId after rsID resolution. "
            "Unresolved GTEx variants are not API errors. Zero QTL results mean no "
            "significant precomputed association returned for the resolved variant IDs."
        ),
    }
    write_json(outdir / "MUSCLE_STAGE2A_SUMMARY.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path,
                   help="Stage 1 MUSCLE_STAGE1_GWAS_AUDIT.tsv")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--dataset", default=DEFAULT_DATASET)
    p.add_argument("--tissue", default=DEFAULT_TISSUE)
    p.add_argument("--pause", type=float, default=0.15)
    return p


def main() -> None:
    args = build_parser().parse_args()
    run(args.input, args.outdir, args.dataset, args.tissue, args.pause)


if __name__ == "__main__":
    main()
