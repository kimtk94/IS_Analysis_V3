#!/usr/bin/env python3
"""Read-only, fail-closed publication evidence inventory for IS and CKD.

This tool reads existing gene-level tables. It never re-estimates GWAS, MR, LD,
SuSiE, or colocalization, and it cannot convert a numeric result into a
scientifically validated causal claim. Only derived, non-person-level fields
are exported; source hashes allow independent reproduction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

CKD_GENES = (
    "ACP1", "CPVL", "F12", "GSTA1", "GSTA3",
    "HLA-E", "INHBC", "SDCCAG8", "UMOD",
)
IS_COLUMNS = (
    "gene_id_stable", "gene_symbol", "region_ids",
    "ancestry_context", "qtl_source_status",
    "archived_abf_assays", "alleles_coharmonized",
    "matched_qtl_ld_verified", "coloc_recomputed",
    "causal_claim_status", "scientific_audit_status",
    "limitations",
)
CKD_COLUMNS = (
    "gene_symbol", "eur_abf_h4", "eur_susie_max_h4",
    "eas_egfr_abf_h4", "eas_bun_abf_h4",
    "kidney_pqtl_maintext", "kidney_eqtl_meta686_hits",
    "kidney_eqtl_tubule356_hits", "kidney_eqtl_glomerulus303_hits",
    "top_cell_type", "dual_egfr_status", "dual_bun_status",
    "technical_status", "koges_anchor_present",
    "koges_longitudinal_validated", "scientific_audit_status",
    "limitations",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def source_info(path: Path, nrows: int) -> dict:
    return {"path": str(path.resolve()), "sha256": sha256(path),
            "size_bytes": path.stat().st_size, "data_rows": nrows}


def read_tsv(path: Path, required: tuple[str, ...]) -> list[dict]:
    if not path.is_file():
        raise ValueError(f"Source missing or not a file: {path}")
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        if len(header) != len(set(header)):
            raise ValueError(f"Duplicate TSV columns in {path}")
        missing = set(required) - set(header)
        if missing:
            raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No rows in {path}")
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"Malformed TSV row / column count in {path}")
    return rows


def present(value: str | None) -> bool:
    return (str(value or "").strip().upper()
            not in {"", "NA", "N/A", "NAN", "NULL", "NONE", "MISSING"})


def truth(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "yes", "1"}


def validate_posterior(value: str | None, gene: str, col: str) -> str:
    if not present(value):
        return "NOT_ASSESSED"
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{gene} {col}: not a probability") from exc
    if not 0 <= number <= 1:
        raise ValueError(f"{gene} {col}: posterior outside [0,1]")
    return str(value).strip()


def prepare_is(gene_path: Path, expected: int, region_path: Path | None = None,
               manifest_path: Path | None = None) -> tuple[list[dict], dict]:
    genes = read_tsv(gene_path, ("gene_id_stable",))
    if len(genes) != expected:
        raise ValueError(f"IS gene count {len(genes)} != expected {expected}")
    ids = [row["gene_id_stable"].strip() for row in genes]
    if len(set(ids)) != len(ids) or any(not re.fullmatch(r"ENSG[0-9]{11}", x) for x in ids):
        raise ValueError("Invalid/duplicate IS stable Ensembl gene IDs")

    regions = None
    counts = None
    if region_path is not None:
        regions = read_tsv(region_path, ("gene_id_stable",))
        if expected == 2225 and len(regions) != 2425:
            raise ValueError(f"IS region pairs {len(regions)} != frozen 2425")
        counts = Counter(r["gene_id_stable"].strip() for r in regions)
        if set(counts) != set(ids):
            raise ValueError("IS regional table omits or adds gene IDs")
        for g in genes:
            if present(g.get("n_positional_regions")) and int(g["n_positional_regions"]) != counts[g["gene_id_stable"]]:
                raise ValueError(f"IS regional multiplicity mismatch: {g['gene_id_stable']}")

    if manifest_path is not None:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if int(manifest.get("stable_gene_ids", -1)) != expected:
            raise ValueError("IS manifest gene count mismatch")
        if regions is not None and int(manifest.get("gene_region_pairs", -1)) != len(regions):
            raise ValueError("IS manifest gene-region count mismatch")

    output = []
    for g in genes:
        gene = g["gene_id_stable"].strip()
        eas = int(g.get("n_EAS_Japanese_regions") or 0)
        eur = int(g.get("n_EUR_regions") or 0)
        if eas < 0 or eur < 0:
            raise ValueError(f"Negative ancestry count: {gene}")
        ancestry = ("EAS_AND_EUR" if eas and eur else
                    "EAS_ONLY" if eas else "EUR_ONLY" if eur else "NOT_ASSESSED")
        flags = []
        allele = truth(g.get("GWAS_and_QTL_alleles_coharmonized"))
        ld = truth(g.get("matched_QTL_LD_verified"))
        coloc = truth(g.get("coloc_recomputed"))
        if not allele:
            flags.append("GWAS_QTL_ALLELE_LINK_NOT_VERIFIED")
        if not ld:
            flags.append("QTL_DONOR_LD_NOT_VERIFIED")
        if not coloc:
            flags.append("COLOC_NOT_RECOMPUTED_OR_NOT_ATTESTED")
        if counts is None:
            flags.append("GENE_REGION_MULTIPLICITY_UNVERIFIED")
        if g.get("ischemic_stroke_causal_gene", "").upper() not in ("", "NOT_ESTABLISHED"):
            flags.append("SOURCE_CAUSAL_LABEL_REQUIRES_INDEPENDENT_REVIEW")
        output.append({
            "gene_id_stable": gene,
            "gene_symbol": g.get("canonical_GENCODE_v19_gene_symbol") or g.get("canonical_gene_symbol") or "",
            "region_ids": g.get("region_ids") or g.get("region_id") or "",
            "ancestry_context": ancestry,
            "qtl_source_status": g.get("new_nominal_QTL_gene_source_status") or "NOT_ASSESSED",
            "archived_abf_assays": g.get("archived_EAS_GTEx_ABF_assays") or g.get("archived_646_coloc_assays") or "0",
            "alleles_coharmonized": str(allele).upper(),
            "matched_qtl_ld_verified": str(ld).upper(),
            "coloc_recomputed": str(coloc).upper(),
            "causal_claim_status": "NOT_ESTABLISHED",
            "scientific_audit_status": "SOURCE_DATA_PRESENT_UNVERIFIED",
            "limitations": ";".join(flags) or "MANUSCRIPT_CLAIM_REVIEW_REQUIRED",
        })
    output.sort(key=lambda x: x["gene_id_stable"])
    info = {"gene_rows": len(genes), "region_pairs": len(regions) if regions is not None else None,
            "allele_attested_rows": sum(x["alleles_coharmonized"] == "TRUE" for x in output),
            "ld_attested_rows": sum(x["matched_qtl_ld_verified"] == "TRUE" for x in output),
            "coloc_recomputed_rows": sum(x["coloc_recomputed"] == "TRUE" for x in output),
            "causal_genes_established_by_this_audit": 0}
    return output, info


def prepare_ckd(path: Path) -> tuple[list[dict], dict]:
    rows = read_tsv(path, ("gene_symbol", "eur_abf_h4",
                           "eas_egfr_abf_h4", "top_cell_type"))
    found = [g["gene_symbol"].strip().upper() for g in rows]
    if len(found) != len(CKD_GENES) or set(found) != set(CKD_GENES) or len(set(found)) != len(found):
        raise ValueError("CKD integrated matrix must have exactly the 9 frozen candidate genes")
    by_gene = {g["gene_symbol"].strip().upper(): g for g in rows}
    out = []
    for gene in CKD_GENES:
        r = by_gene[gene]
        posterior_cols = ("eur_abf_h4", "eur_susie_max_h4",
                          "eas_egfr_abf_h4", "eas_bun_abf_h4")
        probabilities = {key: validate_posterior(r.get(key), gene, key)
                         for key in posterior_cols}
        flags = ["CAUSAL_GENE_NOT_PROVEN", "KOGES_LONGITUDINAL_NOT_VERIFIED"]
        if probabilities["eur_abf_h4"] == "NOT_ASSESSED":
            flags.append("EUR_REGIONAL_COLOC_NOT_ASSESSED")
        if probabilities["eas_egfr_abf_h4"] == "NOT_ASSESSED":
            flags.append("EAS_EGFR_COLOC_NOT_ASSESSED")
        if not present(r.get("top_cell_type")):
            flags.append("KIDNEY_CELL_CONTEXT_NOT_ASSESSED")
        statuses = [r.get("dual_egfr_status", ""), r.get("dual_bun_status", ""),
                    r.get("technical_status", "")]
        if any(("UNRESOLVED" in str(x).upper() or "FAIL" in str(x).upper())
               for x in statuses):
            flags.append("MULTISIGNAL_TECHNICAL_UNRESOLVED")
        out.append({
            "gene_symbol": gene,
            **probabilities,
            "kidney_pqtl_maintext": r.get("kidney_pqtl_maintext") or "NOT_ASSESSED",
            "kidney_eqtl_meta686_hits": r.get("kidney_eqtl_meta686_hits") or "NOT_ASSESSED",
            "kidney_eqtl_tubule356_hits": r.get("kidney_eqtl_tubule356_hits") or "NOT_ASSESSED",
            "kidney_eqtl_glomerulus303_hits": r.get("kidney_eqtl_glomerulus303_hits") or "NOT_ASSESSED",
            "top_cell_type": r.get("top_cell_type") or "NOT_ASSESSED",
            "dual_egfr_status": r.get("dual_egfr_status") or "NOT_ASSESSED",
            "dual_bun_status": r.get("dual_bun_status") or "NOT_ASSESSED",
            "technical_status": r.get("technical_status") or "NOT_ASSESSED",
            "koges_anchor_present": str(present(r.get("koges_rsid"))).upper(),
            "koges_longitudinal_validated": "FALSE",
            "scientific_audit_status": "SOURCE_DATA_PRESENT_UNVERIFIED",
            "limitations": ";".join(flags),
        })
    return out, {"candidate_rows": len(out),
                 "eur_abf_present_rows": sum(x["eur_abf_h4"] != "NOT_ASSESSED" for x in out),
                 "eas_egfr_abf_present_rows": sum(x["eas_egfr_abf_h4"] != "NOT_ASSESSED" for x in out),
                 "controlled_koges_genotype_verified": False}


def write_tsv(path: Path, rows: list[dict], columns: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def freeze(is_input: Path | None, ckd_input: Path | None, out: Path,
           is_regions: Path | None = None, is_manifest: Path | None = None,
           expected_is: int = 2225) -> dict:
    if is_input is None and ckd_input is None:
        raise ValueError("Provide at least one source matrix")
    if is_input is None and (is_regions is not None or is_manifest is not None):
        raise ValueError("IS region/manifest require --is-input")
    inputs = [x for x in (is_input, ckd_input, is_regions, is_manifest) if x]
    for f in inputs:
        if not f.is_file():
            raise ValueError(f"Missing input {f}")
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory already populated: {out}")
    for f in inputs:
        if f.resolve() == out.resolve():
            raise ValueError("Output cannot be a source")
    result, details = {}, {}
    if is_input is not None:
        result["IS"], details["IS"] = prepare_is(is_input, expected_is, is_regions, is_manifest)
    if ckd_input is not None:
        result["CKD"], details["CKD"] = prepare_ckd(ckd_input)
    manifest = {
        "kind": "SOURCE_PROVENANCE_AND_COVERAGE_ONLY",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_validation_completed": False,
        "server_inputs_verified_for_this_run": True,
        "sources": {p.name: source_info(p, len(read_tsv(p, ("gene_id_stable",) if p == is_input or p == is_regions else
                                                   ("gene_symbol",) if p == ckd_input else ())))
                    for p in inputs if p != is_manifest},
        "other_input_hashes": {p.name: sha256(p) for p in inputs if p == is_manifest},
        "coverage": details,
        "limitations": [
            "Source tables and hashes verified locally; scientific assumptions not reanalysed.",
            "Recorded H4 values do not prove causal proteins or genes.",
            "A KoGES anchor rsID never demonstrates longitudinal model validation.",
            "IS operational candidate counts are not causal rankings."
        ]
    }
    out.mkdir(parents=True, exist_ok=True)
    if "IS" in result:
        write_tsv(out / "IS_EVIDENCE_READINESS_FROZEN.tsv", result["IS"], IS_COLUMNS)
    if "CKD" in result:
        write_tsv(out / "CKD_EVIDENCE_READINESS_FROZEN.tsv", result["CKD"], CKD_COLUMNS)
    (out / "EVIDENCE_FREEZE_PROVENANCE.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--is-input", type=Path)
    a.add_argument("--is-regions", type=Path)
    a.add_argument("--is-manifest", type=Path)
    a.add_argument("--ckd-input", type=Path)
    a.add_argument("--expected-is", type=int, default=2225)
    a.add_argument("--out-dir", type=Path, required=True)
    x = a.parse_args()
    if x.expected_is <= 0:
        a.error("--expected-is must be > 0")
    try:
        info = freeze(x.is_input, x.ckd_input, x.out_dir, x.is_regions,
                      x.is_manifest, x.expected_is)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        a.exit(2, f"FAIL_CLOSED: {exc}\n")
    print(json.dumps(info["coverage"], ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
