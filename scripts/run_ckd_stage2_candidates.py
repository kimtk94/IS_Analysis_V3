#!/usr/bin/env python3
"""Build Stage 2 CKD candidate and UKB-PPP full-summary download manifests."""
from __future__ import annotations
import argparse, csv, gzip, json, math, re
from pathlib import Path

def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("r", encoding="utf-8", newline="")

def read_tsv(path: Path):
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def fnum(x):
    try:
        v = float(str(x).strip())
        return v if math.isfinite(v) else None
    except Exception:
        return None

def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def assay_parts(protein_id: str):
    parts = protein_id.split(":")
    if len(parts) < 4:
        raise ValueError(f"unexpected protein_id: {protein_id}")
    gene, accession = parts[0], parts[1]
    oid = next((p for p in parts if re.fullmatch(r"OID\d+", p)), "")
    version = next((p for p in parts if re.fullmatch(r"v\d+", p)), "")
    if not oid or not version:
        raise ValueError(f"OID/version missing in protein_id: {protein_id}")
    return gene, accession, oid, version

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage1-root", type=Path, required=True)
    ap.add_argument("--download-manifest", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--min-support-direction", type=int, default=3)
    ap.add_argument("--min-support-p05-direction", type=int, default=2)
    args = ap.parse_args()

    summary_rows = read_tsv(args.stage1_root / "stage1_protein_summary.tsv")
    iq = read_tsv(args.stage1_root / "instrument_qc.tsv.gz")
    anchor = {(r["protein_id"], r["rsid"]): r for r in iq}

    candidates = []
    for r in summary_rows:
        if str(r.get("screen_fdr05", "")) != "1":
            continue
        if int(float(r.get("support_same_risk_direction_n") or 0)) < args.min_support_direction:
            continue
        if int(float(r.get("support_p05_same_risk_direction_n") or 0)) < args.min_support_p05_direction:
            continue
        if str(r.get("eas_egfr_available", "")) != "1":
            continue
        if str(r.get("eas_top_sign_concordant", "")) != "1":
            continue

        gene, accession, oid, version = assay_parts(r["protein_id"])
        a = anchor.get((r["protein_id"], r["eur_egfr_rsid"]))
        if not a:
            raise SystemExit(f"anchor instrument not found: {r['protein_id']} {r['eur_egfr_rsid']}")
        ep = fnum(r.get("eas_egfr_p"))
        candidates.append({
            "protein_id": r["protein_id"],
            "gene_symbol": gene,
            "uniprot": accession,
            "oid": oid,
            "assay_version": version,
            "anchor_rsid": r["eur_egfr_rsid"],
            "anchor_chr_hg19": a.get("chrom_hg19", ""),
            "anchor_pos_hg19": a.get("pos_hg19", ""),
            "anchor_ref": a.get("ref_allele", ""),
            "anchor_alt": a.get("alt_allele", ""),
            "anchor_f": r["eur_egfr_f"],
            "eur_egfr_beta": r["eur_egfr_beta"],
            "eur_egfr_p": r["eur_egfr_p"],
            "eur_egfr_fdr": r["eur_egfr_fdr"],
            "support_available_n": r["support_available_n"],
            "support_same_risk_direction_n": r["support_same_risk_direction_n"],
            "support_p05_same_risk_direction_n": r["support_p05_same_risk_direction_n"],
            "eas_egfr_beta": r["eas_egfr_beta"],
            "eas_egfr_p": r["eas_egfr_p"],
            "eas_replication_status": "nominal_p05" if ep is not None and ep < 0.05 else "direction_only",
        })
    candidates.sort(key=lambda x: float(x["eur_egfr_fdr"]))
    if not candidates:
        raise SystemExit("no Stage 2 candidates passed the gate")

    write_tsv(args.output_root / "stage2_candidates.tsv", candidates, list(candidates[0]))

    raw_manifest = read_tsv(args.download_manifest)
    plan = []
    for c in candidates:
        needle = f"_{c['uniprot']}_{c['oid']}_{c['assay_version']}_"
        matches = []
        for m in raw_manifest:
            gene = (m.get("gene_symbol") or m.get("gene") or "").upper()
            if gene == c["gene_symbol"].upper() and needle in m.get("source_file", ""):
                matches.append(m)
        by_anc = {}
        for m in matches:
            by_anc.setdefault(m.get("ancestry", "").upper(), []).append(m)
        for anc in ("EUR", "EAS"):
            rows = by_anc.get(anc, [])
            if len(rows) != 1:
                raise SystemExit(f"{c['protein_id']}: expected exactly one {anc} manifest row, got {len(rows)}")
            m = rows[0]
            plan.append({
                "protein_id": c["protein_id"],
                "gene_symbol": c["gene_symbol"],
                "uniprot": c["uniprot"],
                "oid": c["oid"],
                "ancestry": anc,
                "stage2_role": "primary_coloc" if anc == "EUR" else "eas_replication_optional",
                "download_now": 1 if anc == "EUR" else 0,
                "source_file": m.get("source_file", ""),
                "synapse_id": m.get("synapse_id", ""),
                "url": m.get("url") or m.get("source_url", ""),
                "expected_size_bytes": m.get("expected_size_bytes") or m.get("size_bytes", ""),
                "md5": m.get("md5", ""),
                "sha256": m.get("sha256", ""),
            })
    write_tsv(args.output_root / "stage2_pqtl_download_manifest.tsv", plan, list(plan[0]))

    regions = []
    for c in candidates:
        pos = int(float(c["anchor_pos_hg19"]))
        regions.append({
            "protein_id": c["protein_id"],
            "gene_symbol": c["gene_symbol"],
            "anchor_rsid": c["anchor_rsid"],
            "chrom_hg19": c["anchor_chr_hg19"],
            "anchor_pos_hg19": pos,
            "window_bp": 1000000,
            "region_start_hg19": max(1, pos - 1000000),
            "region_end_hg19": pos + 1000000,
        })
    write_tsv(args.output_root / "stage2_anchor_regions_hg19.tsv", regions, list(regions[0]))

    out_summary = {
        "stage": "CKD Stage 2 candidate planning",
        "candidate_rule": {
            "eur_egfr_fdr05": True,
            "min_support_same_risk_direction_n": args.min_support_direction,
            "min_support_p05_same_risk_direction_n": args.min_support_p05_direction,
            "eas_same_anchor_available": True,
            "eas_direction_concordant": True,
        },
        "n_candidates": len(candidates),
        "n_eas_nominal_p05": sum(c["eas_replication_status"] == "nominal_p05" for c in candidates),
        "n_download_rows": len(plan),
        "n_primary_eur_downloads": sum(p["download_now"] == 1 for p in plan),
        "genes": [c["gene_symbol"] for c in candidates],
        "note": "Selection is a follow-up gate, not evidence of causality; colocalization remains required.",
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "STAGE2_PLAN_SUMMARY.json").write_text(json.dumps(out_summary, indent=2) + "\n")
    print(json.dumps(out_summary, indent=2))
    print(f"CKD_STAGE2_PLAN_PASS output={args.output_root}")

if __name__ == "__main__":
    main()
