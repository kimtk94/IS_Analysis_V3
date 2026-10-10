#!/usr/bin/env python3
"""Audit ImmuNexUT rs671-region competing genes with JCTF molecular evidence.

This integrates source-level *marginal* associations only, not coloc, MR, or
independent causal replication. No protected genotypes or bulk archives.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path

VARIANTS = {
    "rs11066015": ("12:112168009:G:A", 111730205),
    "rs671": ("12:112241766:G:A", 111803962),
    "rs11066132": ("12:112468206:C:T", 112030402),
    "rs77768175": ("12:112736118:A:G", 112298314),
}
GENES = ("ALDH2", "BRAP", "RPH3A", "MYL2")
CELLTYPES = {"Neu": "Neutrophils", "CL_Mono": "Classical_monocytes",
             "Plasmablast": "Plasmablasts"}
URL_ROOT = "https://www.immunexut.org/eqtlSnpsLink?snp_id="


def parse_variant_page(payload, expected_rsid):
    if isinstance(payload, bytes):
        payload = payload.decode("utf-8", "replace")
    match = re.search(r"<snps-eqtl-table-component\s+eqtl-data\s*=\s*'([^']*)'", payload)
    if not match:
        raise ValueError(f"Missing browser eQTL component: {expected_rsid}")
    raw = json.loads(match.group(1))
    if not isinstance(raw, list):
        raise ValueError("Expected eQTL association list")
    if not raw:
        raise ValueError("No browser eQTL associations; NOT_TESTED vs true no association unresolved")
    if len(raw) > 1000:
        raise ValueError("Unanticipated variant association display size")
    expected_variant, grch38 = VARIANTS[expected_rsid]
    validated, seen = [], set()
    for item in raw:
        if (item.get("snp_id") != expected_rsid or
            item.get("position") != grch38 or
            item.get("chromosome") != "chr12"):
            raise ValueError(f"Wrong variant identity in {expected_rsid} record")
        gene = item["gene_symbol"]
        cell = item["cell_type"]
        key = (gene, cell)
        if key in seen:
            raise ValueError(f"Repeated gene/cell without documented independent signal: {key}")
        seen.add(key)
        p = float(item["eqtl_pval1"])
        beta = float(item["eqtl_effect_beta1"])
        if not (math.isfinite(p) and 0 < p <= 1 and math.isfinite(beta)):
            raise ValueError(f"Malformed QTL association at {expected_rsid} {gene}")
        validated.append({
            "eas_ais_credible_set_variant_grch37": expected_variant,
            "rsid": expected_rsid, "qtl_grch38_chr": "12",
            "qtl_grch38_pos": grch38, "gene": gene,
            "immunexut_celltype_code": cell,
            "immunexut_celltype": CELLTYPES.get(cell, cell),
            "immunexut_browser_p": p,
            "immunexut_browser_beta_NOT_GWAS_HARMONIZED": beta,
            "source_effect_allele_known": False,
            "source_display_selection": "FDR_LT_0_05_SIGNIFICANT_ONLY",
            "qtl_status": "SUPPORTED_MARGINAL_QTL_NOT_COLOC",
        })
    return validated


def load_jctf(path):
    evidence = {}
    with Path(path).open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = row["gene_symbol"]
            var = row["GWAS_variant_GRCh37"]
            trait = row["qtl_category"]
            if gene not in GENES or var not in {v[0] for v in VARIANTS.values()}:
                continue
            key = (gene, trait)
            j = evidence.setdefault(key, {"genes": gene, "trait": trait, "variants": set(),
                        "min_p": 1.0, "max_susie_pip": 0.0})
            j["variants"].add(var)
            j["min_p"] = min(j["min_p"], float(row["jctf_allele_association_p"]))
            j["max_susie_pip"] = max(j["max_susie_pip"], float(row["jctf_susie_variant_pip"]))
    return evidence


def build_summary(evidences, jctf):
    by_gene = {}
    for gene in GENES:
        var_records = [r for r in evidences if r["gene"] == gene]
        summarized = {
            "gene": gene,
            "immunexut_variants_associated": len({r["rsid"] for r in var_records}),
            "immunexut_celltypes_reported": ";".join(sorted({r["immunexut_celltype"] for r in var_records})),
            "immunexut_min_marginal_p": min((r["immunexut_browser_p"] for r in var_records), default=None),
            "immunexut_beta_unharmonized_same_across_variants": len({
                r["immunexut_browser_beta_NOT_GWAS_HARMONIZED"] for r in var_records}) == 1 if var_records else None,
            "jctf_eqtl_marker_count": len(jctf.get((gene, "eQTL"), {}).get("variants", set())),
            "jctf_eqtl_min_p": jctf.get((gene, "eQTL"), {}).get("min_p"),
            "jctf_pqtl_marker_count": len(jctf.get((gene, "pQTL"), {}).get("variants", set())),
            "jctf_pqtl_min_p": jctf.get((gene, "pQTL"), {}).get("min_p"),
            "jctf_eqtl_max_SuSiE_PIP": jctf.get((gene, "eQTL"), {}).get("max_susie_pip"),
            "jctf_pqtl_max_SuSiE_PIP": jctf.get((gene, "pQTL"), {}).get("max_susie_pip"),
            "same_causal_variant_established": False,
            "GWAS_QTL_full_locus_coloc_valid": False,
            "causal_gene_established": False,
            "source_rank_is_not_causal_rank": True,
        }
        by_gene[gene] = summarized
    return list(by_gene.values()), {
        "audit": "IS_G0022_IMMUNEXUT_VARIANT_CENTRIC_COMPETING_GENES_V1",
        "science_gate": "SUPPORTED_MARGINAL_EQTL_PLUS_JCTF_CROSS_RESOURCE_NOT_INDEPENDENT_CAUSAL_CONFIRMATION",
        "eas_95pct_cs_variants_reviewed": len(VARIANTS),
        "immunexut_variant_gene_rows": len(evidences),
        "immunexut_unique_genes": len({r["gene"] for r in evidences}),
        "immunexut_gene_distribution": dict(Counter(r["gene"] for r in evidences)),
        "gene_evidence_records": by_gene,
        "jctf_release_selection": "P_LT_0_05_OR_PIP_GT_0_001",
        "immunexut_release_selection": "BROWSER_FDR_LT_0_05_NOT_FULL_NOMINAL",
        "primary_mechanistic_competitors": ["ALDH2", "BRAP", "RPH3A", "ACAD10", "NAA25", "HECTD4"],
        "warning": "Lower QTL p across different genes/tissues is not causal ranking; 4 high-LD markers are not 4 independent replicates. JCTF other genes may be absent due to release selection.",
        "independent_AIS_replication_verified": False,
        "valid_colocs": 0, "causal_genes_established": 0,
    }


def acquire_source(path, rsid, fetch, overwrite):
    if fetch:
        if path.exists() and not overwrite:
            raise FileExistsError(f"Existing source: {path}")
        req = urllib.request.Request(URL_ROOT + rsid,
                                     headers={"User-Agent": "Mozilla/5.0 (academic G0022 provenance audit)"})
        with urllib.request.urlopen(req, timeout=22) as r:
            body = r.read(200001)
        if len(body) > 200000 or len(body) < 100:
            raise ValueError("Unexpected HTML source length")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
    if not path.is_file():
        raise FileNotFoundError(f"Source HTML missing: {path}")
    return path.read_bytes()


def build(args):
    all_rows, source_manifest = [], []
    for rsid in VARIANTS:
        path = args.html_dir / f"{rsid}.html"
        source = acquire_source(path, rsid, args.fetch, args.overwrite)
        rows = parse_variant_page(source, rsid)
        all_rows.extend(rows)
        source_manifest.append({
            "rsid":rsid, "html_source_url":URL_ROOT+rsid, "html_path":str(path.resolve()),
            "sha256":hashlib.sha256(source).hexdigest(),
            "parsed_significant_rows":len(rows)
        })
    jctf = load_jctf(args.jctf)
    overview, summary = build_summary(all_rows, jctf)
    summary["sources"] = {
        "immunexut":source_manifest,
        "jctf_filtered_variant_level_tsv":str(args.jctf.resolve()),
        "jctf_source_sha256":hashlib.sha256(args.jctf.read_bytes()).hexdigest(),
        "source_html_entries_reused_if_existing":not args.fetch,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    for filename, rows in (
        ("G0022_IMMUNEXUT_4CS_ALL_GENE_VARIANT_ASSOCIATIONS.tsv", all_rows),
        ("G0022_IMMUNEXUT_JCTF_COMPETING_GENE_MATRIX.tsv", overview)):
        with (args.outdir / filename).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    (args.outdir / "G0022_IMMUNEXUT_JCTF_COMPETING_GENE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False)+"\n")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--html-dir", type=Path, required=True)
    p.add_argument("--jctf", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--fetch", action="store_true",
                   help="Fetch exactly four small public variant browser HTML records")
    p.add_argument("--overwrite", action="store_true",
                   help="Explicitly replace source snapshots; absent means immutable")
    print(json.dumps(build(p.parse_args()), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
