#!/usr/bin/env python3
"""Create an exploratory, noncausal region-level QTL replay queue from public TSV."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

READY = "RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC"

def build(source, out, expected_genes=2225):
    raw = Path(source).read_bytes()
    lines = raw.decode("utf-8").splitlines()
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if len(rows) != expected_genes:
        raise ValueError(f"Expected {expected_genes} genes, found {len(rows)}")
    ids = [r["gene_id_stable"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate stable gene IDs")
    groups = defaultdict(list)
    for r in rows:
        if r["new_nominal_QTL_gene_source_status"] != READY:
            continue
        region = r["region_id"]
        if not region:
            raise ValueError("QTL-ready gene missing region")
        if r["coloc_recomputed"] != "False":
            raise ValueError("Unexpected recomputed coloc claim")
        if r["GWAS_and_QTL_alleles_coharmonized"] != "False":
            raise ValueError("Unexpected harmonization claim")
        groups[region].append(r)
    if not groups:
        raise ValueError("No nominal QTL-ready regions")
    result = []
    for region, genes in sorted(groups.items(), key=lambda x: (-len(x[1]), x[0])):
        sources = sorted({source_id for g in genes for source_id in
                          g["source_nominal_QTL_ids"].split(";") if source_id})
        if not sources:
            raise ValueError(f"Missing QTL source IDs: {region}")
        result.append({
            "region_id": region,
            "n_genes_with_nominal_qtl": len(genes),
            "stable_gene_ids": ";".join(sorted(g["gene_id_stable"] for g in genes)),
            "qtl_source_ids": ";".join(sources),
            "allele_harmonization": "NOT_ASSESSED",
            "ancestry_matched_ld": "NOT_ASSESSED",
            "variant_overlap": "NOT_ASSESSED",
            "colocalization": "NOT_ASSESSED",
            "causal_gene_claim": "NOT_ESTABLISHED",
        })
    output = Path(out)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "IS_PUBLIC_NOMINAL_QTL_REGION_REPLAY_QUEUE.tsv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(result)
    manifest = {
        "status": "EXPLORATORY_QUEUE_ONLY",
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "discovery_denominator": len(rows),
        "qtl_ready_genes": sum(len(v) for v in groups.values()),
        "qtl_ready_regions": len(groups),
        "region_gene_counts": {r["region_id"]: r["n_genes_with_nominal_qtl"] for r in result},
        "note": "Nominal QTL source presence is not evidence of colocalization or causality.",
    }
    (output / "IS_PUBLIC_NOMINAL_QTL_REGION_REPLAY_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--expected-genes", type=int, default=2225)
    args = p.parse_args()
    print(json.dumps(build(args.source, args.out, args.expected_genes), indent=2))
