"""Validate public IS nominal-QTL readiness TSV without claiming causal inference."""
import argparse
import csv
import json
import re
from pathlib import Path

REQUIRED = (
    "gene_id_stable", "canonical_gene_symbol", "region_id",
    "archived_646_coloc_assays", "n_new_QTL_tissues_with_any_nominal_rows",
    "source_nominal_QTL_ids", "n_new_nominal_QTL_gene_variant_rows",
    "min_source_nominal_QTL_p", "new_nominal_QTL_gene_source_status",
    "GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified",
    "coloc_recomputed", "ischemic_stroke_causal_gene",
)
BOOL = ("GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified", "coloc_recomputed")
COUNTS = ("archived_646_coloc_assays", "n_new_QTL_tissues_with_any_nominal_rows",
          "n_new_nominal_QTL_gene_variant_rows")

def validate(path, expected_genes=2225):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Missing or duplicated header")
        missing = set(REQUIRED) - set(reader.fieldnames)
        if missing:
            raise ValueError("Missing columns: " + ", ".join(sorted(missing)))
        seen = set()
        status_counts = {}
        for lineno, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed row {lineno}")
            gene = row["gene_id_stable"]
            if not re.fullmatch(r"ENSG\d{11}", gene):
                raise ValueError(f"Invalid Ensembl ID at row {lineno}")
            if gene in seen:
                raise ValueError(f"Duplicate gene at row {lineno}: {gene}")
            seen.add(gene)
            for key in BOOL:
                if row[key] not in ("True", "False"):
                    raise ValueError(f"Invalid boolean {key} at row {lineno}")
            for key in COUNTS:
                if not row[key].isdigit():
                    raise ValueError(f"Invalid count {key} at row {lineno}")
            count = int(row["n_new_nominal_QTL_gene_variant_rows"])
            p = row["min_source_nominal_QTL_p"]
            if count and not p:
                raise ValueError(f"Missing nominal p-value at row {lineno}")
            if p and not 0 < float(p) <= 1:
                raise ValueError(f"Invalid p-value at row {lineno}")
            if row["ischemic_stroke_causal_gene"] != "NOT_ESTABLISHED":
                raise ValueError(f"Unreviewed causal status at row {lineno}")
            status = row["new_nominal_QTL_gene_source_status"]
            status_counts[status] = status_counts.get(status, 0) + 1
    if len(seen) != expected_genes:
        raise ValueError(f"Expected {expected_genes} genes; found {len(seen)}")
    return {"status": "PASS", "unique_genes": len(seen), "source_status_counts": status_counts,
            "causal_gene_claims": 0}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", type=Path)
    parser.add_argument("--expected-genes", type=int, default=2225)
    args = parser.parse_args()
    print(json.dumps(validate(args.tsv, args.expected_genes), sort_keys=True))
