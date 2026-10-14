#!/usr/bin/env python3
"""Read-only discovery of possible IS GWAS/QTL inputs; no raw data copied."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

EXTENSIONS = (".tsv", ".tsv.gz", ".txt", ".txt.gz", ".csv", ".csv.gz",
              ".parquet", ".vcf", ".vcf.gz", ".bim", ".pvar", ".pvar.zst", ".rds")
TOKENS = ("gwas", "qtl", "eqtl", "sqtl", "pqtl", "nominal", "bbj",
          "gigastroke", "g0025", "g0015", "g0018", "g0023", "g0007")
HEADER_HINTS = ("chr", "chrom", "pos", "position", "ref", "alt", "beta",
                "effect_allele", "other_allele", "pvalue", "variant_id")

def discover(roots, out, max_files=100000):
    found = []
    scanned = 0
    for root in roots:
        root = Path(root).resolve()
        if not root.is_dir():
            continue
        for path in root.rglob("*"):
            if path.is_symlink() or not path.is_file():
                continue
            scanned += 1
            if scanned > max_files:
                raise ValueError("File scan cap exceeded; narrow roots")
            lower = path.name.lower()
            if not lower.endswith(EXTENSIONS):
                continue
            relative = str(path.relative_to(root))
            context = relative.lower()
            if not any(t in context for t in TOKENS):
                continue
            stat = path.stat()
            found.append({"root": str(root), "relative_path": relative,
                          "bytes": stat.st_size, "suffix": "".join(path.suffixes),
                          "matched_tokens": ";".join(t for t in TOKENS if t in context),
                          "header_checked": "NO",
                          "genome_build": "UNKNOWN",
                          "ancestry": "UNKNOWN",
                          "source_verified": "NO",
                          "eligible_for_coloc": "NO"})
    found.sort(key=lambda x: (x["root"], x["relative_path"]))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fieldnames = ["root", "relative_path", "bytes", "suffix", "matched_tokens",
                  "header_checked", "genome_build", "ancestry",
                  "source_verified", "eligible_for_coloc"]
    with (out / "IS_POSSIBLE_SNP_INPUT_INVENTORY.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(found)
    manifest = {"status": "DISCOVERY_ONLY_NOT_VALIDATED",
                "scanned_files": scanned, "possible_inputs": len(found),
                "roots": [str(Path(p).resolve()) for p in roots],
                "warning": "File names and sizes do not verify ancestry, build, columns or scientific eligibility."}
    (out / "IS_POSSIBLE_SNP_INPUT_INVENTORY.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", action="append", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--max-files", type=int, default=100000)
    a = parser.parse_args()
    print(json.dumps(discover(a.root, a.out, a.max_files), indent=2))
