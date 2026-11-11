#!/usr/bin/env python3
"""GTEx v8 median-TPM tissue specificity pilot for prespecified IS candidates.

Uses a complete public aggregate expression atlas (NOT eQTL, differential
expression, or causal locus-to-gene evidence). Prints all missing genes and
records data SHA256, release and tissue count. No selection by association
strength. Source is public GTEx GRCh38 / GENCODE 26.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen

SOURCE = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
          "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
ALTERNATE = ("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
             "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
MAX_BYTES = 30 * 1024 * 1024
CANDIDATES = ("ABO", "ALDH2", "CD40", "F11", "FGF5", "FURIN", "GRK5",
              "IL6RA", "KIAA0319", "MMP12", "PROCR", "SCARA5", "SWAP70",
              "TFPI", "TMPRSS5")
SYMBOL_ALIASES = {"IL6RA": "IL6R"}  # official HGNC gene symbol versus legacy alias
TISSUES = {
    "Artery_Aorta": "Artery - Aorta",
    "Artery_Coronary": "Artery - Coronary",
    "Artery_Tibial": "Artery - Tibial",
    "Heart_Atrial": "Heart - Atrial Appendage",
    "Heart_Ventricle": "Heart - Left Ventricle",
    "Whole_Blood": "Whole Blood",
    "Brain_Cortex": "Brain - Cortex",
    "Brain_Frontal_BA9": "Brain - Frontal Cortex (BA9)",
    "Brain_Cerebellum": "Brain - Cerebellum",
}
GROUPS = {
    "arterial": ("Artery_Aorta", "Artery_Coronary", "Artery_Tibial"),
    "cardiac": ("Heart_Atrial", "Heart_Ventricle"),
    "brain": ("Brain_Cortex", "Brain_Frontal_BA9", "Brain_Cerebellum"),
    "blood": ("Whole_Blood",),
}


def tau(values):
    if not values:
        return None
    vmax = max(values)
    if vmax <= 0:
        return None
    return sum(1 - x / vmax for x in values) / max(1, len(values) - 1)


def fetch(url, dest):
    req = Request(url, headers={"User-Agent": "IS_Analysis_V3_GTEx_Research/1.0"})
    digest = hashlib.sha256()
    count = 0
    with urlopen(req, timeout=40) as response:
        final = urlparse(response.geturl())
        if (final.scheme != "https" or final.hostname != "storage.googleapis.com"
                or response.status != 200):
            raise RuntimeError("Unexpected GTEx download origin or HTTP response")
        declared = response.headers.get("Content-Length", "")
        with open(dest, "wb") as out:
            while True:
                part = response.read(1 << 20)
                if not part:
                    break
                count += len(part)
                if count > MAX_BYTES:
                    raise RuntimeError("GTEx source exceeds safety cap")
                digest.update(part)
                out.write(part)
        if declared.isdigit() and count != int(declared):
            raise RuntimeError("Incomplete GTEx transfer")
    return {"source_url": url, "compressed_bytes": count,
            "local_sha256_not_publisher_verified": digest.hexdigest()}


def parse_gct(path, targets=CANDIDATES):
    matches = {name: [] for name in targets}
    all_gene_names = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        version = stream.readline().strip()
        shape = stream.readline().strip()
        if not version.startswith("#1."):
            raise ValueError("Unexpected GCT version")
        parts = shape.split()
        if len(parts) < 2 or not all(x.isdigit() for x in parts[:2]):
            raise ValueError("Missing GCT matrix dimensions")
        reader = csv.DictReader(stream, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "Name" not in fields or "Description" not in fields:
            raise ValueError("GTEx GCT lacks gene ID and symbol columns")
        tissue_cols = [x for x in fields if x not in ("Name", "Description")]
        if len(tissue_cols) < 40:
            raise ValueError("Unexpectedly few GTEx tissues")
        row_count = 0
        for row in reader:
            row_count += 1
            symbol = (row.get("Description") or "").strip()
            name = next((k for k in matches if SYMBOL_ALIASES.get(k, k) == symbol), None)
            if name is None:
                continue
            values = []
            for tissue in tissue_cols:
                value = float(row[tissue])
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f"Invalid TPM in {symbol}")
                values.append(value)
            matches[name].append((row["Name"], values))
        if int(parts[0]) != row_count or int(parts[1]) != len(tissue_cols):
            raise ValueError("GCT shape/count mismatch")
    return tissue_cols, matches, row_count


def summarize(tissues, matches):
    missing_panel = [v for v in TISSUES.values() if v not in tissues]
    if missing_panel:
        raise ValueError("Missing required GTEx tissues: " + ",".join(missing_panel))
    ix = {x: i for i, x in enumerate(tissues)}
    rows, all_tissues = [], []
    for gene in CANDIDATES:
        records = matches[gene]
        if len(records) != 1:
            rows.append({"gene_symbol": gene, "status":
                        "NOT_FOUND" if not records else "DUPLICATE_GENE_SYMBOL"})
            continue
        ensg, values = records[0]
        order = sorted(range(len(values)), key=lambda i: values[i], reverse=True)
        top_idx = order[0]
        entry = {
            "gene_symbol": gene, "mapped_gene_symbol": SYMBOL_ALIASES.get(gene, gene), "gene_id_gencode_versioned": ensg,
            "status": "EXPRESSION_AVAILABLE_NOT_EQTL", "top_tissue": tissues[top_idx],
            "top_median_tpm": round(values[top_idx], 4),
            "tau_all_tissues": round(tau(values), 4) if tau(values) is not None else None,
            "n_tissues": len(tissues),
        }
        for alias, tissue in TISSUES.items():
            median = values[ix[tissue]]
            entry[alias + "_tpm"] = round(median, 4)
            entry[alias + "_rank_all_tissues"] = order.index(ix[tissue]) + 1
        for group, aliases in GROUPS.items():
            score = max(values[ix[TISSUES[alias]]] for alias in aliases)
            entry[group + "_max_tpm"] = round(score, 4)
            # Purely descriptive ratio; do not treat as formal gene-property P value.
            entry[group + "_fraction_of_gene_tissue_max"] = (
                round(score / max(values), 4) if max(values) > 0 else None
            )
        rows.append(entry)
        for tissue, val in zip(tissues, values):
            all_tissues.append({
                "gene_symbol": gene, "gene_id_gencode_versioned": ensg,
                "tissue": tissue, "median_tpm": round(val, 6),
            })
    return rows, all_tissues


def write_tsv(path, rows):
    if not rows:
        return
    fields = sorted(set().union(*(set(x) for x in rows)))
    with open(path, "w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="artifacts/is_gtex_specificity")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    evidence = {"status": "SOURCE_UNAVAILABLE", "data_type": "GTEx_v8_bulk_gene_median_TPM",
                "tissue_specificity_not_GWAS_eQTL_colocalization": True,
                "candidate_genes_prespecified": list(CANDIDATES),
                "genome_build": "GRCh38",
                "genecode_version": "v26",
                "requires_liftover_for_grch37_gwas": True,
                "source": SOURCE}
    try:
        with tempfile.TemporaryDirectory(prefix="gtex_v8_") as d:
            path = Path(d) / "gtex_median_tpm.gct.gz"
            try:
                evidence.update(fetch(SOURCE, path))
            except (OSError, RuntimeError) as err:
                evidence["primary_attempt_error"] = f"{type(err).__name__}: {str(err)[:140]}"
                evidence.update(fetch(ALTERNATE, path))
            tissues, matches, total_genes = parse_gct(path)
            rows, long_rows = summarize(tissues, matches)
        write_tsv(out / "target_gene_tissue_specificity.tsv", rows)
        write_tsv(out / "target_gene_all_tissues.tsv", long_rows)
        evidence.update({"status": "MEDIAN_TPM_QC_PASS",
                         "tested_gene_rows": total_genes, "tissue_count": len(tissues),
                         "candidate_genes_matched": sum(x["status"] == "EXPRESSION_AVAILABLE_NOT_EQTL" for x in rows),
                         "missing_genes": [x["gene_symbol"] for x in rows if x["status"] != "EXPRESSION_AVAILABLE_NOT_EQTL"]})
        md = ["# GTEx v8 median tissue expression: prespecified IS candidate pilot", "",
              "**This is normal bulk tissue expression, NOT an eQTL / stroke case-control experiment.**",
              "", "| Gene | Top tissue | Top TPM | Tau | Artery max TPM | Heart max TPM | Blood TPM | Brain max TPM |",
              "|---|---|---:|---:|---:|---:|---:|---:|"]
        for x in rows:
            if x["status"] != "EXPRESSION_AVAILABLE_NOT_EQTL":
                md.append(f"| {x['gene_symbol']} | {x['status']} | - | - | - | - | - | - |")
            else:
                md.append("| {gene_symbol} | {top_tissue} | {top_median_tpm} | {tau_all_tissues} | {arterial_max_tpm} | {cardiac_max_tpm} | {Whole_Blood_tpm} | {brain_max_tpm} |".format(**x))
        md += ["", "Tau is a descriptive tissue-specificity index across all GTEx v8 bulk tissue median TPM values.",
               "High expression or Tau does not imply stroke genetic association, colocalization or causal regulation.",
               "GTEx v8 is GRCh38; the verified GIGASTROKE sources are GRCh37; liftover and allele checking remain HOLD."]
        summary = "\n".join(md) + "\n"
        (out / "summary.md").write_text(summary, encoding="utf-8")
        print(summary)
        if os.getenv("GITHUB_STEP_SUMMARY"):
            with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
                stream.write(summary)
    except (OSError, RuntimeError, ValueError, EOFError) as exc:
        evidence["error"] = f"{type(exc).__name__}: {str(exc)[:250]}"
    (out / "source_status.json").write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("SOURCE_STATUS:", evidence["status"])
    # Status written independently; cannot report PASS if inaccessible.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
