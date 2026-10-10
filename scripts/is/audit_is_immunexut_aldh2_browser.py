#!/usr/bin/env python3
"""Audit rs671 + 3 high-LD AIS-CS ALDH2 eQTL entries in Japanese ImmuNexUT.

This is a significant-result browser only, not full nominal QTL or coloc.
No allele orientation is assumed and no beta is compared with AIS GWAS.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import urllib.request
from pathlib import Path

SOURCE_URL = "https://www.immunexut.org/eqtlGenes?gene_symbol=ALDH2"
FOUR_CS = (
    ("12:112168009:G:A", "rs11066015", 111730205),
    ("12:112241766:G:A", "rs671", 111803962),
    ("12:112468206:C:T", "rs11066132", 112030402),
    ("12:112736118:A:G", "rs77768175", 112298314),
)


def parse_page(source):
    if isinstance(source, bytes):
        source = source.decode("utf-8", "replace")
    matches = []
    for item in re.finditer(r'\{"gene_symbol":"ALDH2"[^{}]{20,800}\}', source):
        try:
            record = json.loads(item.group(0))
        except ValueError:
            continue
        if isinstance(record, dict) and "snp_id" in record and "eqtl_pval1" in record:
            matches.append(record)
    if not matches:
        raise ValueError("No parsable ImmuNexUT ALDH2 gene eQTL records")
    if len(matches) > 1000:
        raise ValueError("Unexpected >1000 records; source contract changed")
    observed = []
    for gwas_id, rsid, b38pos in FOUR_CS:
        candidates = [x for x in matches if x.get("snp_id") == rsid and
                      x.get("position") == b38pos and x.get("chromosome") == "chr12" and
                      x.get("cell_type") == "Neu"]
        if len(candidates) != 1:
            raise ValueError(f"Missing/duplicate Neu ALDH2 CS marker {rsid}: {len(candidates)}")
        rec = candidates[0]
        p = float(rec["eqtl_pval1"])
        beta = float(rec["eqtl_effect_beta1"])
        if not (math.isfinite(p) and math.isfinite(beta) and 0 < p <= 1):
            raise ValueError("Invalid eQTL effect or p value")
        observed.append({
            "gwas_cs_variant_grch37": gwas_id, "qtl_rsid": rsid,
            "qtl_grch38_chr": "12", "qtl_grch38_pos": b38pos,
            "molecular_gene": "ALDH2", "cell_type_code": "Neu",
            "cell_type_label": "Neutrophils",
            "qtl_browser_nominal_p": p, "qtl_browser_beta_unharmonized": beta,
            "gwas_qtl_effect_alleles_harmonized": False,
            "source_effect_allele_extracted": False,
            "same_causal_variant_established": False,
            "qtl_browser_publication_gate": "SUPPORTED_MARGINAL_EQTL_ONLY",
        })
    summary = {
        "audit": "IS_G0022_IMMUNEXUT_ALDH2_NEUTROPHIL_BROWSER_V1",
        "source": SOURCE_URL,
        "source_genome_build": "GRCh38",
        "gwas_genome_build": "GRCh37",
        "browser_records_extracted": len(matches),
        "browser_displays_top_n_if_many": 1000,
        "four_gwas_cs_variants_found_as_neutrophil_ALDH2_eQTL": len(observed),
        "qtl_p_values": [x["qtl_browser_nominal_p"] for x in observed],
        "qtl_unharmonized_betas": [x["qtl_browser_beta_unharmonized"] for x in observed],
        "identical_four_p_and_beta": len({(x["qtl_browser_nominal_p"], x["qtl_browser_beta_unharmonized"]) for x in observed}) == 1,
        "gwas_qtl_allele_direction_compared": False,
        "qtl_full_nominal_non_significant_snp_universe_loaded": False,
        "qtl_signed_ld_available": False,
        "new_valid_coloc": 0,
        "validated_causal_gene": False,
        "scientific_state": "FOUR_VARIANT_ASSOCIATION_NEUTROPHIL_SOURCE_SUPPORTED_NOT_COLOCALIZATION",
        "warning": "Distinct source from JCTF, but participants/technical independence not audited. Four correlated variants are not independent functional replications. No RNA association implies disease mediation.",
    }
    return observed, summary


def build(args):
    pagepath = Path(args.source_html)
    if args.fetch:
        if pagepath.exists() and not args.overwrite:
            raise FileExistsError("Existing HTML exists: specify --overwrite for explicit new source version")
        request = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Mozilla/5.0 (G0022 scientific audit)"})
        with urllib.request.urlopen(request, timeout=20) as response:
            content = response.read(1_000_001)
        if len(content) > 1_000_000:
            raise ValueError("Source HTML unexpectedly large")
        pagepath.parent.mkdir(parents=True, exist_ok=True)
        pagepath.write_bytes(content)
    if not pagepath.is_file():
        raise FileNotFoundError("Provide source HTML or opt in to --fetch")
    payload = pagepath.read_bytes()
    rows, summary = parse_page(payload)
    summary["html_sha256"] = hashlib.sha256(payload).hexdigest()
    summary["source_html"] = str(pagepath.resolve())
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "G0022_IMMUNEXUT_ALDH2_FOUR_CS_NEU_EQTL.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (out / "G0022_IMMUNEXUT_ALDH2_NEU_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-html", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--fetch", action="store_true", help="One bounded HTML request (no 43GB archive download)")
    p.add_argument("--overwrite", action="store_true", help="Explicitly replace previous source HTML")
    args = p.parse_args()
    print(json.dumps(build(args), indent=2))


if __name__ == "__main__":
    main()
