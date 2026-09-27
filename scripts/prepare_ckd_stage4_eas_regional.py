#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace",
        )
    return path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    )


def read_tsv(path: Path):
    with open_text(path) as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t",
            )
        )


def norm_chr(x):
    s = str(x).strip()

    if s.lower().startswith("chr"):
        s = s[3:]

    return s.upper()


def norm_a(x):
    return str(x).strip().upper()


def complement(a):
    return {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }.get(a)


def is_snv(a, b):
    return (
        len(a) == 1
        and len(b) == 1
        and a in "ACGT"
        and b in "ACGT"
        and a != b
    )


def is_palindromic(a, b):
    return {a, b} in (
        {"A", "T"},
        {"C", "G"},
    )


def fnum(x):
    try:
        v = float(str(x).strip())
        return v if math.isfinite(v) else None
    except Exception:
        return None


def harmonize(
    pqtl_a0,
    pqtl_a1,
    outcome_ea,
    outcome_oa,
):
    """
    UKB-PPP ALLELE1 = exposure effect allele.

    Return:
        sign to apply to outcome beta
        orientation label
    """

    a0 = norm_a(pqtl_a0)
    a1 = norm_a(pqtl_a1)

    ea = norm_a(outcome_ea)
    oa = norm_a(outcome_oa)

    if not (
        is_snv(a0, a1)
        and is_snv(ea, oa)
    ):
        return None

    # direct:
    # outcome effect allele == pQTL effect allele
    if ea == a1 and oa == a0:
        return 1.0, "direct"

    # swapped
    if ea == a0 and oa == a1:
        return -1.0, "swapped"

    # Palindromic variants cannot safely resolve strand
    # without ancestry-matched allele frequency.
    if (
        is_palindromic(a0, a1)
        or is_palindromic(ea, oa)
    ):
        return None

    cea = complement(ea)
    coa = complement(oa)

    if cea == a1 and coa == a0:
        return 1.0, "strand_direct"

    if cea == a0 and coa == a1:
        return -1.0, "strand_swapped"

    return None


def read_outcome(
    path: Path,
    wanted_positions,
):
    out = defaultdict(list)

    parsed = 0
    malformed = 0

    with open_text(path) as fh:

        header = re.split(
            r"\s+",
            fh.readline().strip(),
        )

        idx = {
            name: i
            for i, name in enumerate(header)
        }

        required = {
            "Allele1",
            "Allele2",
            "Effect",
            "StdErr",
            "P",
            "CHR",
            "BP",
            "SNP",
        }

        missing = required - set(header)

        if missing:
            raise SystemExit(
                f"outcome missing columns: {sorted(missing)}; "
                f"header={header}"
            )

        for line in fh:

            if not line.strip():
                continue

            vals = re.split(
                r"\s+",
                line.strip(),
            )

            if len(vals) != len(header):
                malformed += 1
                continue

            parsed += 1

            chrom = norm_chr(
                vals[idx["CHR"]]
            )

            try:
                pos = int(
                    float(vals[idx["BP"]])
                )
            except Exception:
                continue

            key = (chrom, pos)

            if key not in wanted_positions:
                continue

            ea = norm_a(
                vals[idx["Allele1"]]
            )

            oa = norm_a(
                vals[idx["Allele2"]]
            )

            beta = fnum(
                vals[idx["Effect"]]
            )

            se = fnum(
                vals[idx["StdErr"]]
            )

            p = fnum(
                vals[idx["P"]]
            )

            if (
                beta is None
                or se is None
                or se <= 0
            ):
                continue

            out[key].append({
                "rsid":
                    vals[idx["SNP"]].strip(),

                "effect_allele":
                    ea,

                "other_allele":
                    oa,

                "beta":
                    beta,

                "se":
                    se,

                "p":
                    p,
            })

    return out, parsed, malformed


def choose_outcome(
    pqtl,
    candidates,
):
    hits = []

    for y in candidates:

        h = harmonize(
            pqtl["allele0_pqtl_hg37"],
            pqtl["allele1_pqtl_hg37"],
            y["effect_allele"],
            y["other_allele"],
        )

        if h is None:
            continue

        sign, orient = h

        rsid_exact = int(
            y["rsid"]
            and y["rsid"]
            == pqtl.get("rsid", "")
        )

        p = y["p"]
        score_p = (
            p if p is not None else 1.0
        )

        hits.append(
            (
                -rsid_exact,
                score_p,
                y,
                sign,
                orient,
            )
        )

    if not hits:
        return None

    hits.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    return hits[0][2:]


def write_tsv_gz(
    path: Path,
    rows,
    fields,
):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with gzip.open(
        path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        w.writeheader()
        w.writerows(rows)


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--anchor",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--pqtl-locus-dir",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--outcome",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--phenotype",
        required=True,
    )

    ap.add_argument(
        "--outcome-n",
        type=int,
        required=True,
    )

    ap.add_argument(
        "--output-root",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--min-shared-snps",
        type=int,
        default=50,
    )

    args = ap.parse_args()

    anchors = read_tsv(args.anchor)

    genes = [
        x["gene_symbol"].upper()
        for x in anchors
    ]

    if len(genes) != 9:
        raise SystemExit(
            f"expected 9 genes, got {len(genes)}"
        )

    pqtl_by_gene = {}
    wanted_positions = set()

    for gene in genes:

        path = (
            args.pqtl_locus_dir
            / f"{gene}.tsv.gz"
        )

        if not path.is_file():
            raise SystemExit(
                f"missing pQTL locus: {path}"
            )

        rows = read_tsv(path)

        if not rows:
            raise SystemExit(
                f"{gene}: empty pQTL locus"
            )

        pqtl_by_gene[gene] = rows

        for r in rows:
            wanted_positions.add(
                (
                    norm_chr(r["chr37"]),
                    int(float(r["pos37"])),
                )
            )

    outcome, parsed, malformed = (
        read_outcome(
            args.outcome,
            wanted_positions,
        )
    )

    fields = [
        "protein_id",
        "gene_symbol",

        "snp",
        "chr37",
        "pos37",

        "allele0_pqtl",
        "allele1_pqtl",

        "a1freq_pqtl",
        "maf_pqtl",

        "n_pqtl",
        "beta_pqtl",
        "se_pqtl",
        "p_pqtl",

        "effect_allele_outcome_original",
        "other_allele_outcome_original",

        "n_outcome",

        "beta_outcome_original",
        "beta_outcome_aligned",

        "se_outcome",
        "p_outcome",

        "harmonization",

        "variant_id_grch38",
        "chr38",
        "pos38",
    ]

    qc = []

    args.output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    for gene in genes:

        matched = []

        counts = defaultdict(int)

        for p in pqtl_by_gene[gene]:

            key = (
                norm_chr(p["chr37"]),
                int(float(p["pos37"])),
            )

            choices = outcome.get(
                key,
                [],
            )

            if not choices:
                counts[
                    "no_outcome_position"
                ] += 1
                continue

            m = choose_outcome(
                p,
                choices,
            )

            if m is None:
                counts[
                    "allele_or_palindrome_unresolved"
                ] += 1
                continue

            y, sign, orient = m

            rsid = (
                y["rsid"]
                if y["rsid"].startswith("rs")
                else (
                    f"{key[0]}:"
                    f"{key[1]}:"
                    f"{p['allele0_pqtl_hg37']}:"
                    f"{p['allele1_pqtl_hg37']}"
                )
            )

            matched.append({
                "protein_id":
                    p["protein_id"],

                "gene_symbol":
                    gene,

                "snp":
                    rsid,

                "chr37":
                    key[0],

                "pos37":
                    key[1],

                "allele0_pqtl":
                    p["allele0_pqtl_hg37"],

                "allele1_pqtl":
                    p["allele1_pqtl_hg37"],

                "a1freq_pqtl":
                    p["a1freq_pqtl"],

                "maf_pqtl":
                    p["maf_pqtl"],

                "n_pqtl":
                    p["n_pqtl"],

                "beta_pqtl":
                    p["beta_pqtl"],

                "se_pqtl":
                    p["se_pqtl"],

                "p_pqtl":
                    p["p_pqtl"],

                "effect_allele_outcome_original":
                    y["effect_allele"],

                "other_allele_outcome_original":
                    y["other_allele"],

                "n_outcome":
                    args.outcome_n,

                "beta_outcome_original":
                    y["beta"],

                "beta_outcome_aligned":
                    y["beta"] * sign,

                "se_outcome":
                    y["se"],

                "p_outcome":
                    y["p"],

                "harmonization":
                    orient,

                "variant_id_grch38":
                    p["variant_id_grch38"],

                "chr38":
                    p["chr38"],

                "pos38":
                    p["pos38"],
            })

            counts["matched"] += 1
            counts[
                f"harmonization_{orient}"
            ] += 1

        # Unique SNP ID.
        best = {}

        for r in matched:

            snp = r["snp"]

            pval = (
                fnum(r["p_pqtl"])
                or 1.0
            )

            if (
                snp not in best
                or pval
                < (
                    fnum(
                        best[snp]["p_pqtl"]
                    )
                    or 1.0
                )
            ):
                best[snp] = r

        matched = list(
            best.values()
        )

        matched.sort(
            key=lambda r: (
                int(r["pos37"]),
                r["snp"],
            )
        )

        out = (
            args.output_root
            / "coloc_input"
            / f"{gene}.tsv.gz"
        )

        write_tsv_gz(
            out,
            matched,
            fields,
        )

        qc_row = {
            "gene_symbol":
                gene,

            "pqtl_locus_snvs":
                len(pqtl_by_gene[gene]),

            "shared_unique_snps":
                len(matched),

            "match_fraction":
                (
                    len(matched)
                    / len(pqtl_by_gene[gene])
                ),

            "min_shared_snps_pass":
                int(
                    len(matched)
                    >= args.min_shared_snps
                ),

            **dict(counts),
        }

        qc.append(qc_row)

    qc_fields = sorted(
        {
            k
            for r in qc
            for k in r
        }
    )

    qc_path = (
        args.output_root
        / "REGIONAL_MATCH_QC.tsv"
    )

    with qc_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=qc_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        w.writeheader()
        w.writerows(qc)

    summary = {
        "stage":
            "CKD EAS regional cross-ancestry preparation",

        "phenotype":
            args.phenotype,

        "outcome_n":
            args.outcome_n,

        "outcome_rows_parsed":
            parsed,

        "malformed_outcome_rows":
            malformed,

        "candidate_count":
            len(genes),

        "minimum_shared_snps":
            args.min_shared_snps,

        "matching":
            (
                "GRCh37 chromosome/position "
                "+ allele pair; outcome beta "
                "aligned to UKB-PPP ALLELE1"
            ),

        "palindromic_policy":
            (
                "strand-unresolved palindromic "
                "matches excluded because EAS "
                "summary file has no allele frequency"
            ),

        "qc":
            qc,
    }

    (
        args.output_root
        / "REGIONAL_PREP_SUMMARY.json"
    ).write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        "GENE\tPQTL\tSHARED\tFRACTION"
    )

    for r in qc:
        print(
            r["gene_symbol"],
            r["pqtl_locus_snvs"],
            r["shared_unique_snps"],
            f'{r["match_fraction"]:.4f}',
            sep="\t",
        )

    failed = [
        x
        for x in qc
        if not x["min_shared_snps_pass"]
    ]

    if failed:
        print(
            "WARNING: candidates below "
            "minimum shared SNP threshold:"
        )

        for r in failed:
            print(
                r["gene_symbol"],
                r["shared_unique_snps"],
            )

    print(
        "CKD_EAS_REGIONAL_PREP_PASS"
    )


if __name__ == "__main__":
    main()
