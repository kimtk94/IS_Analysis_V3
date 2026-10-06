#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


ANCHOR = Path(
    "/srv/is-analysis/results/ckd/stage4_koges/"
    "controlled_validation/KOGES_CONTROLLED_ANCHOR_SPEC.tsv"
)

OUT = Path(
    "/srv/is-analysis/results/ckd/stage4_koges/"
    "controlled_validation"
)


def detect_genotype(prefix: Path):

    candidates = {
        "PLINK1": [
            Path(str(prefix) + ".bed"),
            Path(str(prefix) + ".bim"),
            Path(str(prefix) + ".fam"),
        ],

        "PLINK2": [
            Path(str(prefix) + ".pgen"),
            Path(str(prefix) + ".pvar"),
            Path(str(prefix) + ".psam"),
        ],

        "BGEN": [
            Path(str(prefix) + ".bgen"),
            Path(str(prefix) + ".sample"),
        ],
    }

    for fmt, files in candidates.items():

        if all(x.exists() for x in files):
            return fmt, files

    return None, []


def read_anchors():

    with ANCHOR.open() as f:
        rows = list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )

    if len(rows) != 9:
        raise RuntimeError(
            f"Expected 9 anchors, got {len(rows)}"
        )

    return rows


def inspect_variant_index(fmt, files, anchors):

    observed = {}

    if fmt == "PLINK1":

        bim = [
            x for x in files
            if x.suffix == ".bim"
        ][0]

        with bim.open() as f:

            for line in f:

                x = line.split()

                if len(x) < 6:
                    continue

                chrom, vid, cm, pos, a1, a2 = x[:6]

                observed[vid] = {
                    "chrom": chrom,
                    "pos": pos,
                    "a1": a1.upper(),
                    "a2": a2.upper(),
                }

    elif fmt == "PLINK2":

        pvar = [
            x for x in files
            if x.suffix == ".pvar"
        ][0]

        with pvar.open() as f:

            header = None

            for line in f:

                if line.startswith("##"):
                    continue

                if line.startswith("#"):

                    header = (
                        line.rstrip()
                        .lstrip("#")
                        .split("\t")
                    )

                    continue

                x = line.rstrip().split("\t")

                if header and len(x) == len(header):
                    d = dict(zip(header, x))
                else:
                    if len(x) < 5:
                        continue

                    d = {
                        "CHROM": x[0],
                        "POS": x[1],
                        "ID": x[2],
                        "REF": x[3],
                        "ALT": x[4],
                    }

                vid = d["ID"]

                observed[vid] = {
                    "chrom": d["CHROM"],
                    "pos": d["POS"],
                    "ref": d["REF"].upper(),
                    "alt": d["ALT"].upper(),
                }

    # BGEN requires plink2 index/extraction later.
    results = []

    for r in anchors:

        rsid = r["rsid"]

        results.append({
            "gene_symbol": r["gene_symbol"],
            "target_rsid": rsid,
            "target_chr": r["chrom_hg19"],
            "target_pos": r["pos_hg19"],
            "target_ref": r["ref_allele"],
            "target_effect_allele":
                r["effect_allele_alt"],

            "direct_rsid_found":
                int(rsid in observed),

            "observed_variant":
                json.dumps(
                    observed.get(rsid, {}),
                    sort_keys=True,
                ),

            "proxy_required":
                int(rsid not in observed),
        })

    return results


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genotype-prefix",
        required=True,
        help=(
            "PLINK1/PLINK2/BGEN prefix "
            "without extension"
        ),
    )

    parser.add_argument(
        "--phenotype",
        required=True,
        help="Controlled KoGES phenotype file",
    )

    args = parser.parse_args()

    prefix = Path(args.genotype_prefix)
    phenotype = Path(args.phenotype)

    anchors = read_anchors()

    fmt, files = detect_genotype(prefix)

    state = {
        "genotype_prefix": str(prefix),
        "genotype_format": fmt,
        "genotype_files": [
            str(x) for x in files
        ],
        "phenotype": str(phenotype),
        "phenotype_exists":
            phenotype.exists(),
    }

    if fmt is None:

        state["status"] = (
            "FAIL_GENOTYPE_NOT_RECOGNIZED"
        )

        (OUT / "CONTROLLED_INPUT_QC.json").write_text(
            json.dumps(state, indent=2)
        )

        raise SystemExit(
            "No complete PLINK1/PLINK2/BGEN dataset found."
        )

    if not phenotype.exists():

        state["status"] = (
            "FAIL_PHENOTYPE_NOT_FOUND"
        )

        (OUT / "CONTROLLED_INPUT_QC.json").write_text(
            json.dumps(state, indent=2)
        )

        raise SystemExit(
            f"Phenotype not found: {phenotype}"
        )

    if fmt in {"PLINK1", "PLINK2"}:

        results = inspect_variant_index(
            fmt,
            files,
            anchors,
        )

        out = (
            OUT /
            "KOGES_DIRECT_ANCHOR_AVAILABILITY.tsv"
        )

        fields = [
            "gene_symbol",
            "target_rsid",
            "target_chr",
            "target_pos",
            "target_ref",
            "target_effect_allele",
            "direct_rsid_found",
            "observed_variant",
            "proxy_required",
        ]

        with out.open("w", newline="") as f:

            w = csv.DictWriter(
                f,
                fieldnames=fields,
                delimiter="\t",
            )

            w.writeheader()
            w.writerows(results)

        n_direct = sum(
            int(x["direct_rsid_found"])
            for x in results
        )

        state["direct_anchors"] = n_direct
        state["proxy_required"] = 9 - n_direct

    else:

        state["direct_anchors"] = None
        state["proxy_required"] = None
        state["note"] = (
            "BGEN detected; variant availability "
            "must be queried with plink2."
        )

    state["status"] = "INPUTS_RECOGNIZED"

    (OUT / "CONTROLLED_INPUT_QC.json").write_text(
        json.dumps(state, indent=2)
    )

    print("=" * 64)
    print("KoGES CONTROLLED INPUT VALIDATION")
    print("=" * 64)

    print("format =", fmt)
    print(
        "phenotype_exists =",
        phenotype.exists(),
    )

    print(
        "direct_anchors =",
        state["direct_anchors"],
    )

    print(
        "proxy_required =",
        state["proxy_required"],
    )

    print("status =", state["status"])


if __name__ == "__main__":
    main()
