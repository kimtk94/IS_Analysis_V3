#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_tsv(path: Path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def allele(x):
    return (x or "").strip().upper()


def complement(a):
    table = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }
    return table.get(a, a)


def is_palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def harmonize(x_ea, x_oa, y_ea, y_oa):
    x_ea = allele(x_ea)
    x_oa = allele(x_oa)
    y_ea = allele(y_ea)
    y_oa = allele(y_oa)

    if y_ea == x_ea and y_oa == x_oa:
        return 1.0, "direct"

    if y_ea == x_oa and y_oa == x_ea:
        return -1.0, "reverse"

    cy_ea = complement(y_ea)
    cy_oa = complement(y_oa)

    if not is_palindromic(x_ea, x_oa):
        if cy_ea == x_ea and cy_oa == x_oa:
            return 1.0, "strand_direct"

        if cy_ea == x_oa and cy_oa == x_ea:
            return -1.0, "strand_reverse"

    raise ValueError(
        f"cannot harmonize "
        f"exposure={x_ea}/{x_oa} outcome={y_ea}/{y_oa}"
    )


def normal_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def fnum(x):
    if x is None or x == "":
        return float("nan")
    return float(x)


def wald_ratio(beta_x, se_x, beta_y, se_y):
    beta = beta_y / beta_x

    # Conventional first-order Wald-ratio SE.
    se_first = se_y / abs(beta_x)

    # Delta-method version retaining exposure uncertainty.
    se_delta = math.sqrt(
        (se_y ** 2) / (beta_x ** 2)
        +
        (beta_y ** 2 * se_x ** 2) / (beta_x ** 4)
    )

    z = beta / se_first
    p = normal_p(z)

    return {
        "wald_beta": beta,
        "wald_se_first_order": se_first,
        "wald_se_delta": se_delta,
        "wald_z": z,
        "wald_p": p,
        "ci95_low": beta - 1.96 * se_first,
        "ci95_high": beta + 1.96 * se_first,
    }


def kidney_direction(phenotype, beta):
    if phenotype == "eGFRcrea":
        return "favorable" if beta > 0 else "adverse"

    if phenotype == "BUN":
        return "favorable" if beta < 0 else "adverse"

    return "unknown"


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--anchor", type=Path, required=True)
    ap.add_argument("--egfr", type=Path, required=True)
    ap.add_argument("--bun", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)

    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    anchors = read_tsv(args.anchor)

    if len(anchors) != 9:
        raise SystemExit(
            f"expected 9 pre-specified anchors, found {len(anchors)}"
        )

    anchor_by_rsid = {}

    for row in anchors:
        rsid = row["rsid"]

        if rsid in anchor_by_rsid:
            raise SystemExit(f"duplicate anchor rsID: {rsid}")

        anchor_by_rsid[rsid] = row

    outcome_sets = {
        "eGFRcrea": read_tsv(args.egfr),
        "BUN": read_tsv(args.bun),
    }

    long_rows = []

    for phenotype, outcome_rows in outcome_sets.items():

        outcome_by_rsid = {}

        for r in outcome_rows:
            rsid = r["rsid"]

            if rsid in outcome_by_rsid:
                raise SystemExit(
                    f"duplicate outcome rsID {phenotype}: {rsid}"
                )

            outcome_by_rsid[rsid] = r

        for rsid, a in anchor_by_rsid.items():

            if rsid not in outcome_by_rsid:
                long_rows.append({
                    "gene_symbol": a["gene_symbol"],
                    "protein_id": a["protein_id"],
                    "rsid": rsid,
                    "phenotype": phenotype,
                    "status": "missing_outcome",
                })
                continue

            y = outcome_by_rsid[rsid]

            x_ea = allele(a["effect_allele_alt"])
            x_oa = allele(a["ref_allele"])

            y_ea = allele(y["effect_allele"])
            y_oa = allele(y["other_allele"])

            try:
                sign, orientation = harmonize(
                    x_ea,
                    x_oa,
                    y_ea,
                    y_oa,
                )
            except ValueError as exc:
                long_rows.append({
                    "gene_symbol": a["gene_symbol"],
                    "protein_id": a["protein_id"],
                    "rsid": rsid,
                    "phenotype": phenotype,
                    "status": "harmonization_failed",
                    "error": str(exc),
                })
                continue

            bx = fnum(a["beta_cond_alt"])
            sx = fnum(a["se_cond"])

            raw_by = fnum(y["beta"])
            sy = fnum(y["se"])

            by = raw_by * sign

            mr = wald_ratio(
                beta_x=bx,
                se_x=sx,
                beta_y=by,
                se_y=sy,
            )

            row = {
                "gene_symbol": a["gene_symbol"],
                "protein_id": a["protein_id"],
                "rsid": rsid,
                "phenotype": phenotype,
                "status": "ok",

                "exposure_effect_allele": x_ea,
                "exposure_other_allele": x_oa,
                "exposure_beta": bx,
                "exposure_se": sx,
                "exposure_f_stat": a.get("f_stat", ""),

                "outcome_effect_allele_raw": y_ea,
                "outcome_other_allele_raw": y_oa,
                "outcome_beta_raw": raw_by,
                "outcome_beta_harmonized": by,
                "outcome_se": sy,
                "outcome_p_raw": y["p_value"],

                "harmonization": orientation,
                "harmonization_sign": sign,

                **mr,

                "kidney_direction":
                    kidney_direction(phenotype, mr["wald_beta"]),

                "nominal_p_lt_0_05":
                    int(mr["wald_p"] < 0.05),

                "bonferroni_9_p_lt_0_05":
                    int(mr["wald_p"] < (0.05 / 9.0)),

                "stage2_class":
                    a.get("stage2_class", ""),

                "localization_pattern":
                    a.get("localization_pattern", ""),

                "top_compartment":
                    a.get("top_compartment", ""),

                "top_cell_type":
                    a.get("top_cell_type", ""),
            }

            long_rows.append(row)

    ok_rows = [x for x in long_rows if x.get("status") == "ok"]

    long_fields = [
        "gene_symbol",
        "protein_id",
        "rsid",
        "phenotype",
        "status",

        "exposure_effect_allele",
        "exposure_other_allele",
        "exposure_beta",
        "exposure_se",
        "exposure_f_stat",

        "outcome_effect_allele_raw",
        "outcome_other_allele_raw",
        "outcome_beta_raw",
        "outcome_beta_harmonized",
        "outcome_se",
        "outcome_p_raw",

        "harmonization",
        "harmonization_sign",

        "wald_beta",
        "wald_se_first_order",
        "wald_se_delta",
        "wald_z",
        "wald_p",
        "ci95_low",
        "ci95_high",

        "kidney_direction",
        "nominal_p_lt_0_05",
        "bonferroni_9_p_lt_0_05",

        "stage2_class",
        "localization_pattern",
        "top_compartment",
        "top_cell_type",
    ]

    long_path = args.outdir / "EAS_ANCHOR_WALD_MR.tsv"

    with long_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=long_fields,
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(long_rows)

    by_gene = {}

    for r in ok_rows:
        by_gene.setdefault(
            r["gene_symbol"],
            {
                "gene_symbol": r["gene_symbol"],
                "protein_id": r["protein_id"],
                "rsid": r["rsid"],
                "stage2_class": r["stage2_class"],
                "top_compartment": r["top_compartment"],
                "top_cell_type": r["top_cell_type"],
            },
        )

        prefix = (
            "egfr"
            if r["phenotype"] == "eGFRcrea"
            else "bun"
        )

        by_gene[r["gene_symbol"]].update({
            f"{prefix}_wald_beta": r["wald_beta"],
            f"{prefix}_wald_se": r["wald_se_first_order"],
            f"{prefix}_wald_p": r["wald_p"],
            f"{prefix}_direction": r["kidney_direction"],
            f"{prefix}_harmonization": r["harmonization"],
        })

    integrated = []

    for a in anchors:

        gene = a["gene_symbol"]
        x = by_gene.get(gene, {
            "gene_symbol": gene,
            "protein_id": a["protein_id"],
            "rsid": a["rsid"],
        })

        eb = x.get("egfr_wald_beta")
        bb = x.get("bun_wald_beta")

        if eb is None or bb is None:
            coherence = "incomplete"

        elif eb > 0 and bb < 0:
            coherence = "coherent_favorable"

        elif eb < 0 and bb > 0:
            coherence = "coherent_adverse"

        else:
            coherence = "discordant"

        ep = x.get("egfr_wald_p")
        bp = x.get("bun_wald_p")

        if ep is None or bp is None:
            evidence = "incomplete"

        elif coherence.startswith("coherent") and ep < 0.05 and bp < 0.05:
            evidence = "coherent_both_nominal"

        elif coherence.startswith("coherent") and (ep < 0.05 or bp < 0.05):
            evidence = "coherent_one_nominal"

        elif coherence.startswith("coherent"):
            evidence = "coherent_no_nominal"

        else:
            evidence = "discordant"

        x["cross_trait_coherence"] = coherence
        x["eas_anchor_evidence_class"] = evidence

        integrated.append(x)

    integrated_fields = [
        "gene_symbol",
        "protein_id",
        "rsid",
        "stage2_class",
        "top_compartment",
        "top_cell_type",

        "egfr_harmonization",
        "egfr_wald_beta",
        "egfr_wald_se",
        "egfr_wald_p",
        "egfr_direction",

        "bun_harmonization",
        "bun_wald_beta",
        "bun_wald_se",
        "bun_wald_p",
        "bun_direction",

        "cross_trait_coherence",
        "eas_anchor_evidence_class",
    ]

    integrated_path = (
        args.outdir /
        "EAS_ANCHOR_INTEGRATED_EVIDENCE.tsv"
    )

    with integrated_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=integrated_fields,
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(integrated)

    summary = {
        "stage": "CKD EAS anchor transportability MR",
        "pre_specified_anchor_count": len(anchors),
        "expected_outcome_rows": len(anchors) * 2,
        "successful_harmonized_rows": len(ok_rows),
        "failed_or_missing_rows": len(long_rows) - len(ok_rows),

        "phenotype_counts": {
            p: sum(
                r.get("status") == "ok"
                and r.get("phenotype") == p
                for r in long_rows
            )
            for p in ("eGFRcrea", "BUN")
        },

        "coherence_counts": {},
        "important_interpretation": [
            "This is EUR pQTL anchor -> EAS kidney outcome transportability.",
            "It is not a pure EAS protein MR because the exposure effect estimates are EUR-derived.",
            "Do not select replacement SNPs using EAS outcome significance.",
            "Regional EAS colocalization and EAS/Chinese pQTL replication remain required.",
        ],
    }

    for r in integrated:
        k = r["cross_trait_coherence"]
        summary["coherence_counts"][k] = (
            summary["coherence_counts"].get(k, 0) + 1
        )

    summary_path = args.outdir / "EAS_ANCHOR_MR_SUMMARY.json"

    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(
        "GENE\tRSID\t"
        "eGFR_BETA\teGFR_P\t"
        "BUN_BETA\tBUN_P\t"
        "COHERENCE\tCLASS"
    )

    for r in integrated:
        print(
            r.get("gene_symbol", ""),
            r.get("rsid", ""),
            f'{r.get("egfr_wald_beta", float("nan")):.6g}',
            f'{r.get("egfr_wald_p", float("nan")):.4g}',
            f'{r.get("bun_wald_beta", float("nan")):.6g}',
            f'{r.get("bun_wald_p", float("nan")):.4g}',
            r["cross_trait_coherence"],
            r["eas_anchor_evidence_class"],
            sep="\t",
        )

    print()
    print(json.dumps(summary, indent=2))

    if len(ok_rows) != len(anchors) * 2:
        raise SystemExit(
            "EAS_ANCHOR_MR_INCOMPLETE"
        )

    print("CKD_EAS_ANCHOR_MR_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
