#!/usr/bin/env python3
"""Stage 1 CKD proteome-wide MR screen from compact matched summary data.

Primary screen:
- UKB-PPP Sun 2023 ST16 cis conditional signals as exposure instruments.
- Exposure effect allele = ALT in Variant ID (REF:ALT); beta_cond is ALT effect.
- F = beta_cond^2 / se_cond^2; exclude F < threshold.
- Explicit allele/strand harmonization; conservative palindromic handling.
- Pre-specified strongest-F cis instrument per protein -> Wald ratio; no outcome-specific SNP substitution.
- Fixed-effect IVW across all harmonized ST16 signals is sensitivity only.
- EUR eGFRcrea is primary; BH-FDR across tested proteins.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PAL = {frozenset(("A", "T")), frozenset(("C", "G"))}
OUTCOMES = [
    ("EUR", "eGFRcrea"), ("EUR", "CKD"), ("EUR", "BUN"),
    ("EUR", "eGFRcys"), ("EUR", "UACR"),
    ("EAS", "eGFRcrea"), ("EAS", "BUN"),
]
SUPPORT = [("EUR", "CKD"), ("EUR", "BUN"), ("EUR", "eGFRcys"), ("EUR", "UACR")]
RISK_SIGN = {"eGFRcrea": -1.0, "eGFRcys": -1.0, "CKD": 1.0, "BUN": 1.0, "UACR": 1.0}


def fnum(x):
    if x is None:
        return None
    s = str(x).strip()
    if not s or s.lower() in {"na", "nan", ".", "none"}:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def norm_allele(x):
    return "" if x is None else str(x).strip().upper()


def complement(a):
    return COMP.get(a) if len(a) == 1 else None


def is_pal(a, b):
    return len(a) == 1 and len(b) == 1 and frozenset((a, b)) in PAL


def read_tsv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


def p_norm(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def bh_adjust(ps):
    valid = [(i, p) for i, p in enumerate(ps) if p is not None and math.isfinite(p)]
    out = [None] * len(ps)
    if not valid:
        return out
    valid.sort(key=lambda x: x[1])
    m = len(valid)
    prev = 1.0
    for rank_rev, (i, p) in enumerate(reversed(valid), start=1):
        rank = m - rank_rev + 1
        q = min(prev, p * m / rank, 1.0)
        out[i] = q
        prev = q
    return out


def harmonize(exp, out, maf_threshold, freq_tolerance):
    """Return (status, sign), orienting outcome beta to exposure ALT allele."""
    alt, ref = norm_allele(exp["alt_allele"]), norm_allele(exp["ref_allele"])
    ea, oa = norm_allele(out["effect_allele"]), norm_allele(out["other_allele"])
    if not alt or not ref or not ea or not oa:
        return "missing_allele", None

    if is_pal(alt, ref):
        # Palindromic variants cannot be strand-resolved from allele letters alone,
        # but the outcome allele pair must still be the same biological pair.
        if frozenset((ea, oa)) != frozenset((alt, ref)):
            return "palindrome_allele_mismatch", None
        afx, afo = fnum(exp.get("alt_freq")), fnum(out.get("eaf"))
        if afx is None or afo is None:
            return "palindrome_no_frequency", None
        if min(afx, 1.0 - afx) > maf_threshold:
            return "palindrome_ambiguous_maf", None
        d_keep = abs(afo - afx)
        d_flip = abs(afo - (1.0 - afx))
        if min(d_keep, d_flip) > freq_tolerance:
            return "palindrome_frequency_mismatch", None
        if abs(d_keep - d_flip) < 0.05:
            return "palindrome_frequency_ambiguous", None
        return ("palindrome_freq_keep", 1) if d_keep < d_flip else ("palindrome_freq_flip", -1)

    if (ea, oa) == (alt, ref):
        return "direct", 1
    if (ea, oa) == (ref, alt):
        return "swapped", -1
    cea, coa = complement(ea), complement(oa)
    if cea is not None and coa is not None:
        if (cea, coa) == (alt, ref):
            return "strand", 1
        if (cea, coa) == (ref, alt):
            return "strand_swapped", -1
    return "allele_mismatch", None


def load_instruments(path, f_threshold):
    rows = []
    by_rsid = defaultdict(list)
    anchors = {}
    for r in read_tsv(path):
        bx, sx = fnum(r.get("beta_cond")), fnum(r.get("se_cond"))
        status, fstat = "ok", None
        if bx is None or sx is None or sx <= 0 or bx == 0:
            status = "invalid_beta_se"
        else:
            fstat = (bx / sx) ** 2
            if fstat < f_threshold:
                status = "weak_f"
        x = dict(r)
        x.update(beta_exposure=bx, se_exposure=sx, f_stat=fstat, instrument_qc=status)
        rows.append(x)
        rsid = r.get("rsid", "").strip()
        if rsid and status == "ok":
            by_rsid[rsid].append(x)
            protein = r["protein_id"]
            if protein not in anchors or fstat > anchors[protein]["f_stat"]:
                anchors[protein] = x
    return rows, by_rsid, anchors


def load_outcome(path):
    d = {}
    for r in read_tsv(path):
        rsid = r.get("rsid", "").strip()
        if rsid:
            d[rsid] = r
    return d


def wald(inst):
    bx, by, sy = inst["beta_exposure"], inst["beta_outcome"], inst["se_outcome"]
    beta = by / bx
    se = sy / abs(bx)
    return beta, se, p_norm(beta / se) if se > 0 else 0.0


def ivw(insts):
    good = [x for x in insts if x["se_outcome"] > 0 and x["beta_exposure"] != 0]
    if not good:
        return None, None, None, None
    den = sum((x["beta_exposure"] ** 2) / (x["se_outcome"] ** 2) for x in good)
    if den <= 0:
        return None, None, None, None
    num = sum((x["beta_exposure"] * x["beta_outcome"]) / (x["se_outcome"] ** 2) for x in good)
    beta = num / den
    se = math.sqrt(1.0 / den)
    p = p_norm(beta / se)
    q = sum(
        ((x["beta_outcome"] - beta * x["beta_exposure"]) ** 2) / (x["se_outcome"] ** 2)
        for x in good
    )
    return beta, se, p, q


def outcome_result(ancestry, phenotype, outcome, by_rsid, anchors, maf_threshold, freq_tolerance):
    harmonized, reasons = [], Counter()
    per_protein = defaultdict(list)

    for rsid, inst_rows in by_rsid.items():
        o = outcome.get(rsid)
        if o is None:
            continue
        by, sy = fnum(o.get("beta")), fnum(o.get("se"))
        if by is None or sy is None or sy <= 0:
            reasons["invalid_outcome_beta_se"] += len(inst_rows)
            continue
        for exp in inst_rows:
            status, sign = harmonize(exp, o, maf_threshold, freq_tolerance)
            reasons[status] += 1
            if sign is None:
                continue
            h = {
                "protein_id": exp["protein_id"], "gene_symbol": exp["gene_symbol"],
                "rsid": rsid, "variant_id": exp["variant_id"],
                "chrom_hg19": exp["chrom_hg19"], "pos_hg19": exp["pos_hg19"],
                "exposure_effect_allele": norm_allele(exp["alt_allele"]),
                "exposure_other_allele": norm_allele(exp["ref_allele"]),
                "beta_exposure": exp["beta_exposure"], "se_exposure": exp["se_exposure"],
                "f_stat": exp["f_stat"], "pip": exp.get("pip", ""),
                "alt_freq": exp.get("alt_freq", ""),
                "phenotype": phenotype, "ancestry": ancestry,
                "outcome_effect_allele_original": norm_allele(o.get("effect_allele")),
                "outcome_other_allele_original": norm_allele(o.get("other_allele")),
                "beta_outcome_original": by, "beta_outcome": by * sign,
                "se_outcome": sy, "p_outcome": fnum(o.get("p_value")),
                "outcome_eaf": fnum(o.get("eaf")), "outcome_n": o.get("n", ""),
                "effect_scale": o.get("effect_scale", ""),
                "harmonization": status, "outcome_sign_flip": int(sign == -1),
            }
            harmonized.append(h)
            per_protein[h["protein_id"]].append(h)

    results = []
    for protein, insts in per_protein.items():
        insts.sort(key=lambda x: (x["f_stat"], fnum(x.get("pip")) or -1), reverse=True)
        anchor = anchors.get(protein)
        anchor_rsid = anchor["rsid"] if anchor else ""
        anchor_match = next((x for x in insts if x["rsid"] == anchor_rsid), None)

        # Do not substitute a different SNP when the pre-specified strongest
        # exposure anchor is absent from an outcome. This keeps the primary
        # Wald estimate for a protein on the same IV across phenotypes/ancestries.
        if anchor_match is None:
            continue

        top = anchor_match
        b, s, p = wald(top)
        bivw, sivw, pivw, qivw = ivw(insts)
        results.append({
            "protein_id": protein, "gene_symbol": top["gene_symbol"],
            "ancestry": ancestry, "phenotype": phenotype,
            "n_harmonized_instruments": len(insts),
            "top_rsid": top["rsid"], "top_f_stat": top["f_stat"], "top_pip": top["pip"],
            "wald_beta": b, "wald_se": s, "wald_p": p,
            "wald_or": math.exp(b) if top["effect_scale"] == "log_odds" else "",
            "anchor_rsid": anchor_rsid, "anchor_available": 1,
            "anchor_wald_beta": b, "anchor_wald_se": s, "anchor_wald_p": p,
            "ivw_beta_sensitivity": bivw, "ivw_se_sensitivity": sivw,
            "ivw_p_sensitivity": pivw, "ivw_q_sensitivity": qivw,
            "ivw_df_sensitivity": max(len(insts) - 1, 0),
        })

    results.sort(key=lambda r: (r["wald_p"], r["protein_id"]))
    qs = bh_adjust([r["wald_p"] for r in results])
    m = len(results)
    for r, q in zip(results, qs):
        r["wald_fdr_bh"] = q
        r["bonferroni_threshold"] = 0.05 / m if m else None
        r["bonferroni_significant"] = int(m > 0 and r["wald_p"] < 0.05 / m)
    return harmonized, results, reasons


def build_screen_summary(results_by_key):
    lookup = {k: {r["protein_id"]: r for r in rows} for k, rows in results_by_key.items()}
    primary = lookup.get(("EUR", "eGFRcrea"), {})
    eas = lookup.get(("EAS", "eGFRcrea"), {})
    out = []

    for protein, p in primary.items():
        row = {
            "protein_id": protein, "gene_symbol": p["gene_symbol"],
            "eur_egfr_rsid": p["top_rsid"], "eur_egfr_f": p["top_f_stat"],
            "eur_egfr_beta": p["wald_beta"], "eur_egfr_se": p["wald_se"],
            "eur_egfr_p": p["wald_p"], "eur_egfr_fdr": p["wald_fdr_bh"],
            "eur_egfr_bonf": p["bonferroni_significant"],
        }
        primary_risk_sign = (
            math.copysign(1.0, RISK_SIGN["eGFRcrea"] * p["wald_beta"])
            if p["wald_beta"] != 0 else 0.0
        )
        support_available = support_same = support_p05_same = 0
        for anc, pheno in SUPPORT:
            r = lookup.get((anc, pheno), {}).get(protein)
            prefix = "eur_" + pheno.lower()
            if r:
                row[prefix + "_beta"] = r["wald_beta"]
                row[prefix + "_p"] = r["wald_p"]
                row[prefix + "_rsid"] = r["top_rsid"]
                if pheno == "CKD":
                    row[prefix + "_or"] = r["wald_or"]
                support_available += 1
                risk = RISK_SIGN[pheno] * r["wald_beta"]
                same = (risk == 0 and primary_risk_sign == 0) or (risk * primary_risk_sign > 0)
                if same:
                    support_same += 1
                    if r["wald_p"] < 0.05:
                        support_p05_same += 1
            else:
                row[prefix + "_beta"] = ""
                row[prefix + "_p"] = ""
                row[prefix + "_rsid"] = ""
                if pheno == "CKD":
                    row[prefix + "_or"] = ""
        row["support_available_n"] = support_available
        row["support_same_risk_direction_n"] = support_same
        row["support_p05_same_risk_direction_n"] = support_p05_same

        er = eas.get(protein)
        row["eas_egfr_available"] = int(er is not None)
        row["eas_egfr_beta"] = er["wald_beta"] if er else ""
        row["eas_egfr_p"] = er["wald_p"] if er else ""
        row["eas_egfr_rsid"] = er["top_rsid"] if er else ""
        row["eas_top_sign_concordant"] = int(
            er is not None and er["wald_beta"] * p["wald_beta"] > 0
        )
        row["eas_same_top_snp"] = int(er is not None and er["top_rsid"] == p["top_rsid"])
        row["screen_fdr05"] = int(p["wald_fdr_bh"] < 0.05)
        out.append(row)

    out.sort(key=lambda r: (r["eur_egfr_fdr"], r["eur_egfr_p"], r["protein_id"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ready-root", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--f-threshold", type=float, default=10.0)
    ap.add_argument("--palindrome-maf-threshold", type=float, default=0.42)
    ap.add_argument("--freq-tolerance", type=float, default=0.10)
    args = ap.parse_args()

    inst_path = args.ready_root / "instruments" / "UKBPPP_ST16_cis_independent.tsv.gz"
    if not inst_path.exists():
        alt = args.ready_root / "UKBPPP_ST16_cis_independent.tsv.gz"
        if alt.exists():
            inst_path = alt
        else:
            raise SystemExit(f"instrument file not found: {inst_path} or {alt}")

    inst_rows, by_rsid, anchors = load_instruments(inst_path, args.f_threshold)
    qc_fields = [
        "protein_id","gene_symbol","variant_id","rsid","chrom_hg19","pos_hg19",
        "ref_allele","alt_allele","beta_cond","se_cond","alt_freq","pip",
        "f_stat","instrument_qc",
    ]
    write_tsv(args.output_root / "instrument_qc.tsv.gz", inst_rows, qc_fields)

    results_by_key, harmonization_summary = {}, {}
    h_fields = [
        "protein_id","gene_symbol","rsid","variant_id","chrom_hg19","pos_hg19",
        "exposure_effect_allele","exposure_other_allele","beta_exposure","se_exposure",
        "f_stat","pip","alt_freq","phenotype","ancestry",
        "outcome_effect_allele_original","outcome_other_allele_original",
        "beta_outcome_original","beta_outcome","se_outcome","p_outcome","outcome_eaf",
        "outcome_n","effect_scale","harmonization","outcome_sign_flip",
    ]
    mr_fields = [
        "protein_id","gene_symbol","ancestry","phenotype","n_harmonized_instruments",
        "top_rsid","top_f_stat","top_pip","wald_beta","wald_se","wald_p","wald_or",
        "wald_fdr_bh","bonferroni_threshold","bonferroni_significant",
        "anchor_rsid","anchor_available","anchor_wald_beta","anchor_wald_se","anchor_wald_p",
        "ivw_beta_sensitivity","ivw_se_sensitivity","ivw_p_sensitivity",
        "ivw_q_sensitivity","ivw_df_sensitivity",
    ]

    for ancestry, phenotype in OUTCOMES:
        path = args.ready_root / ancestry / f"{phenotype}.tsv.gz"
        if not path.exists():
            raise SystemExit(f"missing outcome: {path}")
        outcome = load_outcome(path)
        harmonized, results, reasons = outcome_result(
            ancestry, phenotype, outcome, by_rsid, anchors,
            args.palindrome_maf_threshold, args.freq_tolerance,
        )
        results_by_key[(ancestry, phenotype)] = results
        write_tsv(
            args.output_root / "harmonized" / ancestry / f"{phenotype}.tsv.gz",
            harmonized, h_fields,
        )
        write_tsv(
            args.output_root / "mr" / ancestry / f"{phenotype}.tsv",
            results, mr_fields,
        )
        harmonization_summary[f"{ancestry}:{phenotype}"] = {
            "source_rows": len(outcome),
            "harmonized_instrument_rows": len(harmonized),
            "proteins_tested": len(results),
            "harmonization_status_counts": dict(sorted(reasons.items())),
            "fdr05_proteins": sum(r["wald_fdr_bh"] < 0.05 for r in results),
            "bonferroni_proteins": sum(r["bonferroni_significant"] for r in results),
        }

    screen = build_screen_summary(results_by_key)
    screen_fields = list(screen[0].keys()) if screen else []
    write_tsv(args.output_root / "stage1_protein_summary.tsv", screen, screen_fields)
    hits = [r for r in screen if r["screen_fdr05"] == 1]
    write_tsv(args.output_root / "stage1_eur_egfr_fdr05_hits.tsv", hits, screen_fields)

    summary = {
        "stage": "CKD Stage 1 proteome-wide MR screen",
        "primary_ancestry": "EUR",
        "primary_outcome": "eGFRcrea",
        "primary_method": "pre-specified strongest-F cis anchor per protein; Wald ratio; no outcome-specific SNP substitution",
        "sensitivity_method": (
            "fixed-effect IVW across all harmonized Sun ST16 conditional cis signals; "
            "interpret cautiously if residual LD remains"
        ),
        "exposure_effect_allele": "ALT allele from ST16 Variant ID; beta_cond oriented to ALT",
        "instrument_f_threshold": args.f_threshold,
        "palindromic_policy": {
            "maf_threshold": args.palindrome_maf_threshold,
            "frequency_tolerance": args.freq_tolerance,
            "no_outcome_frequency": "drop",
        },
        "instrument_rows": len(inst_rows),
        "instrument_qc_counts": dict(sorted(Counter(r["instrument_qc"] for r in inst_rows).items())),
        "proteins_with_strong_anchor": len(anchors),
        "outcomes": harmonization_summary,
        "primary_proteins_tested": len(screen),
        "primary_fdr05_hits": len(hits),
        "notes": [
            "EUR eGFRcrea BH-FDR is the Stage 1 discovery screen.",
            "Primary Wald estimates use the same pre-specified strongest-F exposure anchor across outcomes; proteins are omitted for an outcome if that anchor is unavailable.",
            "CKD/BUN/eGFRcys/UACR are supporting renal phenotypes, not extra discovery endpoints.",
            "EAS direct-rsID coverage is incomplete and palindromic variants lack EAF in current compact files.",
            "Colocalization is required before causal prioritization; Stage 1 MR hits are screening signals only.",
            "Potential exposure-outcome sample overlap should be assessed during interpretation.",
        ],
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "STAGE1_SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(f"CKD_STAGE1_PASS output={args.output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
