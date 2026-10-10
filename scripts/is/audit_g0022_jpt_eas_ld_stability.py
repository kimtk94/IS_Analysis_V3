#!/usr/bin/env python3
"""G0022 rs671 JPT-vs-other-EAS genotype LD sensitivity; NOT cohort-specific GWAS LD."""
import argparse
import csv
import hashlib
import json
import math
import struct
from collections import Counter
from pathlib import Path

import numpy as np

LEAD = "12:112241766:G:A"
POPS = ("EAS", "JPT", "CHB", "CHS", "CDX", "KHV", "EAS_NON_JPT")


def load_panel(panel):
    result = {}
    with Path(panel).open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["super_pop"] == "EAS":
                result[row["sample"]] = row["pop"]
    return result


def get_reference_markers(variants_path, matrix_path, cutoff=0.5, lead=LEAD):
    with Path(variants_path).open(newline="") as f:
        v = list(csv.DictReader(f, delimiter="\t"))
    ids = [r["variant_id"] for r in v]
    if len(ids) != len(set(ids)) or lead not in ids:
        raise ValueError("LD variant alignment invalid")
    n = len(ids)
    path = Path(matrix_path)
    if path.stat().st_size != n * n * 8:
        raise ValueError("LD dimension mismatch")
    with path.open("rb") as f:
        f.seek(ids.index(lead) * n * 8)
        row = struct.unpack(f"<{n}d", f.read(n * 8))
    if not math.isclose(row[ids.index(lead)], 1.0, abs_tol=1e-6):
        raise ValueError("LD lead diagonal mismatch")
    selected = {ids[i]: {"r_eas_original": row[i], "r2_eas_original": row[i] ** 2}
                for i in range(n) if row[i] ** 2 >= cutoff - 1e-12}
    if lead not in selected:
        raise ValueError("lead not selected")
    return selected


def load_dosages(traw, selected, sample_pops):
    matched = {}
    with Path(traw).open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        head = next(reader)
        if head[:6] != ["CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"]:
            raise ValueError("unexpected PLINK .traw format")
        samples = [x.split("_")[0] for x in head[6:]]
        if len(set(samples)) != len(samples):
            raise ValueError("duplicate sample IDs in .traw")
        if len(samples) != 504 and len(sample_pops) > 100:
            raise ValueError("expected 504 EAS samples")
        if any(s not in sample_pops for s in samples):
            raise ValueError("a .traw sample is not from documented EAS panel")
        for row in reader:
            var = row[1]
            if var not in selected:
                continue
            bits = var.split(":")
            if len(bits) != 4 or row[4].upper() != bits[2].upper() or row[5].upper() != bits[3].upper():
                raise ValueError(f"allele mismatch in genotype traw for {var}")
            if var in matched:
                raise ValueError(f"duplicate .traw marker {var}")
            values = [float(x) if x not in ("NA", ".", "") else math.nan for x in row[6:]]
            if len(values) != len(samples):
                raise ValueError("missing genotype dosage column")
            # PLINK2 traw COUNTED allele is REF in this input. ALT = 2 - COUNTED.
            matched[var] = 2.0 - np.asarray(values, dtype=float)
    if set(matched) != set(selected):
        raise ValueError(f"unmatched genotype variants: {len(set(selected) - set(matched))}")
    return samples, matched


def corr2(a, b):
    use = np.isfinite(a) & np.isfinite(b)
    n = int(use.sum())
    if n < 4:
        return None, n
    x, y = a[use], b[use]
    if np.std(x) <= 1e-12 or np.std(y) <= 1e-12:
        return None, n
    corr = float(np.corrcoef(x, y)[0, 1])
    return corr * corr, n


def summarize(samples, dosages, selected, sample_pops, lead=LEAD):
    groups = {pop: [] for pop in POPS}
    for i, sample in enumerate(samples):
        pop = sample_pops[sample]
        groups["EAS"].append(i)
        groups[pop].append(i)
        if pop != "JPT":
            groups["EAS_NON_JPT"].append(i)
    summary = {"audit": "IS_G0022_RS671_JPT_EAS_LD_STABILITY_V1",
               "status": "GENOTYPE_SUBGROUP_LD_SENSITIVITY_NOT_STUDY_GWAS_FINE_MAPPING",
               "lead": lead, "marker_count": len(selected),
               "samples": {pop: len(index) for pop, index in groups.items()},
               "assumption": "PLINK .traw counted allele is REF, ALT dosage = 2 - COUNTED",
               "groups": {},
               "scientific_boundaries": [
                 "JPT104 and other EAS reference samples are nested in EAS504, not independent replication.",
                 "Population-specific LD variation does not replace GIGASTROKE cohort LD, N_eff, or SNP-level QC.",
                 "Monomorphic genotype markers are NA LD, not r2=0.",
               ]}
    rows = []
    max_mismatch = 0.0
    lead_all = dosages[lead]
    for pop in POPS:
        idx = groups[pop]
        if not idx:
            continue
        a = lead_all[idx]
        finite = a[np.isfinite(a)]
        summary["groups"][pop] = {"n": len(idx),
                                    "rs671_alt_af": float(finite.mean() / 2) if len(finite) else None,
                                    "rs671_nonmissing_n": len(finite)}
    for var, orig in selected.items():
        rec = {"variant_id": var, "eas_original_r2": round(orig["r2_eas_original"], 9)}
        for pop in POPS:
            idx = groups[pop]
            if not idx:
                continue
            pair_r2, n = corr2(lead_all[idx], dosages[var][idx])
            freqs = dosages[var][idx]
            valid = freqs[np.isfinite(freqs)]
            rec[f"{pop}_r2"] = round(pair_r2, 9) if pair_r2 is not None else ""
            rec[f"{pop}_n_pair"] = n
            rec[f"{pop}_alt_af"] = round(float(valid.mean() / 2), 9) if len(valid) else ""
            if pop == "EAS" and pair_r2 is not None:
                max_mismatch = max(max_mismatch, abs(pair_r2 - orig["r2_eas_original"]))
        rows.append(rec)
    summary["eas_full_n_original_ld_max_abs_r2_diff"] = max_mismatch
    if max_mismatch > 0.005:
        raise ValueError(f"source dosage vs signed LD disagreement (max delta r2={max_mismatch})")
    nonlead = [x for x in rows if x["variant_id"] != lead]
    valid = [x for x in nonlead if x.get("JPT_r2") != "" and x.get("EAS_r2") != ""]
    summary["n_nonlead_with_jpt_and_eas_ld_defined"] = len(valid)
    summary["jpt_vs_eas_max_abs_r2_diff"] = max((abs(x["JPT_r2"] - x["EAS_r2"]) for x in valid), default=None)
    summary["jpt_vs_eas_median_abs_r2_diff"] = float(np.median([abs(x["JPT_r2"] - x["EAS_r2"]) for x in valid])) if valid else None
    return rows, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traw", type=Path, required=True)
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--variants", type=Path, required=True)
    ap.add_argument("--eas-signed-ld", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--minimum-eas-r2", type=float, default=0.50)
    a = ap.parse_args()
    if not (0 < a.minimum_eas_r2 <= 1):
        raise ValueError("invalid R2 cutoff")
    panel = load_panel(a.panel)
    selected = get_reference_markers(a.variants, a.eas_signed_ld, a.minimum_eas_r2)
    samples, genotypes = load_dosages(a.traw, selected, panel)
    rows, summary = summarize(samples, genotypes, selected, panel)
    summary["source_paths"] = {k: str(v.resolve()) for k, v in (
        ("traw", a.traw), ("panel", a.panel), ("variants", a.variants), ("eas_signed_ld", a.eas_signed_ld))}
    a.outdir.mkdir(parents=True, exist_ok=True)
    out = a.outdir / "G0022_RS671_SUBGROUP_LD_SENSITIVITY.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (a.outdir / "G0022_RS671_SUBGROUP_LD_SENSITIVITY_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
