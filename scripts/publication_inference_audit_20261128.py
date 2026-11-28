"""Reproduce sensitivity diagnostics from two pinned public aggregate readbacks.

Not a reanalysis of raw sequencing or a new GWAS. No clinical or causal claims.
"""
import argparse
import hashlib
import json
from math import comb, floor
from pathlib import Path
from statistics import NormalDist

OCULAR_BLOB = "373083ed26abe0c5330bb706f16e7b2f67d915d5"
IS_BLOB = "e505329c7325c9778088ad426a635a57c4748c3b"


def read_pinned(path, expected_blob):
    raw = Path(path).read_bytes()
    actual = hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()
    if actual != expected_blob:
        raise ValueError("pinned aggregate blob mismatch")
    return json.loads(raw)["result"]


def sign_p(positive, negative):
    n = positive + negative
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(comb(n, k) for k in range(min(positive, negative) + 1)) / 2**n)


def holm_two(a, b):
    if a > b:
        x, y = holm_two(b, a)
        return y, x
    return min(1.0, 2*a), min(1.0, max(2*a, b))


def ocular(result):
    panels = ("all_241", "ECM_16")
    for key in panels:
        row = result[key]
        if (row["pairs"], row["positive"], row["negative"], row["ties"]) != (5, 5, 0, 0):
            raise ValueError("unexpected observed pair denominators")
        if row["two_sided_p"] != sign_p(5, 0) or row["holm_p"] != 0.125:
            raise ValueError("observed sign-test parity failed")
    rows = []
    for pos in range(4):
        for neg in range(4-pos):
            p = sign_p(5+pos, neg)
            rows.append({"additional_positive": pos, "additional_negative": neg,
                         "additional_ties": 3-pos-neg, "p": p,
                         "holm_if_both_panels_same": holm_two(p, p)[0]})
    joint = []
    for a in rows:
        for b in rows:
            pa, pb = holm_two(a["p"], b["p"])
            joint.append({"all_241": a, "ECM_16": b, "holm": [pa, pb],
                          "both_below_05": pa < 0.05 and pb < 0.05})
    return {"status": "HYPOTHETICAL_MISSINGNESS_SENSITIVITY_NOT_OBSERVED",
            "observed_pairs": 5, "missing_pairs": 3, "total_pairs": 8,
            "single_panel_scenarios": rows,
            "joint_count_scenarios": len(joint),
            "joint_count_scenarios_both_significant": sum(x["both_below_05"] for x in joint),
            "joint_scenarios": joint,
            "note": "Scenarios are not equally probable. Ties omitted from exact test. Missing nonresponders are not random; no imputation or clinical claim."}


def stroke(result):
    family = result["family"]
    if family != 192 or result["matched_contrasts"] != 158 or result["matched_cells"] != 170:
        raise ValueError("unexpected fixed-family denominators")
    norm = NormalDist()
    z = norm.inv_cdf(1 - 0.05/(2*family))
    rows = []
    hits = result["significant_max_variance"]
    if len(hits) != 4:
        raise ValueError("unexpected conservative-hit count")
    for x in hits:
        delta = x["delta"]
        lo, hi = x["ci95_max_variance"]
        se = (hi-lo)/(2*1.959964)
        if abs((hi+lo)/2-delta) > 1e-7 or se <= 0:
            raise ValueError("invalid archived interval")
        p = 2*norm.cdf(-abs(delta)/se)
        if abs(p-x["p_max_variance"]) > 1e-8:
            raise ValueError("archived p parity failed")
        rows.append({"rsid": x["rsid"], "first": x["first"], "second": x["second"],
                     "delta": delta, "se_reconstructed": se,
                     "simultaneous_interval": [delta-z*se, delta+z*se],
                     "bonferroni_192": x["bonferroni_192_max"],
                     "max_family_size_at_05_diagnostic": floor(0.05/x["p_max_variance"])})
    return {"status": "RETROSPECTIVE_FIXED_FAMILY_CONDITIONAL_DIAGNOSTIC",
            "family": family, "matched_contrasts": result["matched_contrasts"],
            "z": z, "rows": rows,
            "note": "Not a project-wide multiple-testing correction, new association, independent replication or causal estimate. Max-family count is a fragility diagnostic only."}


def self_test():
    assert sign_p(5, 0) == 0.0625
    assert sign_p(8, 0) == 0.0078125
    assert sign_p(7, 1) == 0.0703125
    assert sign_p(6, 2) == 0.2890625
    assert sign_p(5, 3) == 0.7265625
    assert sign_p(0, 0) == 1
    assert holm_two(0.0078125, 0.03125) == (0.015625, 0.03125)
    assert holm_two(0.03125, 0.0078125) == (0.03125, 0.015625)
    assert holm_two(0.03125, 0.03125) == (0.0625, 0.0625)
    assert 3.65 < NormalDist().inv_cdf(1-0.05/(2*192)) < 3.66


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ocular")
    parser.add_argument("--is-readback")
    parser.add_argument("--out")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test and not any((args.ocular, args.is_readback, args.out)):
        print("Focused numerical self-test passed")
        return
    if not all((args.ocular, args.is_readback, args.out)):
        parser.error("both pinned inputs and output path are required")
    result = {"ocular": ocular(read_pinned(args.ocular, OCULAR_BLOB)),
              "is": stroke(read_pinned(args.is_readback, IS_BLOB))}
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"ocular": {k:v for k,v in result["ocular"].items() if k != "joint_scenarios"},
                      "is": result["is"]}, indent=2))


if __name__ == "__main__":
    main()
