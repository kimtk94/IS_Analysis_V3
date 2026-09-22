#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import date
from pathlib import Path
from statistics import mean

MISSING_STRINGS = {"", "NA", "N/A", ".", "66666", "77777", "88888", "99999"}

FILES = [
    ("TRAIN_F0", "follow_01_base__FOLLOW_01_DATA.txt", "T01_EDATE", "T01_AGE", "T01_CREATININE"),
    ("TRAIN_F1", "follow_02_F1__FOLLOW_02_DATA.txt", "T02_EDATE", "T02_AGE", "T02_CREATININE"),
    ("TRAIN_F2", "follow_03_F2__FOLLOW_03_DATA.txt", "T03_EDATE", "T03_AGE", "T03_CREATININE"),
    ("TRAIN_F3", "follow_04_F3__FOLLOW_04_DATA.txt", "T04_EDATE", "T04_AGE", "T04_CREATININE"),
    ("TRAIN_F4", "follow_05_F4__FOLLOW_05_DATA.txt", "T05_EDATE", "T05_AGE", "T05_CREATININE"),
]


def fnum(x):
    s = str(x or "").strip()
    if s in MISSING_STRINGS:
        return None
    try:
        v = float(s)
    except Exception:
        return None
    if not math.isfinite(v):
        return None
    return v


def parse_yyyymm(x):
    s = str(x or "").strip()
    if not re_digits(s, 6):
        return None
    y, m = int(s[:4]), int(s[4:6])
    if y < 1900 or y > 2100 or m < 1 or m > 12:
        return None
    return date(y, m, 15)


def re_digits(s, n):
    return len(s) == n and s.isdigit()


def decimal_year(d):
    if d is None:
        return None
    start = date(d.year, 1, 1)
    nxt = date(d.year + 1, 1, 1)
    return d.year + (d - start).days / (nxt - start).days


def egfr_ckdepi_2021(scr_mg_dl, age, sex_code):
    if scr_mg_dl is None or age is None or sex_code not in {1, 2}:
        return None
    if scr_mg_dl <= 0 or age < 18 or age > 120:
        return None
    female = sex_code == 2
    kappa = 0.7 if female else 0.9
    alpha = -0.241 if female else -0.302
    ratio = scr_mg_dl / kappa
    val = (
        142.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.200
        * 0.9938 ** age
        * (1.012 if female else 1.0)
    )
    return val


def read_csv(path):
    with Path(path).open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
        return list(csv.DictReader(f))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as f:
        if not fields:
            return
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def ols_slope(points):
    pts = [(x, y) for x, y in points if x is not None and y is not None]
    if len(pts) < 2:
        return None
    xs = [x for x, _ in pts]
    ys = [y for _, y in pts]
    xm, ym = mean(xs), mean(ys)
    den = sum((x - xm) ** 2 for x in xs)
    if den == 0:
        return None
    return sum((x - xm) * (y - ym) for x, y in pts) / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--public-root", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    args = ap.parse_args()

    sex_by_id = {}
    long_rows = []

    # Sex is present in the first follow file only; 1=male, 2=female per official codebook.
    first_path = args.public_root / FILES[0][1]
    first_rows = read_csv(first_path)
    for r in first_rows:
        pid = (r.get("T00_ID") or "").strip()
        sex = fnum(r.get("T00_SEX"))
        if pid and sex in {1.0, 2.0}:
            sex_by_id[pid] = int(sex)

    for wave_index, (wave, fn, date_col, age_col, scr_col) in enumerate(FILES):
        path = args.public_root / fn
        if not path.is_file():
            raise SystemExit(f"missing public training file: {path}")
        rows = read_csv(path)
        for r in rows:
            pid = (r.get("T00_ID") or "").strip()
            if not pid:
                continue
            age = fnum(r.get(age_col))
            scr = fnum(r.get(scr_col))
            d = parse_yyyymm(r.get(date_col))
            sex = sex_by_id.get(pid)
            egfr = egfr_ckdepi_2021(scr, age, sex)
            long_rows.append({
                "participant_id": pid,
                "wave": wave,
                "wave_index": wave_index,
                "exam_yyyymm": (r.get(date_col) or "").strip(),
                "exam_decimal_year": "" if d is None else decimal_year(d),
                "age": "" if age is None else age,
                "sex_code": "" if sex is None else sex,
                "sex_label": "" if sex is None else ("male" if sex == 1 else "female"),
                "creatinine_mg_dl_raw_training": "" if scr is None else scr,
                "egfr_ckdepi2021_raw_training": "" if egfr is None else egfr,
            })

    by_id = {}
    for r in long_rows:
        by_id.setdefault(r["participant_id"], []).append(r)

    subject_rows = []
    for pid, rows in sorted(by_id.items()):
        rows.sort(key=lambda r: r["wave_index"])
        valid = []
        for r in rows:
            x = fnum(r["exam_decimal_year"])
            y = fnum(r["egfr_ckdepi2021_raw_training"])
            if x is not None and y is not None:
                valid.append((x, y, r))
        slope = ols_slope([(x, y) for x, y, _ in valid])

        # Primary baseline is fixed to TRAIN_F0. Do not silently shift baseline
        # to a later wave when F0 eGFR is missing.
        f0 = next(
            (
                (x, y, r)
                for x, y, r in valid
                if r["wave"] == "TRAIN_F0"
            ),
            None,
        )
        baseline_f0 = None if f0 is None else f0[1]
        followup_years = (
            valid[-1][0] - valid[0][0]
            if len(valid) >= 2 else None
        )
        primary_slope_eligible = int(
            len(valid) >= 3
            and followup_years is not None
            and followup_years >= 2.0
        )

        later = [
            (x, y, r)
            for x, y, r in valid
            if r["wave_index"] > 0
        ]
        persistent_low = False
        for i in range(len(later) - 1):
            x1, y1, _ = later[i]
            x2, y2, _ = later[i + 1]
            # KDIGO chronicity requires >=3 months. With this training cohort
            # visits are usually years apart, but retain an explicit >=90-day
            # equivalent threshold for production portability.
            if (
                y1 < 60.0
                and y2 < 60.0
                and (x2 - x1) >= (90.0 / 365.25)
            ):
                persistent_low = True
                break

        incident_primary = int(
            baseline_f0 is not None
            and baseline_f0 >= 60.0
            and persistent_low
        )
        prevalent_ckd_f0 = int(
            baseline_f0 is not None and baseline_f0 < 60.0
        )
        any_later_low = int(
            baseline_f0 is not None
            and baseline_f0 >= 60.0
            and any(y < 60.0 for _, y, _ in later)
        )

        subject_rows.append({
            "participant_id": pid,
            "sex_code": rows[0]["sex_code"] if rows else "",
            "n_waves_total": len(rows),
            "n_egfr_valid": len(valid),
            "baseline_f0_available": int(baseline_f0 is not None),
            "baseline_f0_egfr": "" if baseline_f0 is None else baseline_f0,
            "prevalent_ckd_f0_lt60": prevalent_ckd_f0,
            "last_egfr": "" if not valid else valid[-1][1],
            "followup_years_observed": "" if followup_years is None else followup_years,
            "egfr_slope_ml_min_1.73m2_per_year": "" if slope is None else slope,
            "primary_slope_eligible_ge3_measures_ge2y": primary_slope_eligible,
            "incident_ckd_primary_f0_ge60_persistent_lt60": incident_primary,
            "any_later_egfr_lt60_sensitivity": any_later_low,
        })

    args.output_root.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.output_root / "STAGE4_PUBLIC_EGFR_LONG.tsv",
        long_rows,
    )
    write_tsv(
        args.output_root / "STAGE4_PUBLIC_EGFR_SUBJECT_SUMMARY.tsv",
        subject_rows,
    )

    valid_egfr_n = sum(
        fnum(r["egfr_ckdepi2021_raw_training"]) is not None for r in long_rows
    )
    slope_n = sum(
        fnum(r["egfr_slope_ml_min_1.73m2_per_year"]) is not None
        for r in subject_rows
    )
    primary_slope_n = sum(
        int(r["primary_slope_eligible_ge3_measures_ge2y"])
        for r in subject_rows
    )
    baseline_f0_n = sum(
        int(r["baseline_f0_available"]) for r in subject_rows
    )
    prevalent_f0_n = sum(
        int(r["prevalent_ckd_f0_lt60"]) for r in subject_rows
    )
    strict_events = sum(
        int(r["incident_ckd_primary_f0_ge60_persistent_lt60"])
        for r in subject_rows
    )
    meta = {
        "stage": "CKD Stage 4 public KoGES phenotype prototype",
        "purpose": "workflow QA only; not manuscript-level inference",
        "public_training_root": str(args.public_root),
        "waves": [x[0] for x in FILES],
        "participants": len(subject_rows),
        "long_rows": len(long_rows),
        "valid_egfr_rows": valid_egfr_n,
        "subjects_with_any_slope": slope_n,
        "subjects_primary_slope_eligible_ge3_measures_ge2y": primary_slope_n,
        "subjects_with_f0_egfr": baseline_f0_n,
        "prevalent_ckd_f0_lt60": prevalent_f0_n,
        "prototype_incident_ckd_events_primary": strict_events,
        "egfr_equation": (
            "2021 CKD-EPI creatinine equation: race-free; "
            "sex code 1=male, 2=female."
        ),
        "important_limitations": [
            "Public KoGES training data are for education/prototyping only.",
            "The public BASE_* files use a different synthetic ID namespace and are not joined to FOLLOW_*.",
            "KoGES assay/device creatinine conversion guidance has not been applied here.",
            "Manuscript-level genotype association requires controlled-access individual KoGES data via CODA.",
            "Primary incident-CKD prototype fixes baseline to TRAIN_F0, requires F0 eGFR >=60, and requires two consecutive later eGFR values <60 separated by at least 90 days.",
            "Primary slope eligibility requires at least 3 valid eGFR measurements spanning at least 2 years.",
            "These are workflow definitions for QA, not claims of adjudicated clinical CKD."
        ],
    }
    (args.output_root / "STAGE4_PUBLIC_PROTOTYPE.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(meta, indent=2))
    print("CKD_STAGE4_PUBLIC_PROTOTYPE_PASS")


if __name__ == "__main__":
    main()
