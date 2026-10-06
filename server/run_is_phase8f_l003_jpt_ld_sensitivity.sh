#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

P8E="$ROOT/results/is/stage4_cross_eas/phase8e_l003_rs671"
OUT="$ROOT/results/is/stage4_cross_eas/phase8f_l003_jpt_ld"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase8f_l003_jpt_ld/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

PGEN_ROOT="$ROOT/data/is/ld_reference/1kg_eas_grch37/pgen_gt"
PREFIX="$PGEN_ROOT/BBJ_IS_L003.1KG_EAS.GRCh37"

BBJ="$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"

EAS_VARS="$ROOT/results/is/stage3_finemap/japan/bbj/ld_v3/BBJ_IS_L003.unphased.vcor1.bin.vars"

mkdir -p "$OUT" "$LOGDIR"

printf "step\tstatus\trc\n" > "$STATUS"

run_step() {
  local name="$1"
  shift

  local LOG="$LOGDIR/${name}.log"

  echo
  echo "================================================================================"
  echo "$name"
  echo "================================================================================"

  (
    "$@"
  ) > >(tee "$LOG") 2>&1

  local RC=$?

  if [[ $RC -eq 0 ]]; then
    printf "%s\tPASS\t%s\n" "$name" "$RC" >> "$STATUS"
  else
    printf "%s\tFAIL\t%s\n" "$name" "$RC" >> "$STATUS"
  fi

  return 0
}


preflight() {

  echo "===== PHASE8E ====="

  cat "$P8E/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P8E/STEP_STATUS.tsv"
  )"

  echo "PHASE8E_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    return 1
  fi

  echo
  echo "===== 1KG PGEN ====="

  for ext in pgen pvar psam; do
    f="$PREFIX.$ext"

    if [[ -s "$f" ]]; then
      echo "FOUND $f"
    else
      echo "MISSING $f"
      return 1
    fi
  done

  echo
  echo "===== PLINK2 ====="

  PLINK2=""

  for x in \
    "$REPO/tools/plink2_portable/plink2" \
    "$REPO/tools/plink2_latest/plink2" \
    "$(command -v plink2 2>/dev/null)"
  do

    if [[ -n "$x" && -x "$x" ]]; then

      if "$x" --help r-unphased 2>&1 \
        | grep -q 'ref-based'
      then
        PLINK2="$x"
        break
      fi

    fi

  done

  if [[ -z "$PLINK2" ]]; then
    echo "MODERN_PLINK2_NOT_FOUND"
    return 1
  fi

  echo "$PLINK2" > "$OUT/PLINK2_PATH.txt"

  "$PLINK2" --version

  echo "PREFLIGHT=PASS"
}


find_jpt_samples() {

  PANEL="$(
    find "$ROOT/data" \
      -type f \
      -name 'integrated_call_samples_v3.20130502.ALL.panel' \
      -print \
      2>/dev/null \
      | head -1
  )"

  if [[ -z "$PANEL" || ! -s "$PANEL" ]]; then

    echo "1000G_PANEL_NOT_FOUND"
    return 1

  fi

  echo "PANEL=$PANEL"

  python3 - \
    "$PANEL" \
    "$PREFIX.psam" \
    "$OUT/JPT.keep" <<'PY'

from pathlib import Path
import sys

panel = Path(sys.argv[1])
psam = Path(sys.argv[2])
out = Path(sys.argv[3])

jpt = set()

with panel.open(errors="replace") as f:

    header = f.readline().split()

    lower = [
        x.lower()
        for x in header
    ]

    sample_idx = lower.index("sample")
    pop_idx = lower.index("pop")

    for line in f:

        x = line.split()

        if len(x) <= max(
            sample_idx,
            pop_idx
        ):
            continue

        if x[pop_idx] == "JPT":
            jpt.add(
                x[sample_idx]
            )


with psam.open(errors="replace") as f:

    header = f.readline().split()

    rows = [
        x.split()
        for x in f
        if x.strip()
    ]


clean_header = [
    x.lstrip("#")
    for x in header
]

iid_idx = clean_header.index(
    "IID"
)

fid_idx = (
    clean_header.index("FID")
    if "FID" in clean_header
    else None
)


selected = []

for x in rows:

    iid = x[iid_idx]

    if iid not in jpt:
        continue

    if fid_idx is None:

        selected.append(
            (iid,)
        )

    else:

        selected.append(
            (
                x[fid_idx],
                iid,
            )
        )


with out.open("w") as f:

    if fid_idx is None:

        f.write("#IID\n")

        for (iid,) in selected:
            f.write(iid + "\n")

    else:

        f.write("#FID\tIID\n")

        for fid, iid in selected:
            f.write(
                f"{fid}\t{iid}\n"
            )


print(
    "PANEL_JPT=",
    len(jpt)
)

print(
    "PGEN_JPT_INTERSECTION=",
    len(selected)
)


if len(selected) < 80:

    raise SystemExit(
        "Too few JPT samples"
    )

PY

}


build_jpt_ld() {

  PLINK2="$(
    cat "$OUT/PLINK2_PATH.txt"
  )"

  JPT_OUT="$OUT/L003.JPT"

  "$PLINK2" \
    --pfile "$PREFIX" \
    --keep "$OUT/JPT.keep" \
    --extract "$EAS_VARS" \
    --r-unphased \
      square \
      bin \
      ref-based \
    --out "$JPT_OUT"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    echo "JPT_LD_FAIL=$RC"
    return 1
  fi

  MATRIX="$JPT_OUT.unphased.vcor1.bin"

  VARS="$JPT_OUT.unphased.vcor1.bin.vars"

  if [[ ! -s "$MATRIX" || ! -s "$VARS" ]]; then
    echo "JPT_LD_OUTPUT_MISSING"
    return 1
  fi

  echo "JPT_MATRIX=$MATRIX"
  echo "JPT_VARS=$VARS"
}


compare_ld_and_condition() {

python3 - \
  "$OUT" \
  "$BBJ" <<'PY'

from pathlib import Path
import csv
import gzip
import math
import sys

import numpy as np


OUT = Path(sys.argv[1])
BBJ = Path(sys.argv[2])

EAS_INPUT = Path(
    "/srv/is-analysis/results/is/stage4_cross_eas/"
    "phase8e_l003_rs671/"
    "L003_RS671_CONDITIONAL_INPUT.tsv"
)

JPT_MATRIX = (
    OUT /
    "L003.JPT.unphased.vcor1.bin"
)

JPT_VARS = Path(
    str(JPT_MATRIX) + ".vars"
)


with EAS_INPUT.open() as f:

    eas_rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


eas = {
    x["variant_id"]: x
    for x in eas_rows
}


with JPT_VARS.open() as f:

    vids = [
        x.strip().split()[0]
        for x in f
        if x.strip()
    ]


size = JPT_MATRIX.stat().st_size

n2 = size // 8
n = int(round(
    math.sqrt(n2)
))


if n*n*8 != size:
    raise RuntimeError(
        "JPT LD is not float64 square"
    )


if len(vids) == n + 1:

    if vids[0].lstrip("#").upper() in {
        "ID",
        "VARIANT",
        "VARIANT_ID",
    }:
        vids = vids[1:]


if len(vids) != n:
    raise RuntimeError(
        f"vars={len(vids)} matrix={n}"
    )


R = np.memmap(
    JPT_MATRIX,
    dtype="<f8",
    mode="r",
    shape=(n,n),
)


RS671 = "12:112241766:G:A"


if RS671 not in vids:
    raise RuntimeError(
        "rs671 absent from JPT LD"
    )


k = vids.index(
    RS671
)


# BBJ z-scores
wanted = set(vids)
summary = {}

with gzip.open(
    BBJ,
    "rt",
    errors="replace"
) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for x in reader:

        vid = x["variant_id"]

        if vid not in wanted:
            continue

        summary[vid] = {
            "z":
                float(x["beta"])
                /
                float(x["se"]),
            "beta":
                float(x["beta"]),
            "se":
                float(x["se"]),
            "p":
                float(x["p"]),
            "rsid":
                x["rsid"],
        }


z671 = summary[
    RS671
]["z"]


rows = []


for i, vid in enumerate(
    vids
):

    if vid == RS671:
        continue

    if vid not in summary:
        continue

    r = float(
        R[i,k]
    )

    r2 = r*r

    denom2 = 1.0 - r2

    if (
        not np.isfinite(r)
        or denom2 <= 1e-6
    ):
        continue

    zraw = summary[
        vid
    ]["z"]

    zcond = (
        zraw
        - r*z671
    ) / math.sqrt(
        denom2
    )

    pcond = math.erfc(
        abs(zcond)
        / math.sqrt(2)
    )

    eas_row = eas.get(
        vid,
        {}
    )

    rows.append({
        "variant_id":
            vid,

        "rsid":
            summary[vid][
                "rsid"
            ],

        "z_raw":
            zraw,

        "jpt_r_to_rs671":
            r,

        "jpt_r2_to_rs671":
            r2,

        "jpt_z_cond":
            zcond,

        "jpt_p_cond":
            pcond,

        "eas_r_to_rs671":
            eas_row.get(
                "r_to_rs671",
                ""
            ),

        "eas_r2_to_rs671":
            eas_row.get(
                "r2_to_rs671",
                ""
            ),

        "eas_z_cond":
            eas_row.get(
                "z_cond_rs671",
                ""
            ),

        "eas_p_cond":
            eas_row.get(
                "p_cond_rs671",
                ""
            ),
    })


rows.sort(
    key=lambda x:
        x["jpt_p_cond"]
)


with (
    OUT /
    "L003_JPT_CONDITIONAL_VARIANTS.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            rows[0].keys()
        )
    )

    w.writeheader()
    w.writerows(rows)


targets = [
    "12:112337924:G:A",
    "12:112110489:T:C",
    "12:112574616:C:T",
]


print(
    "JPT_LD_N=",
    n
)

print(
    "RS671_Z=",
    z671
)


print()
print(
    "===== KEY PROXY COMPARISON ====="
)

for vid in targets:

    x = next(
        (
            r
            for r in rows
            if r["variant_id"] == vid
        ),
        None
    )

    print(
        vid,
        x
    )


print()
print(
    "===== JPT TOP 20 ====="
)

for x in rows[:20]:

    print(
        x["variant_id"],
        x["rsid"],
        "JPT_r2=",
        x["jpt_r2_to_rs671"],
        "JPT_p=",
        x["jpt_p_cond"],
        "EAS_r2=",
        x["eas_r2_to_rs671"],
        "EAS_p=",
        x["eas_p_cond"],
    )


gws = [
    x
    for x in rows
    if x["jpt_p_cond"] < 5e-8
]


low_ld_gws = [
    x
    for x in gws
    if x["jpt_r2_to_rs671"] < 0.8
]


summary_out = [{
    "jpt_variants":
        len(rows),

    "jpt_conditional_gws_n":
        len(gws),

    "jpt_conditional_low_ld_gws_n":
        len(low_ld_gws),

    "jpt_top_variant":
        rows[0]["variant_id"],

    "jpt_top_p":
        rows[0]["jpt_p_cond"],

    "jpt_top_r2_to_rs671":
        rows[0]["jpt_r2_to_rs671"],

    "interpretation":
        (
            "ROBUST_SECONDARY_CANDIDATE"
            if low_ld_gws
            else
            "NO_LOW_LD_GENOMEWIDE_SECONDARY_SIGNAL"
        ),
}]


with (
    OUT /
    "L003_JPT_CONDITIONING_SUMMARY.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            summary_out[0].keys()
        )
    )

    w.writeheader()
    w.writerows(
        summary_out
    )


print()
print(
    "===== SUMMARY ====="
)

print(
    summary_out[0]
)

PY

}


build_readiness() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])

p = (
    OUT /
    "L003_JPT_CONDITIONING_SUMMARY.tsv"
)


if not p.exists():

    status = (
        "JPT_ANALYSIS_INCOMPLETE"
    )

    secondary = (
        "UNRESOLVED"
    )

else:

    with p.open() as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    x = rows[0]

    n_low = int(
        x[
            "jpt_conditional_low_ld_gws_n"
        ]
    )

    if n_low > 0:

        secondary = (
            "SUPPORTED_BY_JPT_SENSITIVITY"
        )

    else:

        secondary = (
            "NOT_SUPPORTED_AS_LOW_LD_INDEPENDENT_SIGNAL"
        )

    status = (
        "COMPUTATIONAL_COMPLETE"
    )


rows = [
    (
        "PHASE8F_JPT_LD",
        status
    ),

    (
        "PRIMARY_SIGNAL",
        "RS671"
    ),

    (
        "EAS504_CONDITIONAL_RESULT",
        "GWS_RESIDUAL_IN_NEAR_PERFECT_RS671_PROXIES"
    ),

    (
        "JPT_LD_SENSITIVITY",
        secondary
    ),

    (
        "INDEPENDENT_SECONDARY_SIGNAL",
        secondary
    ),

    (
        "L003_PAPER_GRADE",
        "PENDING_COHORT_MATCHED_OR_INTERNAL_BBJ_LD"
    ),

    (
        "ALDH2_INTERPRETATION",
        "RS671_DOMINANT_EAS_LOCUS_REMAINS_PRIMARY"
    ),
]


with (
    OUT /
    "PHASE8F_READINESS.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "component",
        "status"
    ])

    w.writerows(rows)


for a,b in rows:
    print(
        f"{a} = {b}"
    )

PY

}


run_step \
  "01_preflight" \
  preflight

run_step \
  "02_find_jpt_samples" \
  find_jpt_samples

run_step \
  "03_build_jpt_ld" \
  build_jpt_ld

run_step \
  "04_compare_ld_and_condition" \
  compare_ld_and_condition

run_step \
  "05_build_readiness" \
  build_readiness


echo
echo "================================================================================"
echo "PHASE8F COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== JPT CONDITIONING SUMMARY ====="

cat \
  "$OUT/L003_JPT_CONDITIONING_SUMMARY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== KEY VARIANTS ====="

python3 - \
  "$OUT/L003_JPT_CONDITIONAL_VARIANTS.tsv" <<'PY'

import csv
import sys

p = sys.argv[1]

targets = {
    "12:112337924:G:A",
    "12:112110489:T:C",
    "12:112574616:C:T",
}

try:
    with open(p) as f:

        for x in csv.DictReader(
            f,
            delimiter="\t"
        ):

            if x[
                "variant_id"
            ] in targets:

                print(x)

except FileNotFoundError:
    pass

PY


echo
echo "===== READINESS ====="

cat \
  "$OUT/PHASE8F_READINESS.tsv" \
  2>/dev/null \
  || true


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
