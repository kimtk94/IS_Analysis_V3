#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

R2="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit_r2"

OUT="$ROOT/results/is/stage4_cross_eas/phase8e_l003_rs671"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase8e_l003_rs671/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

BBJ="$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"

LD="$ROOT/results/is/stage3_finemap/japan/bbj/ld_v3/BBJ_IS_L003.unphased.vcor1.bin"
VARS="$ROOT/results/is/stage3_finemap/japan/bbj/ld_v3/BBJ_IS_L003.unphased.vcor1.bin.vars"

GIGA="$ROOT/results/is/stage4_cross_eas/GIGASTROKE_BBJ_HARMONIZED_VARIANTS.tsv"

RS671="12:112241766:G:A"

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


# =============================================================================
# 01. PREFLIGHT
# =============================================================================

preflight() {

  echo "===== PHASE8D-R2 ====="

  if [[ ! -s "$R2/STEP_STATUS_R2.tsv" ]]; then
    echo "R2_STATUS_MISSING"
    return 1
  fi

  cat "$R2/STEP_STATUS_R2.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$R2/STEP_STATUS_R2.tsv"
  )"

  echo
  echo "R2_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    echo "R2_NOT_CLEAN"
    return 1
  fi

  echo
  echo "===== INPUTS ====="

  for f in \
    "$BBJ" \
    "$LD" \
    "$VARS"
  do

    if [[ -s "$f" ]]; then
      echo "FOUND $f $(stat -c%s "$f")"
    else
      echo "MISSING $f"
      return 1
    fi

  done

  echo
  echo "===== DEPENDENCIES ====="

  for exe in \
    python3 \
    Rscript
  do

    if command -v "$exe" >/dev/null 2>&1; then
      echo "FOUND $exe $(command -v "$exe")"
    else
      echo "MISSING $exe"
      return 1
    fi

  done

  Rscript - <<'RS'
if (!requireNamespace("susieR", quietly=TRUE)) {
    stop("susieR missing")
}
cat(
    "susieR=",
    as.character(packageVersion("susieR")),
    "\n",
    sep=""
)
RS

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. BUILD FRESH L003 INPUT + EXACT rs671 CONDITIONAL STATISTICS
# =============================================================================

build_conditioned_input() {

python3 - \
  "$BBJ" \
  "$LD" \
  "$VARS" \
  "$OUT" \
  "$RS671" <<'PY'

from pathlib import Path
import csv
import gzip
import math
import sys

import numpy as np


BBJ = Path(sys.argv[1])
LD = Path(sys.argv[2])
VARS = Path(sys.argv[3])
OUT = Path(sys.argv[4])
RS671 = sys.argv[5]

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# -------------------------------------------------------------------
# Resolve matrix dimensions independently from metadata.
# -------------------------------------------------------------------

actual_bytes = LD.stat().st_size

if actual_bytes % 8 != 0:
    raise RuntimeError(
        f"LD bytes not divisible by 8: {actual_bytes}"
    )

n2 = actual_bytes // 8
n = int(round(math.sqrt(n2)))

if n * n != n2:
    raise RuntimeError(
        f"LD is not square float64: bytes={actual_bytes}"
    )

print("LD_N=", n)
print("LD_BYTES=", actual_bytes)


# -------------------------------------------------------------------
# Parse matrix variable order.
# -------------------------------------------------------------------

with VARS.open(
    errors="replace"
) as f:

    raw_vars = [
        x.strip().split()[0]
        for x in f
        if x.strip()
    ]


if len(raw_vars) == n + 1:

    first = raw_vars[0].lstrip("#").upper()

    if first in {
        "ID",
        "VARIANT",
        "VARIANT_ID",
    }:
        raw_vars = raw_vars[1:]


if len(raw_vars) != n:

    raise RuntimeError(
        f"VARS count={len(raw_vars)} != LD n={n}"
    )


if len(set(raw_vars)) != len(raw_vars):

    raise RuntimeError(
        "duplicate variant IDs in LD vars"
    )


if RS671 not in raw_vars:

    raise RuntimeError(
        f"rs671 variant absent from LD vars: {RS671}"
    )


rs_index = raw_vars.index(
    RS671
)

print(
    "RS671_LD_INDEX=",
    rs_index
)


# -------------------------------------------------------------------
# Structural QC again on exact analysis matrix.
# NO symmetrisation.
# -------------------------------------------------------------------

R = np.memmap(
    LD,
    dtype="<f8",
    mode="r",
    shape=(n, n),
)

diag = np.asarray(
    R[
        np.arange(n),
        np.arange(n)
    ],
    dtype=np.float64,
)

max_diag_dev = float(
    np.max(
        np.abs(
            diag - 1.0
        )
    )
)

max_asym = 0.0
nonfinite = 0

BLOCK = 256

for i0 in range(
    0,
    n,
    BLOCK
):

    i1 = min(
        n,
        i0 + BLOCK
    )

    A = np.asarray(
        R[i0:i1, :],
        dtype=np.float64,
    )

    B = np.asarray(
        R[:, i0:i1],
        dtype=np.float64,
    ).T

    nonfinite += int(
        A.size
        - np.isfinite(A).sum()
    )

    d = np.abs(
        A - B
    )

    if np.isfinite(d).any():

        max_asym = max(
            max_asym,
            float(
                np.nanmax(d)
            )
        )


print(
    "LD_MAX_ASYM=",
    max_asym
)

print(
    "LD_MAX_DIAG_DEV=",
    max_diag_dev
)

print(
    "LD_NONFINITE=",
    nonfinite
)


if (
    nonfinite != 0
    or max_asym > 1e-10
    or max_diag_dev > 1e-10
):
    raise RuntimeError(
        "exact L003 LD structural QC failed"
    )


# -------------------------------------------------------------------
# Read only BBJ rows present in LD matrix.
# -------------------------------------------------------------------

wanted = set(
    raw_vars
)

rows = {}

with gzip.open(
    BBJ,
    "rt",
    errors="replace",
) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t",
    )

    required = {
        "variant_id",
        "chr",
        "pos",
        "rsid",
        "ref",
        "alt",
        "effect_allele",
        "other_allele",
        "beta",
        "se",
        "p",
        "eaf",
    }

    missing = required - set(
        reader.fieldnames or []
    )

    if missing:
        raise RuntimeError(
            f"BBJ missing columns: {sorted(missing)}"
        )

    for r in reader:

        vid = r["variant_id"]

        if vid not in wanted:
            continue

        if vid in rows:
            raise RuntimeError(
                f"duplicate BBJ variant {vid}"
            )

        rows[vid] = r


missing_summary = [
    x
    for x in raw_vars
    if x not in rows
]

print(
    "SUMMARY_MATCHED=",
    len(rows)
)

print(
    "SUMMARY_MISSING=",
    len(missing_summary)
)


if missing_summary:

    with (
        OUT /
        "L003_LD_VARIANTS_MISSING_FROM_BBJ.txt"
    ).open("w") as f:

        for x in missing_summary:
            f.write(x + "\n")

    raise RuntimeError(
        "Not all LD variants were recovered from BBJ canonical"
    )


# -------------------------------------------------------------------
# Exact allele audit.
# -------------------------------------------------------------------

allele_errors = []

beta = np.zeros(
    n,
    dtype=float,
)

se = np.zeros(
    n,
    dtype=float,
)

p_raw = np.zeros(
    n,
    dtype=float,
)

eaf = np.zeros(
    n,
    dtype=float,
)

rsids = []


for i, vid in enumerate(
    raw_vars
):

    r = rows[vid]

    parts = vid.split(
        ":",
        3
    )

    if len(parts) != 4:
        raise RuntimeError(
            f"unexpected variant_id: {vid}"
        )

    chrom, pos, ref, alt = parts

    obs = (
        str(r["chr"]),
        str(
            int(
                float(
                    r["pos"]
                )
            )
        ),
        r["ref"].upper(),
        r["alt"].upper(),
    )

    exp = (
        chrom,
        pos,
        ref.upper(),
        alt.upper(),
    )

    effect = (
        r["effect_allele"]
        .strip()
        .upper()
    )

    other = (
        r["other_allele"]
        .strip()
        .upper()
    )

    if (
        obs != exp
        or effect != exp[3]
        or other != exp[2]
    ):

        allele_errors.append({
            "variant_id": vid,
            "expected": ":".join(exp),
            "observed":
                ":".join(obs),
            "effect_allele": effect,
            "other_allele": other,
        })

    beta[i] = float(
        r["beta"]
    )

    se[i] = float(
        r["se"]
    )

    p_raw[i] = float(
        r["p"]
    )

    eaf[i] = float(
        r["eaf"]
    )

    rsids.append(
        r["rsid"]
    )


print(
    "ALLELE_ERRORS=",
    len(allele_errors)
)


if allele_errors:

    with (
        OUT /
        "L003_ALLELE_ERRORS.tsv"
    ).open(
        "w",
        newline=""
    ) as f:

        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=list(
                allele_errors[0].keys()
            ),
        )

        w.writeheader()
        w.writerows(
            allele_errors
        )

    raise RuntimeError(
        "BBJ allele orientation mismatch"
    )


if np.any(
    ~np.isfinite(beta)
):

    raise RuntimeError(
        "non-finite beta"
    )


if np.any(
    ~np.isfinite(se)
) or np.any(
    se <= 0
):

    raise RuntimeError(
        "invalid SE"
    )


z = beta / se

k = rs_index

z_k = float(
    z[k]
)

p_k = float(
    p_raw[k]
)

beta_k = float(
    beta[k]
)

se_k = float(
    se[k]
)

eaf_k = float(
    eaf[k]
)

r_to_k = np.asarray(
    R[:, k],
    dtype=np.float64,
)

r2_to_k = (
    r_to_k ** 2
)


print()
print(
    "RS671_VARIANT=",
    raw_vars[k]
)

print(
    "RS671_RSID=",
    rsids[k]
)

print(
    "RS671_BETA=",
    beta_k
)

print(
    "RS671_SE=",
    se_k
)

print(
    "RS671_Z=",
    z_k
)

print(
    "RS671_P=",
    p_k
)

print(
    "RS671_EAF=",
    eaf_k
)


# -------------------------------------------------------------------
# Summary-statistic single-SNP conditioning.
#
# Under standardized-score approximation:
#
# z_j | k =
#   (z_j - r_jk z_k) /
#   sqrt(1 - r_jk^2)
#
# This is explicitly exploratory because:
#   - LD = external 1000G EAS 504
#   - GWAS = binary trait
#   - individual-level BBJ data unavailable
#
# -------------------------------------------------------------------

variance_left = (
    1.0
    - r2_to_k
)

keep = np.ones(
    n,
    dtype=bool,
)

keep[k] = False

# Perfect / nearly-perfect proxies cannot be separately conditioned.
keep &= (
    variance_left > 1e-6
)

keep &= np.isfinite(
    variance_left
)

idx = np.where(
    keep
)[0]


print()
print(
    "CONDITIONAL_KEEP=",
    len(idx)
)

print(
    "DROPPED_SELF_OR_R2_NEAR1=",
    n - len(idx)
)


denom = np.sqrt(
    variance_left[idx]
)

z_cond = (
    z[idx]
    - r_to_k[idx] * z_k
) / denom


def normal_two_sided(zv):

    return math.erfc(
        abs(float(zv))
        / math.sqrt(2.0)
    )


p_cond = np.array(
    [
        normal_two_sided(x)
        for x in z_cond
    ],
    dtype=float,
)


# -------------------------------------------------------------------
# Conditional correlation matrix.
#
# R_ij | k =
# (R_ij - r_ik r_jk) /
# sqrt((1-r_ik^2)(1-r_jk^2))
#
# NO symmetrisation.
# -------------------------------------------------------------------

R_sub = np.asarray(
    R[
        np.ix_(
            idx,
            idx
        )
    ],
    dtype=np.float64,
)

rk = r_to_k[
    idx
]

R_cond = (
    R_sub
    - np.outer(
        rk,
        rk
    )
) / (
    denom[:, None]
    * denom[None, :]
)


cond_asym = float(
    np.max(
        np.abs(
            R_cond
            - R_cond.T
        )
    )
)

cond_diag_dev = float(
    np.max(
        np.abs(
            np.diag(
                R_cond
            )
            - 1.0
        )
    )
)

cond_nonfinite = int(
    np.size(R_cond)
    - np.isfinite(
        R_cond
    ).sum()
)


print()
print(
    "COND_R_MAX_ASYM=",
    cond_asym
)

print(
    "COND_R_MAX_DIAG_DEV=",
    cond_diag_dev
)

print(
    "COND_R_NONFINITE=",
    cond_nonfinite
)


if (
    cond_nonfinite != 0
    or cond_asym > 1e-8
    or cond_diag_dev > 1e-7
):
    raise RuntimeError(
        "derived conditional R failed QC"
    )


cond_bin = (
    OUT /
    "L003_RS671_CONDITIONAL_R.float64.bin"
)

R_cond.astype(
    "<f8",
    copy=False
).tofile(
    cond_bin
)


# -------------------------------------------------------------------
# Save ordered SuSiE input.
# -------------------------------------------------------------------

ordered = []

for j, original_i in enumerate(
    idx
):

    vid = raw_vars[
        original_i
    ]

    r = rows[vid]

    ordered.append({
        "matrix_index_1based":
            j + 1,
        "original_ld_index_1based":
            original_i + 1,
        "variant_id":
            vid,
        "rsid":
            r["rsid"],
        "ref":
            r["ref"],
        "alt":
            r["alt"],
        "effect_allele":
            r["effect_allele"],
        "beta":
            beta[original_i],
        "se":
            se[original_i],
        "p_raw":
            p_raw[original_i],
        "eaf":
            eaf[original_i],
        "z_raw":
            z[original_i],
        "r_to_rs671":
            r_to_k[original_i],
        "r2_to_rs671":
            r2_to_k[original_i],
        "z_cond_rs671":
            z_cond[j],
        "p_cond_rs671":
            p_cond[j],
    })


input_tsv = (
    OUT /
    "L003_RS671_CONDITIONAL_INPUT.tsv"
)

with input_tsv.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            ordered[0].keys()
        ),
    )

    w.writeheader()
    w.writerows(
        ordered
    )


# Sorted version for biological interpretation.
ranked = sorted(
    ordered,
    key=lambda x:
        float(
            x["p_cond_rs671"]
        ),
)


with (
    OUT /
    "L003_RS671_CONDITIONAL_VARIANTS_RANKED.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            ranked[0].keys()
        ),
    )

    w.writeheader()
    w.writerows(
        ranked
    )


# -------------------------------------------------------------------
# Summary.
# -------------------------------------------------------------------

raw_best_i = int(
    np.nanargmax(
        np.abs(z)
    )
)

cond_best = ranked[0]

genomewide_n = sum(
    float(x["p_cond_rs671"])
    < 5e-8
    for x in ranked
)

suggestive_n = sum(
    float(x["p_cond_rs671"])
    < 1e-5
    for x in ranked
)

r2_08 = int(
    np.sum(
        r2_to_k >= 0.8
    )
)

r2_05 = int(
    np.sum(
        r2_to_k >= 0.5
    )
)

r2_01 = int(
    np.sum(
        r2_to_k >= 0.1
    )
)


summary = {
    "locus":
        "BBJ_IS_L003",

    "condition_variant":
        RS671,

    "condition_rsid":
        rsids[k],

    "condition_beta":
        beta_k,

    "condition_se":
        se_k,

    "condition_z":
        z_k,

    "condition_p":
        p_k,

    "condition_eaf":
        eaf_k,

    "ld_variants":
        n,

    "conditional_variants":
        len(idx),

    "r2_ge_0.8_to_rs671":
        r2_08,

    "r2_ge_0.5_to_rs671":
        r2_05,

    "r2_ge_0.1_to_rs671":
        r2_01,

    "raw_top_variant":
        raw_vars[
            raw_best_i
        ],

    "raw_top_p":
        p_raw[
            raw_best_i
        ],

    "conditional_top_variant":
        cond_best[
            "variant_id"
        ],

    "conditional_top_rsid":
        cond_best[
            "rsid"
        ],

    "conditional_top_z":
        cond_best[
            "z_cond_rs671"
        ],

    "conditional_top_p":
        cond_best[
            "p_cond_rs671"
        ],

    "conditional_p_lt_5e8_n":
        genomewide_n,

    "conditional_p_lt_1e5_n":
        suggestive_n,

    "conditional_R_max_asym":
        cond_asym,

    "conditional_R_max_diag_dev":
        cond_diag_dev,

    "analysis_status":
        "EXPLORATORY_EXTERNAL_EAS504_SUMMARY_CONDITIONING",
}


with (
    OUT /
    "L003_RS671_CONDITIONING_SUMMARY.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            summary.keys()
        ),
    )

    w.writeheader()
    w.writerow(
        summary
    )


print()
print(
    "===== CONDITIONING SUMMARY ====="
)

for k2, v2 in summary.items():
    print(
        f"{k2}={v2}"
    )


print()
print(
    "===== TOP 20 CONDITIONAL ====="
)

for x in ranked[:20]:

    print(
        x["variant_id"],
        x["rsid"],
        "raw_p=",
        x["p_raw"],
        "r2=",
        x["r2_to_rs671"],
        "z_cond=",
        x["z_cond_rs671"],
        "p_cond=",
        x["p_cond_rs671"],
    )

PY

}


# =============================================================================
# 03. CONDITIONAL SuSiE — N TOTAL / N EFF SENSITIVITY
# =============================================================================

run_conditional_susie() {

Rscript - \
  "$OUT" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

OUT <- args[[1]]

input_file <- file.path(
  OUT,
  "L003_RS671_CONDITIONAL_INPUT.tsv"
)

matrix_file <- file.path(
  OUT,
  "L003_RS671_CONDITIONAL_R.float64.bin"
)

if (!file.exists(input_file)) {
  stop("conditional input missing")
}

if (!file.exists(matrix_file)) {
  stop("conditional matrix missing")
}


D <- read.delim(
  input_file,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

nvar <- nrow(D)

cat(
  "N_VARIANTS=",
  nvar,
  "\n",
  sep=""
)


con <- file(
  matrix_file,
  "rb"
)

vec <- readBin(
  con,
  what="double",
  n=nvar*nvar,
  size=8,
  endian="little"
)

close(con)


if (length(vec) != nvar*nvar) {

  stop(
    "conditional R binary size mismatch"
  )

}


R <- matrix(
  vec,
  nrow=nvar,
  ncol=nvar,
  byrow=TRUE
)


cat(
  "R_MAX_ASYM=",
  max(
    abs(
      R - t(R)
    )
  ),
  "\n",
  sep=""
)


cat(
  "R_MAX_DIAG_DEV=",
  max(
    abs(
      diag(R) - 1
    )
  ),
  "\n",
  sep=""
)


if (!requireNamespace(
      "susieR",
      quietly=TRUE
    )) {
  stop("susieR missing")
}


z <- as.numeric(
  D$z_cond_rs671
)


N_CASE <- 22664
N_CONTROL <- 152022
N_TOTAL <- N_CASE + N_CONTROL

N_EFF <- 4 / (
  1 / N_CASE
  +
  1 / N_CONTROL
)


cat(
  "N_TOTAL=",
  N_TOTAL,
  "\n",
  sep=""
)

cat(
  "N_EFF=",
  N_EFF,
  "\n",
  sep=""
)


run_one <- function(
  mode,
  n_value
) {

  cat(
    "\n========================================\n"
  )

  cat(
    "RUN=",
    mode,
    "\n",
    sep=""
  )

  fn <- susieR::susie_rss

  ff <- names(
    formals(fn)
  )

  A <- list(
    z=z,
    R=R,
    n=n_value,
    L=10
  )

  if (
    "estimate_residual_variance"
    %in% ff
  ) {

    A$estimate_residual_variance <- FALSE

  }


  if (
    "max_iter"
    %in% ff
  ) {

    A$max_iter <- 2000

  }


  fit <- do.call(
    fn,
    A
  )


  saveRDS(
    fit,
    file.path(
      OUT,
      paste0(
        "L003_RS671_CONDITIONAL_SUSIE_",
        mode,
        ".rds"
      )
    )
  )


  pip <- as.numeric(
    fit$pip
  )

  top <- which.max(
    pip
  )


  cs <- list()

  if (
    !is.null(fit$sets)
    &&
    !is.null(fit$sets$cs)
  ) {

    cs <- fit$sets$cs

  }


  cat(
    "CONVERGED=",
    fit$converged,
    "\n",
    sep=""
  )

  cat(
    "N_CS=",
    length(cs),
    "\n",
    sep=""
  )

  cat(
    "TOP_VARIANT=",
    D$variant_id[top],
    "\n",
    sep=""
  )

  cat(
    "TOP_PIP=",
    pip[top],
    "\n",
    sep=""
  )


  summary <- data.frame(
    mode=mode,
    n_used=n_value,
    converged=fit$converged,
    n_variants=nvar,
    n_cs=length(cs),
    top_variant=D$variant_id[top],
    top_rsid=D$rsid[top],
    top_pip=pip[top],
    top_conditional_p=D$p_cond_rs671[top],
    stringsAsFactors=FALSE
  )


  cs_rows <- list()

  if (
    length(cs) > 0
  ) {

    for (
      i in seq_along(cs)
    ) {

      members <- cs[[i]]

      for (
        j in members
      ) {

        cs_rows[[
          length(cs_rows) + 1
        ]] <- data.frame(
          mode=mode,
          cs=i,
          variant_id=
            D$variant_id[j],
          rsid=
            D$rsid[j],
          pip=
            pip[j],
          conditional_p=
            D$p_cond_rs671[j],
          r2_to_rs671=
            D$r2_to_rs671[j],
          stringsAsFactors=FALSE
        )

      }

    }

  }


  list(
    summary=summary,
    cs=cs_rows,
    pip=pip
  )

}


A <- run_one(
  "N_TOTAL",
  N_TOTAL
)

B <- run_one(
  "N_EFF",
  N_EFF
)


SUM <- rbind(
  A$summary,
  B$summary
)


write.table(
  SUM,
  file.path(
    OUT,
    "L003_RS671_CONDITIONAL_SUSIE_SUMMARY.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


all_cs <- c(
  A$cs,
  B$cs
)


if (
  length(all_cs)
  > 0
) {

  CS <- do.call(
    rbind,
    all_cs
  )

} else {

  CS <- data.frame(
    mode=character(),
    cs=integer(),
    variant_id=character(),
    rsid=character(),
    pip=numeric(),
    conditional_p=numeric(),
    r2_to_rs671=numeric()
  )

}


write.table(
  CS,
  file.path(
    OUT,
    "L003_RS671_CONDITIONAL_SUSIE_CS.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


comparison <- data.frame(
  metric=c(
    "top_same",
    "pip_correlation",
    "max_abs_pip_difference"
  ),
  value=c(
    as.character(
      A$summary$top_variant
      ==
      B$summary$top_variant
    ),
    as.character(
      cor(
        A$pip,
        B$pip,
        use="complete.obs"
      )
    ),
    as.character(
      max(
        abs(
          A$pip
          -
          B$pip
        )
      )
    )
  )
)


write.table(
  comparison,
  file.path(
    OUT,
    "L003_RS671_N_SENSITIVITY.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


cat(
  "\n===== SUMMARY =====\n"
)

print(
  SUM
)

cat(
  "\n===== N SENSITIVITY =====\n"
)

print(
  comparison
)

RS

}


# =============================================================================
# 04. LOOK UP CONDITIONAL SECONDARY VARIANTS IN GIGASTROKE
# =============================================================================

lookup_gigastroke() {

python3 - \
  "$OUT" \
  "$GIGA" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])
GIGA = Path(sys.argv[2])

ranked_file = (
    OUT /
    "L003_RS671_CONDITIONAL_VARIANTS_RANKED.tsv"
)


if not ranked_file.exists():

    raise RuntimeError(
        "ranked conditional file missing"
    )


with ranked_file.open() as f:

    ranked = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


top = ranked[:50]

rank_by_variant = {
    x["variant_id"]: i + 1
    for i, x in enumerate(top)
}

cond_by_variant = {
    x["variant_id"]: x
    for x in top
}


out_rows = []


if GIGA.exists():

    with GIGA.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t",
        )

        for r in reader:

            if (
                r.get("locus")
                != "BBJ_IS_L003"
            ):
                continue

            vid = r.get(
                "variant_id"
            )

            if vid not in rank_by_variant:
                continue

            c = cond_by_variant[vid]

            out_rows.append({
                "conditional_rank":
                    rank_by_variant[vid],

                "variant_id":
                    vid,

                "rsid":
                    c.get(
                        "rsid",
                        ""
                    ),

                "bbj_conditional_p":
                    c.get(
                        "p_cond_rs671",
                        ""
                    ),

                "r2_to_rs671":
                    c.get(
                        "r2_to_rs671",
                        ""
                    ),

                "phenotype":
                    r.get(
                        "phenotype",
                        ""
                    ),

                "giga_p":
                    r.get(
                        "giga_p",
                        ""
                    ),

                "giga_beta_aligned":
                    r.get(
                        "giga_beta_aligned",
                        ""
                    ),

                "harmonization":
                    r.get(
                        "harmonization",
                        ""
                    ),

                "direction_concordant":
                    r.get(
                        "direction_concordant",
                        ""
                    ),
            })


out_rows.sort(
    key=lambda x: (
        int(
            x[
                "conditional_rank"
            ]
        ),
        x["phenotype"],
    )
)


outfile = (
    OUT /
    "L003_SECONDARY_GIGASTROKE_LOOKUP.tsv"
)


fields = [
    "conditional_rank",
    "variant_id",
    "rsid",
    "bbj_conditional_p",
    "r2_to_rs671",
    "phenotype",
    "giga_p",
    "giga_beta_aligned",
    "harmonization",
    "direction_concordant",
]


with outfile.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(
        out_rows
    )


print(
    "GIGASTROKE_LOOKUP_ROWS=",
    len(out_rows)
)


for x in out_rows[:50]:

    print(
        x["conditional_rank"],
        x["variant_id"],
        x["phenotype"],
        "cond_p=",
        x["bbj_conditional_p"],
        "giga_p=",
        x["giga_p"],
        "beta=",
        x["giga_beta_aligned"],
    )

PY

}


# =============================================================================
# 05. BUILD READINESS
# =============================================================================

build_readiness() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


def read_one(name):

    p = OUT / name

    if not p.exists():
        return None

    with p.open() as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )

    return (
        rows[0]
        if rows
        else None
    )


conditioning = read_one(
    "L003_RS671_CONDITIONING_SUMMARY.tsv"
)


susie_file = (
    OUT /
    "L003_RS671_CONDITIONAL_SUSIE_SUMMARY.tsv"
)

susie_rows = []

if susie_file.exists():

    with susie_file.open() as f:

        susie_rows = list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


condition_ok = (
    conditioning is not None
)


susie_ok = (
    len(susie_rows) == 2
    and all(
        str(x.get("converged", ""))
        .upper()
        in {
            "TRUE",
            "T",
            "1",
        }
        for x in susie_rows
    )
)


secondary_gws = False

if conditioning is not None:

    try:

        secondary_gws = (
            int(
                conditioning[
                    "conditional_p_lt_5e8_n"
                ]
            )
            > 0
        )

    except Exception:
        pass


rows = [

    (
        "PHASE8E_L003_RS671",
        (
            "COMPUTATIONAL_COMPLETE"
            if condition_ok
            else "INCOMPLETE"
        )
    ),

    (
        "RS671_PRESERVED_AS_CONDITION_VARIANT",
        (
            "YES"
            if condition_ok
            else "UNKNOWN"
        )
    ),

    (
        "CONDITIONING_METHOD",
        "APPROX_SUMMARY_CONDITIONAL"
    ),

    (
        "LD_REFERENCE",
        "1000G_EAS504_EXTERNAL_PROXY"
    ),

    (
        "BBJ_REFERENCE_ALLELES",
        "GRCH37_VERIFIED"
    ),

    (
        "CONDITIONAL_SUSIE",
        (
            "COMPLETE_N_TOTAL_AND_N_EFF_SENSITIVITY"
            if susie_ok
            else "PENDING_OR_REVIEW"
        )
    ),

    (
        "SECONDARY_GENOMEWIDE_SIGNAL_AFTER_RS671",
        (
            "PRESENT"
            if secondary_gws
            else "NOT_DETECTED_AT_P_LT_5E8"
        )
    ),

    (
        "L003_PAPER_GRADE",
        "NO_EXTERNAL_EAS504_LD_AND_SUMMARY_CONDITIONING"
    ),

    (
        "L003_NEXT",
        (
            "SECONDARY_SIGNAL_CHARACTERISATION"
            if secondary_gws
            else "RS671_PRIMARY_SIGNAL_INTERPRETATION"
        )
    ),

    (
        "ALDH2_EQLT_COLOC",
        "DO_NOT_REINTERPRET_EUR_INTERSECTION_AS_NEGATIVE"
    ),

    (
        "PHASE11B_L003",
        "RETAIN_AS_EXPLORATORY_INCONCLUSIVE"
    ),

]


with (
    OUT /
    "PHASE8E_READINESS.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
    )

    w.writerow([
        "component",
        "status",
    ])

    w.writerows(rows)


for a, b in rows:
    print(
        f"{a} = {b}"
    )

PY

}


# =============================================================================
# RUN
# =============================================================================

run_step \
  "01_preflight" \
  preflight

run_step \
  "02_build_rs671_conditioned_input" \
  build_conditioned_input

run_step \
  "03_conditional_susie" \
  run_conditional_susie

run_step \
  "04_gigastroke_secondary_lookup" \
  lookup_gigastroke

run_step \
  "05_build_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE8E L003 / rs671 COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== CONDITIONING SUMMARY ====="

cat \
  "$OUT/L003_RS671_CONDITIONING_SUMMARY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== CONDITIONAL SuSiE ====="

cat \
  "$OUT/L003_RS671_CONDITIONAL_SUSIE_SUMMARY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== N SENSITIVITY ====="

cat \
  "$OUT/L003_RS671_N_SENSITIVITY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== TOP 20 CONDITIONAL VARIANTS ====="

head -21 \
  "$OUT/L003_RS671_CONDITIONAL_VARIANTS_RANKED.tsv" \
  2>/dev/null \
  || true


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE8E_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE8E_READINESS.tsv"


echo
echo "===== OUTPUT FILES ====="

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -printf '%f\t%s bytes\n' \
  | sort


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
