#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

BASE="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit"
OUT="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit_r1"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_eas_phase8d_r1/$RUN_ID"
STATUS="$OUT/STEP_STATUS_R1.tsv"

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
# 1. susieR API AUDIT REPAIR
# =============================================================================

audit_susier_api() {

Rscript - "$OUT" <<'RS'

args <- commandArgs(trailingOnly=TRUE)
OUT <- args[[1]]

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

F <- file.path(
  OUT,
  "SUSIER_API_AUDIT_R1.txt"
)

STATUS <- file.path(
  OUT,
  "SUSIER_API_STATUS_R1.txt"
)

con <- file(
  F,
  open="wt"
)

sink(con)
sink(con, type="message")

api_status <- "AUDITED"

tryCatch({

  cat("===== susieR API AUDIT R1 =====\n\n")

  if (!requireNamespace(
        "susieR",
        quietly=TRUE
      )) {
    stop("susieR is not installed")
  }

  cat(
    "susieR_version=",
    as.character(
      packageVersion("susieR")
    ),
    "\n\n",
    sep=""
  )

  targets <- c(
    "susie_rss",
    "estimate_s_rss",
    "kriging_rss"
  )

  for (fn in targets) {

    cat(
      "\n============================================================\n"
    )

    cat(
      "FUNCTION=",
      fn,
      "\n",
      sep=""
    )

    cat(
      "============================================================\n"
    )

    obj <- tryCatch(
      getExportedValue(
        "susieR",
        fn
      ),
      error=function(e) NULL
    )

    if (is.null(obj)) {

      cat(
        "STATUS=FUNCTION_NOT_EXPORTED\n"
      )

      api_status <<- "PENDING_API_AUDIT"
      next
    }

    cat(
      "\n===== FORMALS =====\n"
    )

    print(
      formals(obj)
    )

    cat(
      "\n===== BODY PREVIEW =====\n"
    )

    body_text <- tryCatch(
      deparse(
        body(obj)
      ),
      error=function(e) character()
    )

    if (length(body_text)) {

      cat(
        paste(
          head(
            body_text,
            120
          ),
          collapse="\n"
        ),
        "\n"
      )

    } else {

      cat(
        "BODY_NOT_AVAILABLE\n"
      )

    }

    cat(
      "\n===== GETANYWHERE =====\n"
    )

    ga <- tryCatch(
      getAnywhere(fn),
      error=function(e) NULL
    )

    if (!is.null(ga)) {
      print(ga)
    }

  }

}, error=function(e) {

  api_status <<- "PENDING_API_AUDIT"

  cat(
    "\nAPI_AUDIT_ERROR=",
    conditionMessage(e),
    "\n",
    sep=""
  )

})

sink(type="message")
sink()

close(con)

writeLines(
  api_status,
  STATUS
)

cat(
  "SUSIER_API_STATUS=",
  api_status,
  "\n",
  sep=""
)

RS

}


# =============================================================================
# 2. LD STRUCTURAL QC REPAIR
# =============================================================================

audit_ld_structural() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import math
import sys

import numpy as np


OUT = Path(sys.argv[1])
ROOT = Path("/srv/is-analysis")

OUT.mkdir(
    parents=True,
    exist_ok=True
)


SEARCH_ROOTS = [
    ROOT / "results/is/stage3_finemap",
    ROOT / "results/is/stage4_cross_eas",
    ROOT / "data/is/ld_reference",
]


def relevant_matrix(p: Path) -> bool:

    low = str(p).lower()

    if not low.endswith(".bin"):
        return False

    return any(
        x in low
        for x in [
            "bbj_is_l001",
            "bbj_is_l002",
            "bbj_is_l003",
            "bbj_is_l004",
        ]
    )


def vars_candidates(p: Path):

    s = str(p)

    candidates = [
        Path(s + ".vars"),
    ]

    if s.endswith(".vcor1.bin"):

        candidates += [
            Path(
                s.replace(
                    ".vcor1.bin",
                    ".vcor1.bin.vars"
                )
            ),
            Path(
                s.replace(
                    ".vcor1.bin",
                    ".vcor1.vars"
                )
            ),
        ]

    if s.endswith(".vcor.bin"):

        candidates += [
            Path(
                s.replace(
                    ".vcor.bin",
                    ".vcor.bin.vars"
                )
            ),
            Path(
                s.replace(
                    ".vcor.bin",
                    ".vcor.vars"
                )
            ),
        ]

    return list(
        dict.fromkeys(candidates)
    )


bins = []

for root in SEARCH_ROOTS:

    if not root.exists():
        continue

    for p in root.rglob("*.bin"):

        if relevant_matrix(p):
            bins.append(p)

bins = sorted(
    set(bins)
)


rows = []


for p in bins:

    vp = next(
        (
            x
            for x in vars_candidates(p)
            if x.exists()
        ),
        None
    )

    row = {
        "matrix_path": str(p),
        "vars_path": "",
        "n_vars_lines": "",
        "n": "",
        "actual_bytes": p.stat().st_size,
        "expected_f64": "",
        "expected_f32": "",
        "dtype": "",
        "format_status": "",
        "max_asym": "",
        "diag_min": "",
        "diag_max": "",
        "max_diag_absdev1": "",
        "matrix_min": "",
        "matrix_max": "",
        "nonfinite": "",
        "structural_status": "",
    }

    if vp is None:

        row["format_status"] = (
            "VARS_NOT_FOUND"
        )

        row["structural_status"] = (
            "UNRESOLVED"
        )

        rows.append(row)
        continue

    row["vars_path"] = str(vp)

    with vp.open(
        errors="replace"
    ) as f:

        lines = [
            x.strip()
            for x in f
            if x.strip()
        ]

    row["n_vars_lines"] = len(lines)

    possible_n = [
        len(lines)
    ]

    if lines:

        first = lines[0].upper()

        if (
            first.startswith("#")
            or first in {
                "ID",
                "VARIANT",
                "VARIANT_ID",
            }
        ):
            possible_n.append(
                len(lines) - 1
            )

    possible_n = [
        n
        for n in dict.fromkeys(
            possible_n
        )
        if n > 0
    ]

    actual = p.stat().st_size

    resolved = None

    for n in possible_n:

        if actual == n * n * 8:

            resolved = (
                n,
                np.dtype("<f8"),
                "FLOAT64_EXACT"
            )

            break

        if actual == n * n * 4:

            resolved = (
                n,
                np.dtype("<f4"),
                "FLOAT32_EXACT"
            )

            break

    if resolved is None:

        row["format_status"] = (
            "FORMAT_UNRESOLVED"
        )

        row["structural_status"] = (
            "UNRESOLVED"
        )

        rows.append(row)
        continue

    n, dtype, fmt = resolved

    row["n"] = n
    row["dtype"] = str(dtype)
    row["format_status"] = fmt
    row["expected_f64"] = n * n * 8
    row["expected_f32"] = n * n * 4

    M = np.memmap(
        p,
        dtype=dtype,
        mode="r",
        shape=(n, n)
    )

    diag = np.asarray(
        M[
            np.arange(n),
            np.arange(n)
        ],
        dtype=np.float64
    )

    row["diag_min"] = float(
        np.nanmin(diag)
    )

    row["diag_max"] = float(
        np.nanmax(diag)
    )

    row["max_diag_absdev1"] = float(
        np.nanmax(
            np.abs(
                diag - 1.0
            )
        )
    )

    matrix_min = math.inf
    matrix_max = -math.inf
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
            M[i0:i1, :],
            dtype=np.float64
        )

        finite = np.isfinite(A)

        nonfinite += int(
            A.size
            - finite.sum()
        )

        if finite.any():

            matrix_min = min(
                matrix_min,
                float(
                    np.nanmin(A)
                )
            )

            matrix_max = max(
                matrix_max,
                float(
                    np.nanmax(A)
                )
            )

        B = np.asarray(
            M[:, i0:i1],
            dtype=np.float64
        ).T

        D = np.abs(
            A - B
        )

        if np.isfinite(D).any():

            max_asym = max(
                max_asym,
                float(
                    np.nanmax(D)
                )
            )

    row["matrix_min"] = matrix_min
    row["matrix_max"] = matrix_max
    row["max_asym"] = max_asym
    row["nonfinite"] = nonfinite

    tol = (
        1e-5
        if dtype.itemsize == 4
        else 1e-10
    )

    structural_pass = (
        nonfinite == 0
        and max_asym <= tol
        and row["max_diag_absdev1"] <= tol
        and matrix_min >= -1.00001
        and matrix_max <= 1.00001
    )

    row["structural_status"] = (
        "PASS"
        if structural_pass
        else "REVIEW"
    )

    rows.append(row)


outfile = (
    OUT /
    "LD_STRUCTURAL_QC_R1.tsv"
)

fields = [
    "matrix_path",
    "vars_path",
    "n_vars_lines",
    "n",
    "actual_bytes",
    "expected_f64",
    "expected_f32",
    "dtype",
    "format_status",
    "max_asym",
    "diag_min",
    "diag_max",
    "max_diag_absdev1",
    "matrix_min",
    "matrix_max",
    "nonfinite",
    "structural_status",
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
    w.writerows(rows)


print(
    "LD_MATRICES=",
    len(rows)
)

for r in rows:

    print(
        r["structural_status"],
        r["format_status"],
        "n=",
        r["n"],
        "asym=",
        r["max_asym"],
        "diagdev=",
        r["max_diag_absdev1"],
        r["matrix_path"]
    )

PY

}


# =============================================================================
# 3. EXACT SUBSET / PHASE8B PROVENANCE AUDIT
# =============================================================================

audit_exact_subset() {

python3 - "$OUT" "$REPO" <<'PY'

from pathlib import Path
import csv
import hashlib
import re
import sys


OUT = Path(sys.argv[1])
REPO = Path(sys.argv[2])
ROOT = Path("/srv/is-analysis")

script = (
    REPO /
    "server/run_is_eas_master_phase8b_repair.sh"
)

rows = []

report = (
    OUT /
    "PHASE8B_SUBSET_SCRIPT_AUDIT.txt"
)


if not script.exists():

    rows.append({
        "target": "PHASE8B_SCRIPT",
        "script": str(script),
        "candidate_path": "",
        "status": "SCRIPT_NOT_FOUND",
        "evidence": "",
    })

    report.write_text(
        "PHASE8B SCRIPT NOT FOUND\n"
    )

else:

    text = script.read_text(
        errors="replace"
    )

    lines = text.splitlines()

    sha = hashlib.sha256(
        script.read_bytes()
    ).hexdigest()

    terms = [
        "subset",
        "common",
        "match",
        "idx",
        "index",
        "ld[",
        "r[",
        "matrix",
        "l003",
        "l004",
        "possible_switch",
        "kriging",
        "estimate_s_rss",
        "susie_rss",
        "writebin",
        "save",
        "saverds",
        "write.table",
        "fwrite",
        "np.save",
        "tofile",
    ]

    evidence = []

    for i, line in enumerate(
        lines,
        start=1
    ):

        low = line.lower()

        if any(
            t in low
            for t in terms
        ):

            lo = max(
                1,
                i - 2
            )

            hi = min(
                len(lines),
                i + 2
            )

            block = [
                f"{j}: {lines[j-1]}"
                for j in range(
                    lo,
                    hi + 1
                )
            ]

            evidence.append(
                "\n".join(block)
            )

    explicit_write_terms = [
        "writebin(",
        "saverds(",
        "np.save(",
        ".tofile(",
        "write.table(",
        "fwrite(",
    ]

    explicit_matrix_write = any(
        x in text.lower()
        for x in explicit_write_terms
    )

    report.write_text(
        "SCRIPT="
        + str(script)
        + "\n"
        + "SHA256="
        + sha
        + "\n"
        + "EXPLICIT_MATRIX_WRITE_TOKEN="
        + str(
            explicit_matrix_write
        )
        + "\n\n"
        + "\n\n"
          "============================================================\n"
          "\n\n".join(evidence)
        + "\n"
    )

    # ---------------------------------------------------------
    # Candidate persisted matrices
    # ---------------------------------------------------------

    targets = {
        "L003_BBJ": [],
        "L004_AIS": [],
    }

    search_roots = [
        ROOT /
        "results/is/stage4_cross_eas",
        ROOT /
        "results/is/stage3_finemap",
    ]

    for base in search_roots:

        if not base.exists():
            continue

        for p in base.rglob("*"):

            if not p.is_file():
                continue

            low = str(p).lower()

            if p.suffix.lower() not in {
                ".bin",
                ".npy",
                ".rds",
                ".rdata",
                ".mtx",
            }:
                continue

            if (
                "l003" in low
                and "bbj" in low
            ):
                targets[
                    "L003_BBJ"
                ].append(p)

            if (
                "l004" in low
                and "ais" in low
            ):
                targets[
                    "L004_AIS"
                ].append(p)

    for target, paths in targets.items():

        unique = sorted(
            set(paths)
        )

        if not unique:

            rows.append({
                "target": target,
                "script": str(script),
                "candidate_path": "",
                "status":
                    "NO_PERSISTED_CANDIDATE_FOUND",
                "evidence":
                    "No exact subset matrix candidate found by filename/path.",
            })

            continue

        for p in unique:

            rows.append({
                "target": target,
                "script": str(script),
                "candidate_path": str(p),
                "status":
                    "CANDIDATE_FOUND_PROVENANCE_NOT_ASSUMED",
                "evidence":
                    "Candidate exists; exact Phase8B diagnostic provenance must be demonstrated from script/path lineage.",
            })


outfile = (
    OUT /
    "EXACT_SUBSET_PROVENANCE_R1.tsv"
)

fields = [
    "target",
    "script",
    "candidate_path",
    "status",
    "evidence",
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
    w.writerows(rows)


print(
    "SUBSET_PROVENANCE_ROWS=",
    len(rows)
)

for r in rows:

    print(
        r["target"],
        r["status"],
        r["candidate_path"]
    )

print(
    "SCRIPT_AUDIT=",
    report
)

PY

}


# =============================================================================
# 4. REFRESH READINESS
# =============================================================================

refresh_readiness() {

python3 - "$BASE" "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


BASE = Path(sys.argv[1])
OUT = Path(sys.argv[2])


# ------------------------------------------------------------
# susieR API
# ------------------------------------------------------------

api = "PENDING_API_AUDIT"

p = OUT / "SUSIER_API_STATUS_R1.txt"

if p.exists():

    x = p.read_text().strip()

    if x == "AUDITED":
        api = "AUDITED"


# ------------------------------------------------------------
# LD structural result
# ------------------------------------------------------------

ld_status = "PENDING"

ldf = OUT / "LD_STRUCTURAL_QC_R1.tsv"

if ldf.exists():

    with ldf.open() as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    tested = [
        x
        for x in rows
        if x[
            "format_status"
        ] in {
            "FLOAT64_EXACT",
            "FLOAT32_EXACT",
        }
    ]

    failures = [
        x
        for x in tested
        if x[
            "structural_status"
        ] != "PASS"
    ]

    unresolved = [
        x
        for x in rows
        if x[
            "structural_status"
        ] == "UNRESOLVED"
    ]

    if (
        tested
        and not failures
        and not unresolved
    ):

        ld_status = (
            "PASS_ALL_DISCOVERED_MATRICES"
        )

    elif (
        tested
        and not failures
    ):

        ld_status = (
            "PASS_TESTED_WITH_UNRESOLVED_CANDIDATES"
        )

    else:

        ld_status = (
            "REVIEW_REQUIRED"
        )


# ------------------------------------------------------------
# Harmonization
# ------------------------------------------------------------

harmon_status = "PENDING"

hf = (
    BASE /
    "HARMONIZED_INPUT_AUDIT.tsv"
)

if hf.exists():

    with hf.open() as f:

        hrows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    complement = 0

    alignment_files = 0

    for x in hrows:

        if x.get(
            "alignment_column"
        ):

            alignment_files += 1

        try:
            complement += int(
                x.get(
                    "contains_complement_status",
                    0
                )
            )
        except Exception:
            pass

    if (
        alignment_files > 0
        and complement == 0
    ):

        harmon_status = (
            "FINAL_AUDITED_FILES_MATCH_SWAP_ONLY_NO_COMPLEMENT"
        )

    else:

        harmon_status = (
            "REVIEW_REQUIRED"
        )


# ------------------------------------------------------------
# Exact subset provenance
# ------------------------------------------------------------

subset_status = (
    "PENDING_EXACT_PHASE8B_PROVENANCE"
)

sp = (
    OUT /
    "EXACT_SUBSET_PROVENANCE_R1.tsv"
)

if sp.exists():

    with sp.open() as f:

        srows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    if srows:

        subset_status = (
            "CANDIDATES_INVENTORIED_PROVENANCE_NOT_ASSUMED"
        )


# ------------------------------------------------------------
# Current readiness
# ------------------------------------------------------------

rows = [

    (
        "PHASE8D_R1",
        "AUDIT_COMPLETE"
    ),

    (
        "GIGASTROKE_ANCESTRY",
        "VERIFIED_LOCAL_METADATA_EAST_ASIAN"
    ),

    (
        "GIGASTROKE_SAMPLE_SIZE",
        "VERIFIED_LOCAL_METADATA_TOTAL_N"
    ),

    (
        "GIGASTROKE_CASE_CONTROL_N",
        "NOT_AVAILABLE_IN_CURRENT_LOCAL_METADATA"
    ),

    (
        "GIGASTROKE_REF_ALT",
        "PENDING_GRCH37_REFERENCE_AUDIT"
    ),

    (
        "BBJ_REF_FASTA",
        "PENDING_GRCH37_FASTA"
    ),

    (
        "HARMONIZATION",
        harmon_status
    ),

    (
        "LD_ORIENTATION",
        "DOCUMENTED_SIGNED_REF_BASED_FROM_PHASE5_PROVENANCE"
    ),

    (
        "LD_STRUCTURAL_QC",
        ld_status
    ),

    (
        "SUSIER_RSS_API",
        api
    ),

    (
        "EXACT_SUBSET_MATRIX_PROVENANCE",
        subset_status
    ),

    (
        "BBJ_PRIMARY_FINEMAP",
        "EXPLORATORY_EAS504"
    ),

    (
        "PHASE7B_GIGASTROKE_SUSIE",
        "EXPLORATORY_SAMPLE_SIZE_MISSING"
    ),

    (
        "PHASE8_GIGASTROKE_NAWARE",
        "INVALID_N_ASSIGNMENT"
    ),

    (
        "PHASE8B_GIGASTROKE_NAWARE",
        "EXPLORATORY_CORRECT_N_BUT_VALIDATION_PENDING"
    ),

    (
        "PHASE8C_MASTER_REBUILD",
        "BOOKKEEPING_COMPLETE"
    ),

    (
        "BBJ_GIGASTROKE_CS_COMPARE",
        "EXPLORATORY_REF_PENDING"
    ),

    (
        "PHASE9A_OUTPUTS",
        "EXPLORATORY_FROZEN_CORE_VALIDATION_PENDING"
    ),

    (
        "PHASE9B_OUTPUTS",
        "EXPLORATORY_FROZEN_CORE_VALIDATION_PENDING"
    ),

    (
        "PHASE10_11_OUTPUTS",
        "EXPLORATORY_FROZEN_CORE_VALIDATION_PENDING"
    ),

    (
        "PHASE11B_FGF5",
        "EXPLORATORY_MOLECULAR_COLOC_SUPPORT"
    ),

    (
        "L003_CONDITIONING",
        "PENDING_EAS_SPECIFIC_RS671_PRESERVING_ANALYSIS"
    ),

    (
        "TPMI",
        "WAITING_RATE_LIMIT"
    ),

    (
        "CKB",
        "WAITING_DECRYPTION_KEY"
    ),
]


outfile = (
    OUT /
    "SCIENTIFIC_READINESS_R1.tsv"
)

with outfile.open(
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
  "01_susier_api_repair" \
  audit_susier_api

run_step \
  "02_ld_structural_qc_repair" \
  audit_ld_structural

run_step \
  "03_exact_subset_provenance" \
  audit_exact_subset

run_step \
  "04_refresh_readiness" \
  refresh_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE8D-R1 COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"

echo
echo "===== READINESS R1 ====="

column -t -s $'\t' \
  "$OUT/SCIENTIFIC_READINESS_R1.tsv" \
  2>/dev/null \
  || cat "$OUT/SCIENTIFIC_READINESS_R1.tsv"

echo
echo "===== LD QC SUMMARY ====="

if [[ -s "$OUT/LD_STRUCTURAL_QC_R1.tsv" ]]; then

  python3 - "$OUT/LD_STRUCTURAL_QC_R1.tsv" <<'PY'

import csv
import sys

p = sys.argv[1]

with open(p) as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

print(
    "matrices =",
    len(rows)
)

for x in rows:

    print(
        x["structural_status"],
        x["format_status"],
        "n=" + str(x["n"]),
        "asym=" + str(x["max_asym"]),
        "diagdev=" + str(
            x["max_diag_absdev1"]
        ),
        x["matrix_path"]
    )

PY

fi


echo
echo "===== EXACT SUBSET PROVENANCE ====="

if [[ -s "$OUT/EXACT_SUBSET_PROVENANCE_R1.tsv" ]]; then
  cat "$OUT/EXACT_SUBSET_PROVENANCE_R1.tsv"
fi


echo
echo "===== OUTPUTS ====="

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -printf '%f\t%s bytes\n' \
  | sort

echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"

