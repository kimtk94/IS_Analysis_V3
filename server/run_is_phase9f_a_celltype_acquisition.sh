#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

P9E="$ROOT/results/is/stage5_functional/phase9e_celltype_resource_audit"

DATA="$ROOT/data/is/celltype"
OUT="$ROOT/results/is/stage5_functional/phase9f_a_acquisition"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9f_a/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

ENV="$ROOT/envs/is_singlecell"

mkdir -p \
  "$DATA" \
  "$OUT" \
  "$LOGDIR" \
  "$ROOT/envs"

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

  echo "===== PHASE9E ====="

  if [[ ! -s "$P9E/STEP_STATUS.tsv" ]]; then
    echo "PHASE9E_STATUS_MISSING"
    return 1
  fi

  cat "$P9E/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P9E/STEP_STATUS.tsv"
  )"

  echo
  echo "PHASE9E_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    return 1
  fi


  echo
  echo "===== DEPENDENCIES ====="

  for exe in \
    python3 \
    curl \
    tar \
    gzip \
    sha256sum
  do

    if command -v "$exe" >/dev/null 2>&1; then
      echo "FOUND $exe $(command -v "$exe")"
    else
      echo "MISSING $exe"
      return 1
    fi

  done


  echo
  echo "===== DISK ====="

  df -h "$ROOT"

  AVAIL_KB="$(
    df -Pk "$ROOT" \
    | awk 'NR==2 {print $4}'
  )"

  echo "AVAILABLE_KB=$AVAIL_KB"

  # Require >= 20 GB free.
  if [[ "$AVAIL_KB" -lt 20971520 ]]; then
    echo "INSUFFICIENT_DISK"
    return 1
  fi

  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. PYTHON SINGLE-CELL ENVIRONMENT
# =============================================================================

prepare_environment() {

  echo "===== ENV ====="
  echo "ENV=$ENV"

  if [[ ! -x "$ENV/bin/python" ]]; then

    echo "Creating venv..."

    python3 -m venv "$ENV"

    RC=$?

    if [[ $RC -ne 0 ]]; then

      echo "VENV_CREATION_FAIL"

      echo "Trying python3-venv install..."

      sudo apt-get update
      sudo apt-get install -y python3-venv

      python3 -m venv "$ENV"

      RC=$?

      if [[ $RC -ne 0 ]]; then
        return 1
      fi

    fi

  fi


  "$ENV/bin/python" \
    -m pip install \
    --upgrade \
    pip \
    setuptools \
    wheel

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  echo
  echo "===== INSTALL SCIENTIFIC STACK ====="

  "$ENV/bin/python" \
    -m pip install \
    numpy \
    pandas \
    scipy \
    h5py \
    anndata \
    scanpy \
    pyarrow

  RC=$?

  if [[ $RC -ne 0 ]]; then
    echo "PYTHON_PACKAGE_INSTALL_FAIL"
    return 1
  fi


  echo
  echo "===== VERSION AUDIT ====="

  "$ENV/bin/python" - <<'PY'
import numpy
import pandas
import scipy
import h5py
import anndata
import scanpy

mods = [
    numpy,
    pandas,
    scipy,
    h5py,
    anndata,
    scanpy,
]

for m in mods:
    print(
        m.__name__,
        getattr(m, "__version__", "")
    )
PY


  "$ENV/bin/python" \
    -m pip freeze \
    > "$OUT/PYTHON_ENV_FREEZE.txt"

  echo "ENVIRONMENT=PASS"
}


# =============================================================================
# Generic resumable downloader
# =============================================================================

download_file() {

  local URL="$1"
  local DEST="$2"

  mkdir -p "$(dirname "$DEST")"

  if [[ -s "$DEST" ]]; then

    echo "EXISTS $DEST $(stat -c%s "$DEST")"

    return 0
  fi


  local PART="${DEST}.part"

  echo
  echo "URL=$URL"
  echo "DEST=$DEST"

  curl \
    -L \
    --fail \
    --retry 5 \
    --retry-delay 5 \
    --connect-timeout 30 \
    -C - \
    -o "$PART" \
    "$URL"

  local RC=$?

  if [[ $RC -ne 0 ]]; then

    echo "DOWNLOAD_FAIL=$RC"

    return 1
  fi


  if [[ ! -s "$PART" ]]; then

    echo "EMPTY_DOWNLOAD"

    return 1
  fi


  mv "$PART" "$DEST"

  echo "DOWNLOAD_PASS $(stat -c%s "$DEST")"

  return 0
}


# =============================================================================
# 03. GSE189432 — STROKE MYELOID
# =============================================================================

download_gse189432() {

  D="$DATA/GSE189432"
  mkdir -p "$D"

  BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189432/suppl"

  download_file \
    "$BASE/GSE189432_RAW.tar" \
    "$D/GSE189432_RAW.tar"

  RC1=$?

  download_file \
    "$BASE/GSE189432_annotations.csv.gz" \
    "$D/GSE189432_annotations.csv.gz"

  RC2=$?


  if [[ $RC1 -ne 0 || $RC2 -ne 0 ]]; then
    return 1
  fi


  echo
  echo "===== TAR VALIDATION ====="

  tar -tf \
    "$D/GSE189432_RAW.tar" \
    > "$D/ARCHIVE_CONTENTS.txt"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  echo "FILES_IN_TAR=$(wc -l < "$D/ARCHIVE_CONTENTS.txt")"

  head -30 \
    "$D/ARCHIVE_CONTENTS.txt"

  echo
  echo "===== EXTRACT ====="

  mkdir -p "$D/processed"

  tar -xf \
    "$D/GSE189432_RAW.tar" \
    -C "$D/processed"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  find "$D/processed" \
    -maxdepth 2 \
    -type f \
    -printf '%f\t%s\n' \
    | sort \
    > "$D/EXTRACTED_FILES.tsv"

  echo "GSE189432=PASS"
}


# =============================================================================
# 04. GSE234052 — STROKE PERICYTE
# =============================================================================

download_gse234052() {

  D="$DATA/GSE234052"
  mkdir -p "$D"

  BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234052/suppl"

  FILES=(
    "GSE234052_barcodes.tsv.gz"
    "GSE234052_features.tsv.gz"
    "GSE234052_matrix.mtx.gz"
  )

  FAIL=0

  for fn in "${FILES[@]}"
  do

    download_file \
      "$BASE/$fn" \
      "$D/$fn"

    RC=$?

    if [[ $RC -ne 0 ]]; then
      FAIL=$((FAIL + 1))
    fi

  done


  if [[ $FAIL -ne 0 ]]; then
    echo "DOWNLOAD_FAILURES=$FAIL"
    return 1
  fi


  echo
  echo "===== GZIP VALIDATION ====="

  for fn in "${FILES[@]}"
  do

    gzip -t "$D/$fn"

    RC=$?

    echo "$fn gzip_rc=$RC"

    if [[ $RC -ne 0 ]]; then
      return 1
    fi

  done


  echo "GSE234052=PASS"
}


# =============================================================================
# 05. GSE225948 — STROKE IMMUNE / ENDOTHELIAL
# =============================================================================

download_gse225948() {

  D="$DATA/GSE225948"
  mkdir -p "$D"

  BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225948/suppl"

  download_file \
    "$BASE/GSE225948_RAW.tar" \
    "$D/GSE225948_RAW.tar"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  tar -tf \
    "$D/GSE225948_RAW.tar" \
    > "$D/ARCHIVE_CONTENTS.txt"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  mkdir -p "$D/processed"

  tar -xf \
    "$D/GSE225948_RAW.tar" \
    -C "$D/processed"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  find "$D/processed" \
    -maxdepth 2 \
    -type f \
    -printf '%f\t%s\n' \
    | sort \
    > "$D/EXTRACTED_FILES.tsv"


  echo
  echo "FILES_IN_TAR=$(wc -l < "$D/ARCHIVE_CONTENTS.txt")"

  echo "GSE225948=PASS"
}


# =============================================================================
# 06. GSE256490 — HUMAN BRAIN VASCULAR REFERENCE
# =============================================================================

download_gse256490() {

  D="$DATA/GSE256490"
  mkdir -p "$D"

  BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256490/suppl"

  download_file \
    "$BASE/GSE256490_RAW.tar" \
    "$D/GSE256490_RAW.tar"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  echo
  echo "===== ARCHIVE VALIDATION ====="

  tar -tf \
    "$D/GSE256490_RAW.tar" \
    > "$D/ARCHIVE_CONTENTS.txt"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  echo "FILES_IN_TAR=$(wc -l < "$D/ARCHIVE_CONTENTS.txt")"

  head -40 \
    "$D/ARCHIVE_CONTENTS.txt"


  # Do not blindly extract if expanded footprint is unexpectedly huge.
  EST_KB="$(
    tar -tvf "$D/GSE256490_RAW.tar" \
    | awk '{s += $3} END {printf "%.0f\n", s/1024}'
  )"

  echo "ESTIMATED_EXTRACTED_KB=$EST_KB"


  AVAIL_KB="$(
    df -Pk "$D" \
    | awk 'NR==2 {print $4}'
  )"

  echo "AVAILABLE_KB=$AVAIL_KB"


  if [[ "$EST_KB" -gt "$((AVAIL_KB / 2))" ]]; then

    echo "EXTRACTION_SKIPPED_DISK_SAFETY"

    return 0
  fi


  mkdir -p "$D/processed"

  tar -xf \
    "$D/GSE256490_RAW.tar" \
    -C "$D/processed"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  find "$D/processed" \
    -maxdepth 2 \
    -type f \
    -printf '%f\t%s\n' \
    | sort \
    > "$D/EXTRACTED_FILES.tsv"


  echo "GSE256490=PASS"
}


# =============================================================================
# 07. UNIVERSAL FILE INVENTORY
# =============================================================================

inventory_downloads() {

python3 - "$DATA" "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


DATA = Path(sys.argv[1])
OUT = Path(sys.argv[2])


extensions = [

    ".h5ad",
    ".h5",
    ".loom",

    ".mtx",
    ".mtx.gz",

    ".csv",
    ".csv.gz",

    ".tsv",
    ".tsv.gz",

    ".txt",
    ".txt.gz",

    ".rds",
    ".rda",
    ".rdata",

    ".bed",
    ".bed.gz",

    ".parquet",
    ".parquet.gz",

    ".tar",
]


rows = []


for p in DATA.rglob("*"):

    if not p.is_file():
        continue

    low = p.name.lower()

    matched = next(
        (
            x for x in extensions
            if low.endswith(x)
        ),
        ""
    )

    if not matched:
        continue


    accession = ""

    for part in p.parts:

        if part.startswith("GSE"):
            accession = part
            break


    rows.append({

        "accession":
            accession,

        "path":
            str(p),

        "filename":
            p.name,

        "extension":
            matched,

        "bytes":
            p.stat().st_size,
    })


rows.sort(
    key=lambda x:
        (
            x["accession"],
            x["path"],
        )
)


with (
    OUT /
    "PHASE9F_A_FILE_INVENTORY.tsv"
).open(
    "w",
    newline=""
) as f:

    fields = [
        "accession",
        "path",
        "filename",
        "extension",
        "bytes",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(rows)


print(
    "INVENTORY_FILES=",
    len(rows)
)


for acc in sorted(
    set(
        x["accession"]
        for x in rows
    )
):

    subset = [
        x for x in rows
        if x["accession"] == acc
    ]

    print(
        acc,
        "files=",
        len(subset),
        "bytes=",
        sum(
            int(x["bytes"])
            for x in subset
        ),
    )

PY

}


# =============================================================================
# 08. DATASET ROLE / STATISTICAL VALIDITY MANIFEST
# =============================================================================

build_manifest() {

cat > "$OUT/CELLTYPE_DATASET_ANALYSIS_MANIFEST.tsv" <<'EOF'
accession	species	role	primary_gene_branch	inference_level	planned_analysis
GSE256490	Human	HUMAN_VASCULAR_REFERENCE	FGF5;SH3PXD2A;COL4A2;COL4A1	LOCALIZATION_ONLY	cell-type expression/detection;vascular lineage localization
GSE225948	Mouse	STROKE_IMMUNE_ENDOTHELIAL	SH3PXD2A;FGF5	DONOR_AWARE_DISEASE_VALIDATION	pseudobulk stroke-vs-sham by cell type/time
GSE189432	Mouse	STROKE_MYELOID	SH3PXD2A	DISEASE_STATE_VALIDATION	microglia/myeloid expression and stroke-state validation
GSE234052	Mouse	STROKE_PERICYTE	COL4A2;COL4A1	TEMPORAL_DESCRIPTIVE	pericyte expression across contra/ipisilateral and 1h/12h/24h; avoid overclaiming inferential DE
EOF

cat "$OUT/CELLTYPE_DATASET_ANALYSIS_MANIFEST.tsv"

}


# =============================================================================
# 09. READINESS
# =============================================================================

build_readiness() {

python3 - "$DATA" "$OUT" "$ENV" <<'PY'

from pathlib import Path
import csv
import sys


DATA = Path(sys.argv[1])
OUT = Path(sys.argv[2])
ENV = Path(sys.argv[3])


required = {
    "GSE189432": DATA / "GSE189432",
    "GSE234052": DATA / "GSE234052",
    "GSE225948": DATA / "GSE225948",
    "GSE256490": DATA / "GSE256490",
}


rows = []


all_present = True

for acc, path in required.items():

    files = [
        p for p in path.rglob("*")
        if p.is_file()
        and not p.name.endswith(".part")
    ]

    status = (
        "PRESENT"
        if files
        else "MISSING"
    )

    if not files:
        all_present = False

    rows.append(
        (
            acc,
            status
        )
    )


env_ok = (
    ENV / "bin/python"
).exists()


rows.extend([

    (
        "SINGLECELL_ENV",
        (
            "READY"
            if env_ok
            else "NOT_READY"
        )
    ),

    (
        "PHASE9F_A_ACQUISITION",
        (
            "COMPLETE"
            if (
                all_present
                and env_ok
            )
            else "PARTIAL"
        )
    ),

    (
        "HUMAN_REFERENCE_BRANCH",
        (
            "READY_FOR_PARSING"
            if (
                required[
                    "GSE256490"
                ].exists()
            )
            else "PENDING"
        )
    ),

    (
        "SH3PXD2A_DISEASE_BRANCH",
        (
            "READY_FOR_PARSING"
            if (
                required[
                    "GSE189432"
                ].exists()
                and required[
                    "GSE225948"
                ].exists()
            )
            else "PENDING"
        )
    ),

    (
        "COL4A2_PERICYTE_BRANCH",
        (
            "READY_FOR_PARSING"
            if required[
                "GSE234052"
            ].exists()
            else "PENDING"
        )
    ),

    (
        "NEXT_STAGE",
        "PHASE9F_B_FORMAT_AUDIT_AND_TARGET_EXPRESSION"
    ),
])


with (
    OUT /
    "PHASE9F_A_READINESS.tsv"
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
  "02_prepare_environment" \
  prepare_environment

run_step \
  "03_download_GSE189432" \
  download_gse189432

run_step \
  "04_download_GSE234052" \
  download_gse234052

run_step \
  "05_download_GSE225948" \
  download_gse225948

run_step \
  "06_download_GSE256490" \
  download_gse256490

run_step \
  "07_inventory_downloads" \
  inventory_downloads

run_step \
  "08_analysis_manifest" \
  build_manifest

run_step \
  "09_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9F-A COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== DATASET MANIFEST ====="

column -t -s $'\t' \
  "$OUT/CELLTYPE_DATASET_ANALYSIS_MANIFEST.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/CELLTYPE_DATASET_ANALYSIS_MANIFEST.tsv"


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9F_A_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE9F_A_READINESS.tsv"


echo
echo "===== DISK ====="

du -sh "$DATA"/* 2>/dev/null || true
df -h "$ROOT"


echo
echo "===== OUTPUTS ====="

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -printf '%f\t%s bytes\n' \
  | sort


echo
echo "OUTPUT_ROOT=$OUT"
echo "DATA_ROOT=$DATA"
echo "ENV=$ENV"
echo "LOG_ROOT=$LOGDIR"
