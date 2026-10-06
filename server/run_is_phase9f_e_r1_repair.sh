#!/usr/bin/env bash

ROOT="/srv/is-analysis"

P9FD="$ROOT/results/is/stage5_functional/phase9f_d_validated_pseudobulk"

PB="$P9FD/gse225948_pseudobulk"
DGE="$P9FD/gse225948_dge"

DATA="$ROOT/data/is/celltype"

OUT="$ROOT/results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff_r1"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9f_e_r1/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

RLIB="$ROOT/envs/is_Rlib"

RDS="$DATA/GSE256493/GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz"

REMOTE="gdrive:MASTER_DEGREE/IS_COLAB/phase9f_e_human_vascular/input"
REMOTE_FILE="$REMOTE/$(basename "$RDS")"

mkdir -p \
  "$OUT" \
  "$LOGDIR" \
  "$RLIB" \
  "$OUT/robust_dge"

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
  elif [[ $RC -eq 2 ]]; then
    printf "%s\tSKIP\t%s\n" "$name" "$RC" >> "$STATUS"
  else
    printf "%s\tFAIL\t%s\n" "$name" "$RC" >> "$STATUS"
  fi

  return 0
}


# =============================================================================
# 01. edgeR LIBRARY RECOVERY
# =============================================================================

recover_edger() {

  mkdir -p "$RLIB"

  echo "===== FILESYSTEM ====="

  find \
    "$RLIB" \
    "$HOME/R" \
    -maxdepth 5 \
    -type d \
    -name edgeR \
    -print \
    2>/dev/null \
    | sort -u


  echo
  echo "===== RAW R LIBRARY PATHS ====="

  Rscript - <<'RS'
cat("R=", R.version.string, "\n")

cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n")

cat(".libPaths BEFORE:\n")
print(.libPaths())

cat(
  "edgeR default namespace=",
  requireNamespace(
    "edgeR",
    quietly=TRUE
  ),
  "\n"
)
RS


  echo
  echo "===== EXPLICIT PROJECT LIBRARY ====="

  Rscript - \
    "$RLIB" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

RLIB <- args[[1]]

dir.create(
  RLIB,
  recursive=TRUE,
  showWarnings=FALSE
)

.libPaths(
  unique(
    c(
      RLIB,
      .libPaths()
    )
  )
)

cat(".libPaths AFTER:\n")
print(.libPaths())


project_ok <- requireNamespace(
  "edgeR",
  lib.loc=RLIB,
  quietly=TRUE
)


cat(
  "edgeR in project lib=",
  project_ok,
  "\n"
)


if (!project_ok) {

  cat(
    "edgeR not visible in project lib; installing/recovering...\n"
  )


  if (!requireNamespace(
        "BiocManager",
        quietly=TRUE
      )) {

    install.packages(
      "BiocManager",
      repos="https://cloud.r-project.org",
      lib=RLIB
    )

    .libPaths(
      unique(
        c(
          RLIB,
          .libPaths()
        )
      )
    )

  }


  BiocManager::install(
    "edgeR",
    lib=RLIB,
    ask=FALSE,
    update=FALSE
  )

}


if (!requireNamespace(
      "edgeR",
      lib.loc=RLIB,
      quietly=TRUE
    )) {

  stop(
    "edgeR still unavailable in explicit project library"
  )

}


library(
  edgeR,
  lib.loc=RLIB
)


cat(
  "edgeR_VERSION=",
  as.character(
    packageVersion(
      "edgeR"
    )
  ),
  "\n",
  sep=""
)


cat(
  "edgeR_PATH=",
  find.package(
    "edgeR",
    lib.loc=RLIB
  ),
  "\n",
  sep=""
)


cat(
  "EDGE_R_RECOVERY=PASS\n"
)

RS

}


# =============================================================================
# 02. TARGET TESTABILITY
# =============================================================================

build_target_testability() {

  Rscript - \
    "$RLIB" \
    "$PB" \
    "$OUT" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

RLIB <- args[[1]]
PB <- args[[2]]
OUT <- args[[3]]


.libPaths(
  unique(
    c(
      RLIB,
      .libPaths()
    )
  )
)


library(
  edgeR,
  lib.loc=RLIB
)


targets <- c(
  "Fgf5",
  "Aldh2",
  "Sh3pxd2a",
  "Col4a2",
  "Col4a1",
  "Calhm2",
  "Neurl1",
  "Ina"
)


manifest <- read.delim(
  file.path(
    PB,
    "PSEUDOBULK_MANIFEST.tsv"
  ),
  stringsAsFactors=FALSE,
  check.names=FALSE
)


strict <- tolower(
  as.character(
    manifest$strict_replicated_min3
  )
) %in% c(
  "true",
  "1",
  "t"
)


manifest <- manifest[
  strict,
  ,
  drop=FALSE
]


cat(
  "STRICT_CELLTYPES=",
  paste(
    manifest$celltype,
    collapse=","
  ),
  "\n",
  sep=""
)


rows <- list()


for (
  i in seq_len(
    nrow(manifest)
  )
) {

  celltype <- manifest$celltype[i]


  C <- read.delim(
    gzfile(
      manifest$counts_file[i]
    ),
    stringsAsFactors=FALSE,
    check.names=FALSE
  )


  D <- read.delim(
    manifest$design_file[i],
    stringsAsFactors=FALSE,
    check.names=FALSE
  )


  gene <- C$gene
  C$gene <- NULL


  counts <- as.matrix(C)

  storage.mode(
    counts
  ) <- "numeric"

  rownames(counts) <- gene


  if (
    anyDuplicated(
      rownames(counts)
    )
  ) {

    counts <- rowsum(
      counts,
      group=rownames(counts),
      reorder=FALSE
    )

  }


  counts <- counts[
    ,
    D$sample,
    drop=FALSE
  ]


  group <- factor(
    D$treatment,
    levels=c(
      "Sham",
      "D02",
      "D14"
    )
  )


  y <- DGEList(
    counts=counts,
    group=group
  )


  keep <- filterByExpr(
    y,
    group=group
  )


  kept <- rownames(y)[
    keep
  ]


  for (
    g in targets
  ) {

    present <- (
      g %in%
      rownames(y)
    )


    testable <- (
      g %in%
      kept
    )


    total <- (

      if (present)

        sum(
          y$counts[
            g,
            ,
            drop=TRUE
          ]
        )

      else

        NA_real_

    )


    status <- (

      if (!present)

        "ABSENT_FROM_MATRIX"

      else if (!testable)

        "NOT_TESTED_LOW_EXPRESSION"

      else

        "TESTABLE"

    )


    rows[[
      length(rows)+1
    ]] <- data.frame(

      celltype=
        celltype,

      gene=
        g,

      present_in_pseudobulk=
        present,

      total_count=
        total,

      passed_filterByExpr=
        testable,

      interpretation=
        status,

      stringsAsFactors=FALSE

    )

  }

}


X <- do.call(
  rbind,
  rows
)


write.table(
  X,
  file.path(
    OUT,
    "GSE225948_TARGET_TESTABILITY.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


print(X)

RS

}


# =============================================================================
# 03. ROBUST edgeR — STRICTLY REPLICATED CELL TYPES
# =============================================================================

run_robust_edger() {

  Rscript - \
    "$RLIB" \
    "$PB" \
    "$OUT/robust_dge" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

RLIB <- args[[1]]
PB <- args[[2]]
OUT <- args[[3]]


.libPaths(
  unique(
    c(
      RLIB,
      .libPaths()
    )
  )
)


library(
  edgeR,
  lib.loc=RLIB
)


dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)


targets <- c(
  "Fgf5",
  "Aldh2",
  "Sh3pxd2a",
  "Col4a2",
  "Col4a1",
  "Calhm2",
  "Neurl1",
  "Ina"
)


manifest <- read.delim(
  file.path(
    PB,
    "PSEUDOBULK_MANIFEST.tsv"
  ),
  stringsAsFactors=FALSE,
  check.names=FALSE
)


strict <- tolower(
  as.character(
    manifest$strict_replicated_min3
  )
) %in% c(
  "true",
  "1",
  "t"
)


manifest <- manifest[
  strict,
  ,
  drop=FALSE
]


target_results <- list()
run_rows <- list()


for (
  i in seq_len(
    nrow(manifest)
  )
) {

  celltype <- manifest$celltype[i]


  cat(
    "\n========================================\n"
  )

  cat(
    "CELLTYPE=",
    celltype,
    "\n",
    sep=""
  )


  C <- read.delim(
    gzfile(
      manifest$counts_file[i]
    ),
    stringsAsFactors=FALSE,
    check.names=FALSE
  )


  D <- read.delim(
    manifest$design_file[i],
    stringsAsFactors=FALSE,
    check.names=FALSE
  )


  gene <- C$gene
  C$gene <- NULL


  counts <- as.matrix(C)

  storage.mode(
    counts
  ) <- "numeric"

  rownames(counts) <- gene


  if (
    anyDuplicated(
      rownames(counts)
    )
  ) {

    counts <- rowsum(
      counts,
      group=rownames(counts),
      reorder=FALSE
    )

  }


  counts <- counts[
    ,
    D$sample,
    drop=FALSE
  ]


  group <- factor(
    D$treatment,
    levels=c(
      "Sham",
      "D02",
      "D14"
    )
  )


  y <- DGEList(
    counts=counts,
    group=group
  )


  keep <- filterByExpr(
    y,
    group=group
  )


  y <- y[
    keep,
    ,
    keep.lib.sizes=FALSE
  ]


  y <- calcNormFactors(
    y
  )


  design <- model.matrix(
    ~0 + group
  )


  colnames(
    design
  ) <- levels(
    group
  )


  y <- estimateDisp(
    y,
    design,
    robust=TRUE
  )


  fit <- glmQLFit(
    y,
    design,
    robust=TRUE
  )


  contrast_list <- list(

    D02_vs_Sham=
      makeContrasts(
        D02-Sham,
        levels=design
      ),

    D14_vs_Sham=
      makeContrasts(
        D14-Sham,
        levels=design
      )
  )


  for (
    contrast_name
    in names(
      contrast_list
    )
  ) {

    qlf <- glmQLFTest(
      fit,
      contrast=
        contrast_list[[
          contrast_name
        ]]
    )


    tab <- topTags(
      qlf,
      n=Inf,
      sort.by="PValue"
    )$table


    tab$gene <- rownames(
      tab
    )


    run_rows[[
      length(run_rows)+1
    ]] <- data.frame(

      celltype=
        celltype,

      contrast=
        contrast_name,

      samples=
        ncol(y),

      genes_tested=
        nrow(y),

      significant_fdr05=
        sum(
          tab$FDR < 0.05,
          na.rm=TRUE
        ),

      stringsAsFactors=FALSE
    )


    tt <- tab[
      tab$gene
      %in%
      targets,
      ,
      drop=FALSE
    ]


    if (
      nrow(tt) > 0
    ) {

      tt$gene <- rownames(
        tt
      )

      tt$celltype <- celltype

      tt$contrast <- contrast_name


      target_results[[
        length(target_results)+1
      ]] <- tt

    }

  }

}


RUN <- do.call(
  rbind,
  run_rows
)


write.table(
  RUN,
  file.path(
    OUT,
    "ROBUST_DGE_RUN_SUMMARY.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


if (
  length(
    target_results
  ) > 0
) {

  TAR <- do.call(
    rbind,
    target_results
  )


  TAR <- TAR[
    order(
      TAR$FDR,
      TAR$PValue
    ),
    ,
    drop=FALSE
  ]


  write.table(
    TAR,
    file.path(
      OUT,
      "ROBUST_TARGET_DGE.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )


  print(TAR)

} else {

  write.table(
    data.frame(
      status=
        "NO_PRIORITIZED_TARGET_PASSED_FILTER"
    ),
    file.path(
      OUT,
      "ROBUST_TARGET_DGE.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

}


cat(
  "\n===== RUN SUMMARY =====\n"
)

print(RUN)

RS

}


# =============================================================================
# 04. FREEZE MOUSE FUNCTIONAL LAYER
# =============================================================================

freeze_mouse_layer() {

python3 - \
  "$P9FD" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import sys

import pandas as pd


P9FD = Path(sys.argv[1])
OUT = Path(sys.argv[2])


g189 = pd.read_csv(
    P9FD /
    "GSE189432_TARGET_CLUSTER_SUMMARY_R1.tsv",
    sep="\t"
)


primary = pd.read_csv(
    P9FD /
    "GSE225948_TARGET_DGE_INTERPRETATION.tsv",
    sep="\t"
)


testability = pd.read_csv(
    OUT /
    "GSE225948_TARGET_TESTABILITY.tsv",
    sep="\t"
)


robust_path = (
    OUT /
    "robust_dge" /
    "ROBUST_TARGET_DGE.tsv"
)


try:

    robust = pd.read_csv(
        robust_path,
        sep="\t"
    )

except Exception:

    robust = pd.DataFrame()


sh = (
    g189[
        g189["gene"]
        == "Sh3pxd2a"
    ]
    .sort_values(
        [
            "detection_fraction",
            "cpm",
        ],
        ascending=False
    )
)


top_sh = (
    sh.iloc[0]
    if len(sh)
    else None
)


primary_sig = int(
    (
        primary["FDR"]
        < 0.05
    ).sum()
)


robust_sig = 0

if (
    not robust.empty
    and "FDR"
    in robust.columns
):

    robust_sig = int(
        (
            robust["FDR"]
            < 0.05
        ).sum()
    )


def status_for(gene):

    x = testability[
        testability["gene"]
        == gene
    ]


    if x.empty:

        return (
            "NO_STRICT_CELLTYPE_RECORD"
        )


    passed = (
        x[
            "passed_filterByExpr"
        ]
        .astype(str)
        .str.lower()
        .isin(
            [
                "true",
                "1",
                "t",
            ]
        )
    )


    if passed.any():

        return (
            "TESTED_IN_AT_LEAST_ONE_STRICT_CELLTYPE"
        )


    return (
        "NOT_TESTED_LOW_EXPRESSION_IN_STRICT_CELLTYPES"
    )


rows = [

    (
        "GSE189432_BARCODE_REPAIR",
        "PASS_8_OF_8_100_PERCENT_MATCH"
    ),

    (
        "SH3PXD2A_GSE189432_TOP_CLUSTER",
        (
            str(
                top_sh[
                    "cluster"
                ]
            )
            if top_sh
            is not None
            else "NA"
        )
    ),

    (
        "SH3PXD2A_GSE189432_TOP_DETECTION",
        (
            str(
                top_sh[
                    "detection_fraction"
                ]
            )
            if top_sh
            is not None
            else "NA"
        )
    ),

    (
        "SH3PXD2A_GSE189432_INTERPRETATION",
        "MYELOID_MICROGLIAL_LOCALIZATION_DESCRIPTIVE"
    ),

    (
        "GSE225948_PRIMARY_STRATUM",
        "BRAIN_W8_SHAM3_D02_4_D14_3"
    ),

    (
        "GSE225948_STRICT_PRIMARY_CELLTYPES",
        "EC1;MdC3;Mg1;Tc1"
    ),

    (
        "PRIORITIZED_TARGETS_FDR05_PRIMARY",
        str(
            primary_sig
        )
    ),

    (
        "PRIORITIZED_TARGETS_FDR05_ROBUST_SENSITIVITY",
        str(
            robust_sig
        )
    ),

    (
        "FGF5_DGE_TESTABILITY",
        status_for(
            "Fgf5"
        )
    ),

    (
        "ALDH2_DGE_TESTABILITY",
        status_for(
            "Aldh2"
        )
    ),

    (
        "SH3PXD2A_DGE_TESTABILITY",
        status_for(
            "Sh3pxd2a"
        )
    ),

    (
        "COL4A2_DGE_TESTABILITY",
        status_for(
            "Col4a2"
        )
    ),

    (
        "COL4A1_DGE_TESTABILITY",
        status_for(
            "Col4a1"
        )
    ),

    (
        "COL4A1_COL4A2_EC1_D02",
        "NOMINAL_UPREGULATION_NOT_FDR_SIGNIFICANT"
    ),

    (
        "MOUSE_LAYER_CONCLUSION",
        "CELLTYPE_LOCALIZATION_SUPPORTED_TARGET_SPECIFIC_STROKE_DGE_NOT_ESTABLISHED"
    ),
]


with (
    OUT /
    "MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv"
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
        "status",
    ])

    w.writerows(
        rows
    )


for a, b in rows:

    print(
        f"{a} = {b}"
    )

PY

}


# =============================================================================
# 05. VERIFY / RECOVER HUMAN RDS DRIVE COPY
# =============================================================================

verify_human_upload() {

  if [[ ! -s "$RDS" ]]; then

    echo "LOCAL_RDS_MISSING"

    return 1
  fi


  if ! command -v rclone >/dev/null 2>&1; then

    echo "RCLONE_MISSING"

    return 1
  fi


  LOCAL_SIZE="$(
    stat -c%s "$RDS"
  )"


  echo "LOCAL_SIZE=$LOCAL_SIZE"


  REMOTE_SIZE="$(
    rclone size \
      "$REMOTE_FILE" \
      --json \
      2>/dev/null \
      | python3 -c '
import json,sys
try:
    x=json.load(sys.stdin)
    print(x.get("bytes",0))
except Exception:
    print(0)
'
  )"


  echo "REMOTE_SIZE_INITIAL=$REMOTE_SIZE"


  if [[ "$REMOTE_SIZE" != "$LOCAL_SIZE" ]]; then

    echo
    echo "REMOTE_COPY_INCOMPLETE_OR_MISSING"
    echo "RUNNING_COPYTO"

    rclone copyto \
      "$RDS" \
      "$REMOTE_FILE" \
      --progress

    RC=$?

    if [[ $RC -ne 0 ]]; then

      echo "RCLONE_COPY_FAIL=$RC"

      return 1
    fi

  else

    echo "REMOTE_ALREADY_COMPLETE"

  fi


  REMOTE_SIZE="$(
    rclone size \
      "$REMOTE_FILE" \
      --json \
      | python3 -c '
import json,sys
x=json.load(sys.stdin)
print(x.get("bytes",0))
'
  )"


  echo "REMOTE_SIZE_FINAL=$REMOTE_SIZE"


  if [[ "$REMOTE_SIZE" != "$LOCAL_SIZE" ]]; then

    echo "REMOTE_SIZE_MISMATCH"

    return 1
  fi


  sha256sum \
    "$RDS" \
    > "$OUT/GSE256493_RDS_SHA256.txt"


  cat > "$OUT/GSE256493_DRIVE_HANDOFF.tsv" <<EOF
component	status
LOCAL_RDS	PRESENT
LOCAL_SIZE	$LOCAL_SIZE
REMOTE_RDS	PRESENT_SIZE_MATCH
REMOTE_SIZE	$REMOTE_SIZE
REMOTE_PATH	$REMOTE_FILE
COLAB_STATUS	READY
EOF


  cat \
    "$OUT/GSE256493_DRIVE_HANDOFF.tsv"


  echo "HUMAN_UPLOAD=PASS"
}


# =============================================================================
# 06. READINESS
# =============================================================================

build_readiness() {

python3 - \
  "$STATUS" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


STATUS = Path(sys.argv[1])
OUT = Path(sys.argv[2])


with STATUS.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


required = {

    "01_recover_edgeR",
    "02_target_testability",
    "03_robust_edgeR_sensitivity",
    "04_freeze_mouse_layer",
    "05_verify_human_upload",
}


state = {
    x["step"]:
        x["status"]
    for x in rows
}


ok = all(
    state.get(x)
    == "PASS"
    for x in required
)


readiness = [

    (
        "PHASE9F_E_R1",
        (
            "COMPLETE"
            if ok
            else "PARTIAL"
        )
    ),

    (
        "EDGER_LIBRARY",
        (
            "RECOVERED_EXPLICIT_PROJECT_LIB"
            if state.get(
                "01_recover_edgeR"
            ) == "PASS"
            else "FAILED"
        )
    ),

    (
        "MOUSE_FUNCTIONAL_LAYER",
        (
            "FROZEN"
            if state.get(
                "04_freeze_mouse_layer"
            ) == "PASS"
            else "PENDING"
        )
    ),

    (
        "GSE234052",
        "DEFER_RAW_WHITELIST_CELL_CALLING"
    ),

    (
        "GSE256493",
        (
            "DRIVE_HANDOFF_COMPLETE"
            if state.get(
                "05_verify_human_upload"
            ) == "PASS"
            else "UPLOAD_PENDING"
        )
    ),

    (
        "NEXT_STAGE",
        (
            "PHASE9F_E_COLAB_HUMAN_AUTHOR_ANNOTATION"
            if ok
            else "REPAIR_FAILED_STEP_FIRST"
        )
    ),
]


with (
    OUT /
    "PHASE9F_E_R1_READINESS.tsv"
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
        "status",
    ])

    w.writerows(
        readiness
    )


for a, b in readiness:

    print(
        f"{a} = {b}"
    )

PY

}


# =============================================================================
# RUN
# =============================================================================

run_step \
  "01_recover_edgeR" \
  recover_edger

run_step \
  "02_target_testability" \
  build_target_testability

run_step \
  "03_robust_edgeR_sensitivity" \
  run_robust_edger

run_step \
  "04_freeze_mouse_layer" \
  freeze_mouse_layer

run_step \
  "05_verify_human_upload" \
  verify_human_upload

run_step \
  "06_build_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9F-E R1 COMPLETE"
echo "================================================================================"


echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== EDGE R ====="

Rscript - \
  "$RLIB" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

RLIB <- args[[1]]

.libPaths(
  c(
    RLIB,
    .libPaths()
  )
)

cat(
  "edgeR=",
  as.character(
    packageVersion(
      "edgeR"
    )
  ),
  "\n",
  sep=""
)

cat(
  "path=",
  find.package(
    "edgeR"
  ),
  "\n",
  sep=""
)

RS


echo
echo "===== TARGET TESTABILITY ====="

column -t -s $'\t' \
  "$OUT/GSE225948_TARGET_TESTABILITY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== ROBUST TARGET DGE ====="

column -t -s $'\t' \
  "$OUT/robust_dge/ROBUST_TARGET_DGE.tsv" \
  2>/dev/null \
  || true


echo
echo "===== MOUSE FREEZE ====="

column -t -s $'\t' \
  "$OUT/MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv" \
  2>/dev/null \
  || true


echo
echo "===== HUMAN HANDOFF ====="

cat \
  "$OUT/GSE256493_DRIVE_HANDOFF.tsv" \
  2>/dev/null \
  || true


echo
echo "===== READINESS ====="

cat \
  "$OUT/PHASE9F_E_R1_READINESS.tsv" \
  2>/dev/null \
  || true


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
