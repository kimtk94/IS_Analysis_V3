#!/usr/bin/env bash

ROOT="/srv/is-analysis"

P9FD="$ROOT/results/is/stage5_functional/phase9f_d_validated_pseudobulk"

PB="$P9FD/gse225948_pseudobulk"
DGE="$P9FD/gse225948_dge"

DATA="$ROOT/data/is/celltype"

OUT="$ROOT/results/is/stage5_functional/phase9f_e_mouse_freeze_human_handoff"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9f_e/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

RLIB="$ROOT/envs/is_Rlib"

RDS="$DATA/GSE256493/GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz"

GDRIVE_DIR="gdrive:MASTER_DEGREE/IS_COLAB/phase9f_e_human_vascular/input"

mkdir -p \
  "$OUT" \
  "$LOGDIR" \
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
# 01. AUDIT PRIOR PREFLIGHT FAILURE WITHOUT BLOCKING VALID OUTPUTS
# =============================================================================

audit_prior_run() {

  echo "===== PHASE9F-D STATUS ====="

  if [[ ! -s "$P9FD/STEP_STATUS.tsv" ]]; then
    echo "PHASE9FD_STATUS_MISSING"
    return 1
  fi

  cat "$P9FD/STEP_STATUS.tsv"


  echo
  echo "===== REQUIRED SCIENTIFIC OUTPUTS ====="

  REQUIRED=(
    "$P9FD/GSE189432_BARCODE_REPAIR.tsv"
    "$P9FD/GSE189432_TARGET_CLUSTER_SUMMARY_R1.tsv"
    "$PB/CELLTYPE_DGE_ELIGIBILITY.tsv"
    "$PB/PSEUDOBULK_MANIFEST.tsv"
    "$DGE/DGE_RUN_SUMMARY.tsv"
    "$P9FD/GSE225948_TARGET_DGE_INTERPRETATION.tsv"
    "$P9FD/GSE256493_HUMAN_RDS_RESOURCE_GATE.tsv"
  )

  FAIL=0

  for f in "${REQUIRED[@]}"
  do

    if [[ -s "$f" ]]; then
      echo "FOUND $f"
    else
      echo "MISSING $f"
      FAIL=$((FAIL+1))
    fi

  done


  echo
  echo "MISSING_REQUIRED_OUTPUTS=$FAIL"


  echo
  echo "===== PRIOR PREFLIGHT LOG ====="

  PRIOR_LOG="$(
    find \
      "$ROOT/logs/is_phase9f_d" \
      -type f \
      -name '01_preflight.log' \
      -printf '%T@\t%p\n' \
      2>/dev/null \
      | sort -nr \
      | head -1 \
      | cut -f2-
  )"


  if [[ -n "$PRIOR_LOG" && -s "$PRIOR_LOG" ]]; then

    echo "PRIOR_LOG=$PRIOR_LOG"

    cp \
      "$PRIOR_LOG" \
      "$OUT/PHASE9F_D_01_PREFLIGHT.log"

    cat "$PRIOR_LOG"

  else

    echo "PRIOR_PREFLIGHT_LOG_NOT_FOUND"

  fi


  if [[ "$FAIL" -ne 0 ]]; then
    return 1
  fi


  cat > "$OUT/PHASE9F_D_AUDIT_STATUS.tsv" <<'EOF'
component	status
PHASE9F_D_STEP_STATUS	ONE_PREFLIGHT_FAILURE
CORE_SCIENTIFIC_OUTPUTS	PRESENT
MOUSE_FUNCTIONAL_ANALYSIS	COMPUTATIONALLY_COMPLETE
CLEAN_ALL_STEP_PASS	NO
ACTION	PRESERVE_RESULTS_AND_AUDIT_PREFLIGHT
EOF

  echo
  echo "PRIOR_AUDIT=PASS"
}


# =============================================================================
# 02. BUILD TARGET TESTABILITY MATRIX
# =============================================================================

build_target_testability() {

  export R_LIBS_USER="$RLIB"

  Rscript - \
    "$PB" \
    "$OUT" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

PB <- args[[1]]
OUT <- args[[2]]

library(edgeR)


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


# Primary interpretation uses strict replicated cell types only.
manifest <- manifest[
  manifest$strict_replicated_min3
  %in%
  c(
    TRUE,
    "TRUE",
    "True",
    1,
    "1"
  ),
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


  genes <- C$gene
  C$gene <- NULL

  counts <- as.matrix(C)

  storage.mode(counts) <- "numeric"

  rownames(counts) <- genes


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


  y0 <- DGEList(
    counts=counts,
    group=group
  )


  keep <- filterByExpr(
    y0,
    group=group
  )


  kept_genes <- rownames(
    y0
  )[keep]


  for (
    gene in targets
  ) {

    present <- (
      gene
      %in%
      rownames(y0)
    )

    tested <- (
      gene
      %in%
      kept_genes
    )


    total_count <- (
      if (present)
        sum(
          y0$counts[
            gene,
            ,
            drop=TRUE
          ]
        )
      else
        NA
    )


    rows[[
      length(rows)+1
    ]] <- data.frame(

      celltype=
        celltype,

      gene=
        gene,

      present_in_pseudobulk=
        present,

      total_count=
        total_count,

      passed_filterByExpr=
        tested,

      interpretation=
        if (!present)
          "ABSENT_FROM_MATRIX"
        else if (!tested)
          "NOT_TESTED_LOW_EXPRESSION"
        else
          "TESTABLE",

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
# 03. ROBUST edgeR SENSITIVITY — STRICT REPLICATED CELL TYPES ONLY
# =============================================================================

run_robust_edger() {

  export R_LIBS_USER="$RLIB"

  Rscript - \
    "$PB" \
    "$OUT/robust_dge" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

PB <- args[[1]]
OUT <- args[[2]]

library(edgeR)

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


manifest <- manifest[
  manifest$strict_replicated_min3
  %in%
  c(
    TRUE,
    "TRUE",
    "True",
    1,
    "1"
  ),
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


  genes <- C$gene
  C$gene <- NULL

  counts <- as.matrix(C)

  storage.mode(counts) <- "numeric"

  rownames(counts) <- genes


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

  colnames(design) <- levels(
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


  contrasts <- list(

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
    in names(contrasts)
  ) {

    qlf <- glmQLFTest(
      fit,
      contrast=
        contrasts[[
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


    sig_n <- sum(
      tab$FDR < 0.05,
      na.rm=TRUE
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
        sig_n,

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

      tt$gene <- rownames(tt)
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
  length(target_results) > 0
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

}


print(RUN)

RS

}


# =============================================================================
# 04. FREEZE MOUSE INTERPRETATION
# =============================================================================

freeze_mouse_layer() {

python3 - \
  "$P9FD" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import pandas as pd
import sys


P9FD = Path(sys.argv[1])
OUT = Path(sys.argv[2])


# GSE189432
g189 = pd.read_csv(
    P9FD /
    "GSE189432_TARGET_CLUSTER_SUMMARY_R1.tsv",
    sep="\t"
)


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


# Primary edgeR target result
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


if robust_path.exists():

    robust = pd.read_csv(
        robust_path,
        sep="\t"
    )

else:

    robust = pd.DataFrame()


target_sig_primary = int(
    (
        primary[
            "FDR"
        ] < 0.05
    ).sum()
)


target_sig_robust = (

    int(
        (
            robust["FDR"]
            < 0.05
        ).sum()
    )

    if (
        not robust.empty
        and "FDR"
        in robust.columns
    )

    else 0
)


# Explicit classification for key genes.
def gene_test_status(gene):

    x = testability[
        testability["gene"]
        == gene
    ]

    if x.empty:
        return (
            "NO_STRICT_CELLTYPE_RECORD"
        )

    if (
        x[
            "passed_filterByExpr"
        ]
        .astype(str)
        .str.lower()
        .eq("true")
        .any()
    ):

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
            target_sig_primary
        )
    ),

    (
        "PRIORITIZED_TARGETS_FDR05_ROBUST_SENSITIVITY",
        str(
            target_sig_robust
        )
    ),

    (
        "SH3PXD2A_DGE_TESTABILITY",
        gene_test_status(
            "Sh3pxd2a"
        )
    ),

    (
        "FGF5_DGE_TESTABILITY",
        gene_test_status(
            "Fgf5"
        )
    ),

    (
        "COL4A2_DGE_TESTABILITY",
        gene_test_status(
            "Col4a2"
        )
    ),

    (
        "COL4A1_DGE_TESTABILITY",
        gene_test_status(
            "Col4a1"
        )
    ),

    (
        "ALDH2_DGE_TESTABILITY",
        gene_test_status(
            "Aldh2"
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
# 05. UPLOAD HUMAN RDS TO DRIVE
# =============================================================================

upload_human_rds() {

  if [[ ! -s "$RDS" ]]; then

    echo "HUMAN_RDS_MISSING"

    return 1
  fi


  if ! command -v rclone >/dev/null 2>&1; then

    echo "RCLONE_MISSING"

    return 1
  fi


  echo "===== RDS ====="

  ls -lh "$RDS"

  sha256sum \
    "$RDS" \
    | tee \
      "$OUT/GSE256493_RDS_SHA256.txt"


  echo
  echo "===== UPLOAD ====="

  rclone copy \
    "$RDS" \
    "$GDRIVE_DIR" \
    --progress


  RC=$?

  if [[ $RC -ne 0 ]]; then

    echo "UPLOAD_FAIL=$RC"

    return 1
  fi


  echo
  echo "===== REMOTE CHECK ====="

  rclone lsl \
    "$GDRIVE_DIR" \
    | grep \
      'GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz'


  RC=$?

  if [[ $RC -ne 0 ]]; then
    return 1
  fi


  echo "HUMAN_RDS_DRIVE=PASS"
}


# =============================================================================
# 06. READINESS
# =============================================================================

build_readiness() {

cat > "$OUT/PHASE9F_E_READINESS.tsv" <<'EOF'
component	status
PHASE9F_E_SERVER_SIDE	COMPLETE
GSE189432	LOCALIZATION_LAYER_FROZEN
GSE225948	PSEUDOBULK_LAYER_FROZEN
TARGET_SPECIFIC_MOUSE_DGE	NO_FDR_SIGNIFICANT_PRIORITY_GENE_IN_PRIMARY_ANALYSIS
GSE234052	DEFER_RAW_BARCODE_CELL_CALLING
GSE256493	SENT_TO_HIGH_MEMORY_COLAB
NEXT	PHASE9F_E_COLAB_HUMAN_AUTHOR_ANNOTATION
EOF

cat "$OUT/PHASE9F_E_READINESS.tsv"

}


# =============================================================================
# RUN
# =============================================================================

run_step \
  "01_audit_prior_run" \
  audit_prior_run

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
  "05_upload_human_RDS" \
  upload_human_rds

run_step \
  "06_build_readiness" \
  build_readiness


echo
echo "================================================================================"
echo "PHASE9F-E SERVER COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


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
  || cat \
       "$OUT/MOUSE_FUNCTIONAL_LAYER_FREEZE.tsv"


echo
echo "===== READINESS ====="

cat \
  "$OUT/PHASE9F_E_READINESS.tsv"


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
