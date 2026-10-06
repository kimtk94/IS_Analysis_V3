#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

P9FC="$ROOT/results/is/stage5_functional/phase9f_c_target_expression"

DATA="$ROOT/data/is/celltype"
OUT="$ROOT/results/is/stage5_functional/phase9f_d_validated_pseudobulk"

ENV="$ROOT/envs/is_singlecell"
RLIB="$ROOT/envs/is_Rlib"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9f_d/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

mkdir -p \
  "$OUT" \
  "$LOGDIR" \
  "$RLIB" \
  "$OUT/gse225948_pseudobulk" \
  "$OUT/gse225948_dge"

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
# 01. PREFLIGHT
# =============================================================================

preflight() {

  echo "===== PHASE9F-C STATUS ====="

  if [[ ! -s "$P9FC/STEP_STATUS.tsv" ]]; then
    echo "PHASE9FC_STATUS_MISSING"
    return 1
  fi

  cat "$P9FC/STEP_STATUS.tsv"

  echo
  echo "===== EXPECTED KNOWN FAILURE ====="

  BAD="$(
    awk -F'\t' '
      NR>1 &&
      $2=="FAIL" &&
      $1!="03_GSE189432_target_expression"
      {n++}
      END {print n+0}
    ' "$P9FC/STEP_STATUS.tsv"
  )"

  echo "UNEXPECTED_FAILURES=$BAD"

  if [[ "$BAD" -ne 0 ]]; then
    return 1
  fi


  echo
  echo "===== INPUTS ====="

  for f in \
    "$DATA/GSE189432/GSE189432_annotations.csv.gz" \
    "$P9FC/GSE225948_TARGET_PSEUDOBULK.tsv"
  do

    if [[ -s "$f" ]]; then
      echo "FOUND $f"
    else
      echo "MISSING $f"
      return 1
    fi
  done


  if [[ ! -x "$ENV/bin/python" ]]; then
    echo "PYTHON_ENV_MISSING"
    return 1
  fi


  "$ENV/bin/python" - <<'PY'
import numpy
import pandas
import scipy
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
PY


  echo
  echo "===== MEMORY ====="

  free -h

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. CORRECT PHASE9F-C BOOKKEEPING
# =============================================================================

correct_phase9fc_status() {

cat > "$OUT/PHASE9F_C_CORRECTED_STATUS.tsv" <<'EOF'
component	status
PHASE9F_C	PARTIAL
GSE189432_TARGET_EXPRESSION	FAILED_BARCODE_MATCH_PENDING_REPAIR
GSE225948_TARGET_PSEUDOBULK	COMPLETE
GSE234052_MATRIX	CELL_CALLING_REQUIRED
GSE256493_RDS	ACQUIRED_GZIP_VALIDATED
HUMAN_RDS_PARSE	PENDING_RESOURCE_GATE
EOF

cat "$OUT/PHASE9F_C_CORRECTED_STATUS.tsv"

}


# =============================================================================
# 03. GSE189432 BARCODE REPAIR + TARGET EXPRESSION
# =============================================================================

repair_gse189432() {

"$ENV/bin/python" - \
  "$DATA/GSE189432" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import re
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

ANNOT = ROOT / "GSE189432_annotations.csv.gz"

TARGETS = [
    "Fgf5",
    "Aldh2",
    "Sh3pxd2a",
    "Col4a2",
    "Col4a1",
    "Calhm2",
    "Neurl1",
    "Ina",
]


ann = pd.read_csv(
    ANNOT,
    compression="gzip"
)

ann["barcode"] = (
    ann["barcode"]
    .astype(str)
)

ann["sample"] = (
    ann["sample"]
    .astype(str)
)


# ------------------------------------------------------------------
# Extract canonical 10x 16-nt barcode core.
# No fuzzy matching.
# ------------------------------------------------------------------

DNA16 = re.compile(
    r"([ACGT]{16})",
    flags=re.I
)


def barcode_core(x):

    m = DNA16.search(
        str(x)
    )

    if not m:
        return ""

    return (
        m.group(1)
        .upper()
    )


def read_lines(path):

    with gzip.open(
        path,
        "rt",
        errors="replace"
    ) as f:

        return [
            x.strip()
            for x in f
            if x.strip()
        ]


matrices = sorted(
    ROOT.rglob(
        "*_matrix.mtx.gz"
    )
)


diag_rows = []
expr_rows = []


for matrix_path in matrices:

    prefix = matrix_path.name[
        :-len("_matrix.mtx.gz")
    ]

    sample = re.sub(
        r"^GSM\d+_",
        "",
        prefix
    )

    barcode_path = (
        matrix_path.parent /
        f"{prefix}_barcodes.tsv.gz"
    )

    feature_path = (
        matrix_path.parent /
        f"{prefix}_features.tsv.gz"
    )


    if not (
        barcode_path.exists()
        and feature_path.exists()
    ):
        raise RuntimeError(
            f"Incomplete triplet: {prefix}"
        )


    raw_barcodes = read_lines(
        barcode_path
    )


    meta = ann[
        ann["sample"] == sample
    ].copy()


    raw_core = [
        barcode_core(x)
        for x in raw_barcodes
    ]

    ann_core = [
        barcode_core(x)
        for x in meta["barcode"]
    ]


    raw_counts = {}

    for i, c in enumerate(
        raw_core
    ):

        if not c:
            continue

        raw_counts.setdefault(
            c,
            []
        ).append(i)


    unique_raw = {
        c: idxs[0]
        for c, idxs
        in raw_counts.items()
        if len(idxs) == 1
    }


    ann_counts = {}

    for c in ann_core:

        if not c:
            continue

        ann_counts[c] = (
            ann_counts.get(c, 0)
            + 1
        )


    matched_indices = []

    ambiguous = 0
    unresolved = 0


    for c in ann_core:

        if not c:

            matched_indices.append(
                -1
            )

            unresolved += 1
            continue


        if ann_counts.get(
            c,
            0
        ) != 1:

            matched_indices.append(
                -1
            )

            ambiguous += 1
            continue


        idx = unique_raw.get(
            c,
            -1
        )

        matched_indices.append(
            idx
        )

        if idx < 0:
            unresolved += 1


    matched_n = sum(
        x >= 0
        for x in matched_indices
    )


    rate = (
        matched_n /
        len(meta)
        if len(meta)
        else 0
    )


    diag_rows.append({

        "sample":
            sample,

        "raw_cells":
            len(raw_barcodes),

        "annotated_cells":
            len(meta),

        "matched_cells":
            matched_n,

        "match_rate":
            rate,

        "ambiguous_annotation_barcodes":
            ambiguous,

        "unresolved_annotation_barcodes":
            unresolved,

        "raw_example":
            ";".join(
                raw_barcodes[:5]
            ),

        "annotation_example":
            ";".join(
                meta[
                    "barcode"
                ]
                .head(5)
                .astype(str)
            ),

        "strategy":
            "UNIQUE_16BP_10X_CORE",

        "status":
            (
                "PASS"
                if rate >= 0.95
                else "FAIL"
            ),
    })


    print(
        sample,
        "raw=",
        len(raw_barcodes),
        "annotated=",
        len(meta),
        "matched=",
        matched_n,
        "rate=",
        rate,
    )


    if rate < 0.95:

        continue


    meta["matrix_col"] = (
        matched_indices
    )

    meta = meta[
        meta["matrix_col"]
        >= 0
    ].copy()


    # ------------------------------------------------------------------
    # Feature names
    # ------------------------------------------------------------------

    genes = []

    with gzip.open(
        feature_path,
        "rt",
        errors="replace"
    ) as f:

        for line in f:

            x = line.rstrip(
                "\n"
            ).split("\t")

            if len(x) >= 2:
                genes.append(
                    x[1]
                )
            elif x:
                genes.append(
                    x[0]
                )
            else:
                genes.append("")


    gene_index = {}

    for i, g in enumerate(
        genes
    ):

        if (
            g in TARGETS
            and g not in gene_index
        ):
            gene_index[g] = i


    # ------------------------------------------------------------------
    # Matrix
    # ------------------------------------------------------------------

    with gzip.open(
        matrix_path,
        "rb"
    ) as f:

        X = mmread(
            f
        ).tocsr()


    if X.shape[1] != len(
        raw_barcodes
    ):

        raise RuntimeError(
            f"{sample}: matrix-barcode mismatch"
        )


    total_umi = np.asarray(
        X.sum(
            axis=0
        )
    ).ravel()


    clusters = (
        meta["cluster"]
        .astype(str)
        .to_numpy()
    )

    cols = (
        meta["matrix_col"]
        .astype(int)
        .to_numpy()
    )


    for cluster in sorted(
        set(clusters)
    ):

        select = (
            clusters
            == cluster
        )

        cc = cols[
            select
        ]

        n_cells = len(cc)

        library_size = float(
            total_umi[
                cc
            ].sum()
        )


        for gene in TARGETS:

            if gene not in gene_index:

                expr_rows.append({

                    "sample":
                        sample,

                    "cluster":
                        cluster,

                    "gene":
                        gene,

                    "n_cells":
                        n_cells,

                    "gene_present":
                        False,

                    "sum_counts":
                        "",

                    "detected_cells":
                        "",

                    "detection_fraction":
                        "",

                    "library_size":
                        library_size,

                    "cpm":
                        "",
                })

                continue


            gi = gene_index[
                gene
            ]


            vals = (
                X[
                    gi,
                    cc
                ]
                .toarray()
                .ravel()
            )


            s = float(
                vals.sum()
            )

            detected = int(
                np.sum(
                    vals > 0
                )
            )


            expr_rows.append({

                "sample":
                    sample,

                "cluster":
                    cluster,

                "gene":
                    gene,

                "n_cells":
                    n_cells,

                "gene_present":
                    True,

                "sum_counts":
                    s,

                "detected_cells":
                    detected,

                "detection_fraction":
                    (
                        detected /
                        n_cells
                        if n_cells
                        else 0
                    ),

                "library_size":
                    library_size,

                "cpm":
                    (
                        s /
                        library_size *
                        1e6
                        if library_size > 0
                        else 0
                    ),
            })


diag = pd.DataFrame(
    diag_rows
)

diag.to_csv(
    OUT /
    "GSE189432_BARCODE_REPAIR.tsv",
    sep="\t",
    index=False
)


if not (
    diag["status"]
    == "PASS"
).all():

    print()
    print(
        "BARCODE_REPAIR_INCOMPLETE"
    )

    print(
        diag.to_string(
            index=False
        )
    )

    raise SystemExit(1)


expr = pd.DataFrame(
    expr_rows
)


expr.to_csv(
    OUT /
    "GSE189432_TARGET_CLUSTER_EXPRESSION_R1.tsv",
    sep="\t",
    index=False
)


usable = expr[
    expr[
        "gene_present"
    ] == True
].copy()


summary = (
    usable.groupby(
        [
            "cluster",
            "gene",
        ]
    )
    .agg(
        samples=(
            "sample",
            "nunique"
        ),
        cells=(
            "n_cells",
            "sum"
        ),
        sum_counts=(
            "sum_counts",
            "sum"
        ),
        detected_cells=(
            "detected_cells",
            "sum"
        ),
        library_size=(
            "library_size",
            "sum"
        ),
    )
    .reset_index()
)


summary[
    "detection_fraction"
] = (
    summary[
        "detected_cells"
    ]
    /
    summary[
        "cells"
    ]
)


summary[
    "cpm"
] = (
    summary[
        "sum_counts"
    ]
    /
    summary[
        "library_size"
    ]
    * 1e6
)


summary.to_csv(
    OUT /
    "GSE189432_TARGET_CLUSTER_SUMMARY_R1.tsv",
    sep="\t",
    index=False
)


print()
print(
    "===== SH3PXD2A CLUSTERS ====="
)

print(
    summary[
        summary["gene"]
        == "Sh3pxd2a"
    ]
    .sort_values(
        [
            "detection_fraction",
            "cpm",
        ],
        ascending=False
    )
    .head(20)
    .to_string(
        index=False
    )
)

PY

}


# =============================================================================
# 04. GSE225948 — FULL GENE PSEUDOBULK FOR BRAIN / W8
# =============================================================================

build_gse225948_full_pseudobulk() {

"$ENV/bin/python" - \
  "$DATA/GSE225948/processed" \
  "$OUT/gse225948_pseudobulk" <<'PY'

from pathlib import Path
import csv
import gzip
import re
import sys

import numpy as np
import pandas as pd


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


metadata_files = sorted(
    ROOT.glob(
        "*_metadata.csv.gz"
    )
)


primary = []


for meta_path in metadata_files:

    sample = meta_path.name.replace(
        "_metadata.csv.gz",
        ""
    )


    m = pd.read_csv(
        meta_path,
        compression="gzip",
        index_col=0
    )


    def one(col):

        if col not in m.columns:
            return ""

        x = (
            m[col]
            .dropna()
            .astype(str)
            .unique()
        )

        if len(x) != 1:
            return ";".join(
                sorted(x)
            )

        return str(
            x[0]
        )


    tissue = one(
        "tissue"
    )

    age = one(
        "age"
    )

    treatment = one(
        "treatment"
    )


    if (
        tissue.lower()
        == "brain"
        and age
        == "W8"
        and treatment
        in {
            "Sham",
            "D02",
            "D14",
        }
    ):

        primary.append({

            "sample":
                sample,

            "metadata":
                meta_path,

            "counts":
                ROOT /
                f"{sample}_counts.csv.gz",

            "treatment":
                treatment,
        })


print(
    "PRIMARY_SAMPLES=",
    len(primary)
)


for x in primary:

    print(
        x["sample"],
        x["treatment"]
    )


counts_by_column = {}
cell_design = []
gene_order = None


for info in primary:

    sample = info[
        "sample"
    ]

    treatment = info[
        "treatment"
    ]

    meta = pd.read_csv(
        info["metadata"],
        compression="gzip",
        index_col=0
    )

    meta.index = (
        meta.index
        .astype(str)
    )


    if "sub.celltype" not in meta.columns:

        raise RuntimeError(
            f"{sample}: sub.celltype missing"
        )


    counts_path = info[
        "counts"
    ]


    if not counts_path.exists():

        raise RuntimeError(
            f"{sample}: counts missing"
        )


    with gzip.open(
        counts_path,
        "rt",
        errors="replace"
    ) as f:

        reader = csv.reader(f)

        header = next(
            reader
        )

        cells = [
            str(x)
            for x in header[
                1:
            ]
        ]


        if len(cells) != len(
            meta
        ):

            raise RuntimeError(
                f"{sample}: cell count mismatch"
            )


        if not all(
            x in meta.index
            for x in cells
        ):

            raise RuntimeError(
                f"{sample}: cell IDs do not match metadata"
            )


        meta = meta.reindex(
            cells
        )


        groups = (
            meta[
                "sub.celltype"
            ]
            .fillna(
                "NA"
            )
            .astype(str)
        )


        categories = sorted(
            x
            for x in groups.unique()
            if x != "NA"
        )


        cat_to_code = {
            c: i
            for i, c
            in enumerate(
                categories
            )
        }


        codes = np.array(
            [
                cat_to_code.get(
                    x,
                    -1
                )
                for x in groups
            ],
            dtype=int
        )


        valid = (
            codes >= 0
        )


        ncells = np.bincount(
            codes[
                valid
            ],
            minlength=len(
                categories
            )
        )


        sample_vectors = {
            c: []
            for c in categories
        }


        local_genes = []


        for row_idx, row in enumerate(
            reader
        ):

            if not row:
                continue


            gene = str(
                row[0]
            )

            local_genes.append(
                gene
            )


            vals = np.asarray(
                row[1:],
                dtype=np.float64
            )


            if len(vals) != len(
                cells
            ):

                raise RuntimeError(
                    f"{sample}: count vector length mismatch at {gene}"
                )


            sums = np.bincount(
                codes[
                    valid
                ],
                weights=vals[
                    valid
                ],
                minlength=len(
                    categories
                )
            )


            for j, c in enumerate(
                categories
            ):

                sample_vectors[
                    c
                ].append(
                    sums[j]
                )


    if gene_order is None:

        gene_order = local_genes

    else:

        if local_genes != gene_order:

            raise RuntimeError(
                f"{sample}: gene order differs"
            )


    for j, c in enumerate(
        categories
    ):

        key = (
            sample,
            c
        )

        counts_by_column[
            key
        ] = np.asarray(
            sample_vectors[c],
            dtype=np.float64
        )


        cell_design.append({

            "sample":
                sample,

            "treatment":
                treatment,

            "celltype":
                c,

            "n_cells":
                int(
                    ncells[j]
                ),
        })


design = pd.DataFrame(
    cell_design
)


design.to_csv(
    OUT /
    "CELLTYPE_SAMPLE_DESIGN.tsv",
    sep="\t",
    index=False
)


eligibility = []


for celltype in sorted(
    design[
        "celltype"
    ].unique()
):

    d = design[
        (
            design["celltype"]
            == celltype
        )
        &
        (
            design["n_cells"]
            >= 20
        )
    ].copy()


    n_sham = int(
        (
            d["treatment"]
            == "Sham"
        ).sum()
    )

    n_d02 = int(
        (
            d["treatment"]
            == "D02"
        ).sum()
    )

    n_d14 = int(
        (
            d["treatment"]
            == "D14"
        ).sum()
    )


    eligible = (
        n_sham >= 2
        and n_d02 >= 2
        and n_d14 >= 2
    )


    strict = (
        n_sham >= 3
        and n_d02 >= 3
        and n_d14 >= 3
    )


    eligibility.append({

        "celltype":
            celltype,

        "sham_samples":
            n_sham,

        "d02_samples":
            n_d02,

        "d14_samples":
            n_d14,

        "edgeR_eligible_min2":
            eligible,

        "strict_replicated_min3":
            strict,
    })


elig = pd.DataFrame(
    eligibility
)


elig.to_csv(
    OUT /
    "CELLTYPE_DGE_ELIGIBILITY.tsv",
    sep="\t",
    index=False
)


manifest = []


for _, e in elig.iterrows():

    if not bool(
        e[
            "edgeR_eligible_min2"
        ]
    ):
        continue


    celltype = e[
        "celltype"
    ]


    d = design[
        (
            design[
                "celltype"
            ]
            == celltype
        )
        &
        (
            design[
                "n_cells"
            ]
            >= 20
        )
    ].copy()


    samples = d[
        "sample"
    ].tolist()


    matrix = pd.DataFrame({

        "gene":
            gene_order
    })


    for sample in samples:

        matrix[
            sample
        ] = counts_by_column[
            (
                sample,
                celltype
            )
        ]


    safe = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        celltype
    )


    counts_file = (
        OUT /
        f"{safe}.counts.tsv.gz"
    )

    design_file = (
        OUT /
        f"{safe}.design.tsv"
    )


    matrix.to_csv(
        counts_file,
        sep="\t",
        index=False,
        compression="gzip"
    )


    d[
        [
            "sample",
            "treatment",
            "n_cells",
        ]
    ].to_csv(
        design_file,
        sep="\t",
        index=False
    )


    manifest.append({

        "celltype":
            celltype,

        "counts_file":
            str(
                counts_file
            ),

        "design_file":
            str(
                design_file
            ),

        "n_samples":
            len(samples),

        "strict_replicated_min3":
            bool(
                e[
                    "strict_replicated_min3"
                ]
            ),
    })


pd.DataFrame(
    manifest
).to_csv(
    OUT /
    "PSEUDOBULK_MANIFEST.tsv",
    sep="\t",
    index=False
)


print()
print(
    "===== ELIGIBILITY ====="
)

print(
    elig.to_string(
        index=False
    )
)


print()
print(
    "ELIGIBLE_CELLTYPES=",
    len(manifest)
)


if len(manifest) == 0:

    raise SystemExit(
        "No cell type has sufficient biological replication"
    )

PY

}


# =============================================================================
# 05. edgeR FULL-GENE DGE
# =============================================================================

run_edger() {

  export R_LIBS_USER="$RLIB"

  Rscript - \
    "$RLIB" \
    "$OUT/gse225948_pseudobulk" \
    "$OUT/gse225948_dge" <<'RS'

args <- commandArgs(
  trailingOnly=TRUE
)

RLIB <- args[[1]]
PB <- args[[2]]
OUT <- args[[3]]

dir.create(
  RLIB,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

.libPaths(
  c(
    RLIB,
    .libPaths()
  )
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

}


if (!requireNamespace(
      "edgeR",
      quietly=TRUE
    )) {

  BiocManager::install(
    "edgeR",
    ask=FALSE,
    update=FALSE,
    lib=RLIB
  )

}


if (!requireNamespace(
      "edgeR",
      quietly=TRUE
    )) {
  stop(
    "edgeR installation failed"
  )
}


library(edgeR)


manifest <- read.delim(
  file.path(
    PB,
    "PSEUDOBULK_MANIFEST.tsv"
  ),
  stringsAsFactors=FALSE,
  check.names=FALSE
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


  # Aggregate duplicated symbols if present.
  counts <- as.matrix(C)

  storage.mode(
    counts
  ) <- "numeric"


  rownames(counts) <- gene


  if (anyDuplicated(
        rownames(counts)
      )) {

    counts <- rowsum(
      counts,
      group=rownames(counts),
      reorder=FALSE
    )

  }


  if (!all(
        D$sample
        %in%
        colnames(counts)
      )) {
    stop(
      paste(
        celltype,
        "design/count columns mismatch"
      )
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


  if (any(
        is.na(group)
      )) {
    stop(
      paste(
        celltype,
        "unexpected treatment"
      )
    )
  }


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


  if (nrow(y) < 100) {

    cat(
      "SKIP_LOW_EXPRESSED_GENES=",
      nrow(y),
      "\n",
      sep=""
    )

    next
  }


  y <- calcNormFactors(
    y,
    method="TMM"
  )


  design <- model.matrix(
    ~0 + group
  )


  colnames(design) <- levels(
    group
  )


  y <- estimateDisp(
    y,
    design
  )


  fit <- glmQLFit(
    y,
    design,
    robust=FALSE
  )


  contrasts <- list(

    D02_vs_Sham=
      makeContrasts(
        D02 - Sham,
        levels=design
      ),

    D14_vs_Sham=
      makeContrasts(
        D14 - Sham,
        levels=design
      )
  )


  safe <- gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    celltype
  )


  for (
    contrast_name
    in names(
      contrasts
    )
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

    tab$celltype <- celltype
    tab$contrast <- contrast_name


    tab <- tab[
      ,
      c(
        "gene",
        "celltype",
        "contrast",
        setdiff(
          colnames(tab),
          c(
            "gene",
            "celltype",
            "contrast"
          )
        )
      )
    ]


    outfile <- file.path(
      OUT,
      paste0(
        safe,
        ".",
        contrast_name,
        ".edgeR.tsv.gz"
      )
    )


    con <- gzfile(
      outfile,
      "wt"
    )

    write.table(
      tab,
      con,
      sep="\t",
      quote=FALSE,
      row.names=FALSE
    )

    close(con)


    tt <- tab[
      tab$gene
      %in%
      targets,
      ,
      drop=FALSE
    ]


    if (nrow(tt) > 0) {

      target_results[[
        length(
          target_results
        ) + 1
      ]] <- tt

    }


    run_rows[[
      length(
        run_rows
      ) + 1
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

      strict_replicated_min3=
        manifest$
        strict_replicated_min3[i],

      stringsAsFactors=FALSE
    )

  }

}


if (
  length(
    target_results
  ) > 0
) {

  target_master <- do.call(
    rbind,
    target_results
  )

} else {

  target_master <- data.frame()

}


if (
  nrow(
    target_master
  ) > 0
) {

  write.table(
    target_master,
    file.path(
      OUT,
      "TARGET_DGE_MASTER.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

}


if (
  length(
    run_rows
  ) > 0
) {

  run_summary <- do.call(
    rbind,
    run_rows
  )

  write.table(
    run_summary,
    file.path(
      OUT,
      "DGE_RUN_SUMMARY.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  print(
    run_summary
  )

}


cat(
  "\n===== TARGET RESULTS =====\n"
)


if (
  nrow(
    target_master
  ) > 0
) {

  target_master <- target_master[
    order(
      target_master$FDR,
      target_master$PValue
    ),
    ,
    drop=FALSE
  ]

  print(
    target_master
  )

} else {

  cat(
    "NO_TARGET_RESULTS\n"
  )

}

RS

}


# =============================================================================
# 06. TARGET INTERPRETATION MASTER
# =============================================================================

build_target_dge_master() {

"$ENV/bin/python" - \
  "$OUT/gse225948_dge" \
  "$OUT" <<'PY'

from pathlib import Path
import sys

import pandas as pd


DGE = Path(sys.argv[1])
OUT = Path(sys.argv[2])

p = (
    DGE /
    "TARGET_DGE_MASTER.tsv"
)


if not p.exists():

    print(
        "TARGET_DGE_MASTER_MISSING"
    )

    raise SystemExit(1)


d = pd.read_csv(
    p,
    sep="\t"
)


mapping = {

    "Fgf5":
        "FGF5",

    "Aldh2":
        "ALDH2",

    "Sh3pxd2a":
        "SH3PXD2A",

    "Col4a2":
        "COL4A2",

    "Col4a1":
        "COL4A1",

    "Calhm2":
        "CALHM2",

    "Neurl1":
        "NEURL1",

    "Ina":
        "INA",
}


d[
    "human_gene"
] = d[
    "gene"
].map(
    mapping
)


d[
    "significant_fdr05"
] = (
    d["FDR"]
    < 0.05
)


d[
    "direction"
] = "NEUTRAL"

d.loc[
    d["logFC"] > 0,
    "direction"
] = "UP"

d.loc[
    d["logFC"] < 0,
    "direction"
] = "DOWN"


d = d.sort_values(
    [
        "FDR",
        "PValue",
    ]
)


d.to_csv(
    OUT /
    "GSE225948_TARGET_DGE_INTERPRETATION.tsv",
    sep="\t",
    index=False
)


print(
    d[
        [
            "human_gene",
            "gene",
            "celltype",
            "contrast",
            "logFC",
            "PValue",
            "FDR",
            "direction",
            "significant_fdr05",
        ]
    ]
    .head(100)
    .to_string(
        index=False
    )
)

PY

}


# =============================================================================
# 07. HUMAN RDS SERVER RESOURCE GATE
# =============================================================================

human_rds_resource_gate() {

"$ENV/bin/python" - \
  "$DATA" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import subprocess
import sys


DATA = Path(sys.argv[1])
OUT = Path(sys.argv[2])


RDS = (
    DATA /
    "GSE256493" /
    "GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz"
)


mem_kb = 0

with open(
    "/proc/meminfo"
) as f:

    for line in f:

        if line.startswith(
            "MemTotal:"
        ):

            mem_kb = int(
                line.split()[1]
            )

            break


mem_gb = (
    mem_kb /
    1024 /
    1024
)


def r_package_available(pkg):

    p = subprocess.run(
        [
            "Rscript",
            "-e",
            (
                f'cat(requireNamespace('
                f'"{pkg}",quietly=TRUE))'
            ),
        ],
        capture_output=True,
        text=True,
    )

    return (
        p.returncode == 0
        and "TRUE"
        in p.stdout
    )


seurat = r_package_available(
    "Seurat"
)

seurat_obj = (
    r_package_available(
        "SeuratObject"
    )
)


safe_ram = (
    mem_gb >= 16
)


if (
    not RDS.exists()
):

    decision = (
        "BLOCKED_RDS_MISSING"
    )

elif (
    not safe_ram
):

    decision = (
        "DEFER_TO_HIGH_MEMORY_RUNTIME"
    )

elif not (
    seurat
    or seurat_obj
):

    decision = (
        "INSTALL_SEURATOBJECT_BEFORE_PARSE"
    )

else:

    decision = (
        "SERVER_PARSE_ALLOWED"
    )


row = {

    "rds_path":
        str(RDS),

    "rds_bytes":
        (
            RDS.stat().st_size
            if RDS.exists()
            else ""
        ),

    "server_total_ram_gib":
        mem_gb,

    "minimum_parse_gate_gib":
        16,

    "seurat_installed":
        seurat,

    "seuratobject_installed":
        seurat_obj,

    "decision":
        decision,

    "scientific_status":
        (
            "DATA_ACQUIRED_ANALYSIS_PENDING"
            if RDS.exists()
            else "DATA_PENDING"
        ),
}


with (
    OUT /
    "GSE256493_HUMAN_RDS_RESOURCE_GATE.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            row.keys()
        )
    )

    w.writeheader()
    w.writerow(row)


for k, v in row.items():

    print(
        f"{k}={v}"
    )

PY

}


# =============================================================================
# 08. READINESS
# =============================================================================

build_readiness() {

"$ENV/bin/python" - \
  "$OUT" \
  "$STATUS" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])
STATUS = Path(sys.argv[2])


step_rows = []

if STATUS.exists():

    with STATUS.open() as f:

        step_rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


def passed(step):

    for x in step_rows:

        if x[
            "step"
        ] == step:

            return (
                x["status"]
                == "PASS"
            )

    return False


g189 = passed(
    "03_repair_GSE189432"
)

pb = passed(
    "04_build_GSE225948_full_pseudobulk"
)

dge = passed(
    "05_run_edgeR_DGE"
)

target = passed(
    "06_build_target_DGE_master"
)


gate_file = (
    OUT /
    "GSE256493_HUMAN_RDS_RESOURCE_GATE.tsv"
)

human_status = (
    "UNKNOWN"
)

if gate_file.exists():

    with gate_file.open() as f:

        rr = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    if rr:

        human_status = rr[0][
            "decision"
        ]


rows = [

    (
        "PHASE9F_D",
        (
            "COMPLETE_MOUSE_FUNCTIONAL_LAYER"
            if (
                g189
                and pb
                and dge
                and target
            )
            else "PARTIAL"
        )
    ),

    (
        "GSE189432_BARCODE_REPAIR",
        (
            "COMPLETE"
            if g189
            else "REVIEW_REQUIRED"
        )
    ),

    (
        "GSE189432_INFERENCE",
        "DESCRIPTIVE_CELLTYPE_LOCALIZATION"
    ),

    (
        "GSE225948_PRIMARY_DESIGN",
        "BRAIN_W8_SHAM3_D02_4_D14_3"
    ),

    (
        "GSE225948_FULL_GENE_PSEUDOBULK",
        (
            "COMPLETE"
            if pb
            else "FAILED"
        )
    ),

    (
        "GSE225948_EDGE_R",
        (
            "COMPLETE"
            if dge
            else "FAILED"
        )
    ),

    (
        "GSE225948_TARGET_DGE",
        (
            "COMPLETE"
            if target
            else "FAILED"
        )
    ),

    (
        "GSE234052",
        "DEFER_CELL_CALLING_RAW_WHITELIST_MATRIX"
    ),

    (
        "GSE256493_HUMAN_REFERENCE",
        human_status
    ),

    (
        "HUMAN_REFERENCE_SCIENTIFIC_STATUS",
        "ACQUIRED_NOT_YET_ANALYSED"
    ),

    (
        "NEXT_STAGE",
        (
            "PHASE9F_E_HUMAN_HIGH_MEMORY_PARSE_AND_INTEGRATED_INTERPRETATION"
            if (
                g189
                and dge
            )
            else "REPAIR_FAILED_MOUSE_BRANCH_FIRST"
        )
    ),
]


with (
    OUT /
    "PHASE9F_D_READINESS.tsv"
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
  "02_correct_Phase9FC_status" \
  correct_phase9fc_status

run_step \
  "03_repair_GSE189432" \
  repair_gse189432

run_step \
  "04_build_GSE225948_full_pseudobulk" \
  build_gse225948_full_pseudobulk

run_step \
  "05_run_edgeR_DGE" \
  run_edger

run_step \
  "06_build_target_DGE_master" \
  build_target_dge_master

run_step \
  "07_human_RDS_resource_gate" \
  human_rds_resource_gate

run_step \
  "08_build_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9F-D COMPLETE"
echo "================================================================================"


echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== GSE189432 BARCODE REPAIR ====="

column -t -s $'\t' \
  "$OUT/GSE189432_BARCODE_REPAIR.tsv" \
  2>/dev/null \
  || true


echo
echo "===== GSE189432 SH3PXD2A ====="

"$ENV/bin/python" - \
  "$OUT/GSE189432_TARGET_CLUSTER_SUMMARY_R1.tsv" <<'PY'

import pandas as pd
import sys

try:

    d = pd.read_csv(
        sys.argv[1],
        sep="\t"
    )

    q = (
        d[
            d["gene"]
            == "Sh3pxd2a"
        ]
        .sort_values(
            [
                "detection_fraction",
                "cpm",
            ],
            ascending=False
        )
        .head(20)
    )

    print(
        q.to_string(
            index=False
        )
    )

except Exception as e:

    print(
        "NOT_AVAILABLE",
        type(e).__name__,
        str(e)
    )

PY


echo
echo "===== GSE225948 DGE ELIGIBILITY ====="

column -t -s $'\t' \
  "$OUT/gse225948_pseudobulk/CELLTYPE_DGE_ELIGIBILITY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== TARGET DGE ====="

column -t -s $'\t' \
  "$OUT/GSE225948_TARGET_DGE_INTERPRETATION.tsv" \
  2>/dev/null \
  | head -100 \
  || true


echo
echo "===== HUMAN RDS RESOURCE GATE ====="

column -t -s $'\t' \
  "$OUT/GSE256493_HUMAN_RDS_RESOURCE_GATE.tsv" \
  2>/dev/null \
  || true


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9F_D_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE9F_D_READINESS.tsv"


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
