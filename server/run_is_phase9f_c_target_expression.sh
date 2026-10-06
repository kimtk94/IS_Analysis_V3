#!/usr/bin/env bash

ROOT="/srv/is-analysis"

P9FB="$ROOT/results/is/stage5_functional/phase9f_b_format_audit"

DATA="$ROOT/data/is/celltype"
OUT="$ROOT/results/is/stage5_functional/phase9f_c_target_expression"

ENV="$ROOT/envs/is_singlecell"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9f_c/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

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

  echo "===== PHASE9F-B ====="

  if [[ ! -s "$P9FB/STEP_STATUS.tsv" ]]; then
    echo "P9FB_STATUS_MISSING"
    return 1
  fi

  cat "$P9FB/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P9FB/STEP_STATUS.tsv"
  )"

  echo
  echo "PHASE9FB_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    return 1
  fi


  if [[ ! -x "$ENV/bin/python" ]]; then
    echo "SINGLECELL_ENV_MISSING"
    return 1
  fi


  "$ENV/bin/python" - <<'PY'
import numpy
import pandas
import scipy
import anndata
import scanpy

print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
print("anndata", anndata.__version__)
print("scanpy", scanpy.__version__)
PY


  echo
  echo "===== MEMORY ====="
  free -h

  echo
  echo "===== DISK ====="
  df -h "$ROOT"

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. REPAIR GSE234052 INTERPRETATION
# =============================================================================

repair_gse234052_status() {

"$ENV/bin/python" - \
  "$P9FB/GSE234052_BARCODE_SUFFIX_AUDIT.tsv" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


INPUT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

with INPUT.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


EXPECTED_WHITELIST = 6794880

counts = [
    int(x["n_cells"])
    for x in rows
]


all_whitelist = (
    len(rows) == 6
    and all(
        x == EXPECTED_WHITELIST
        for x in counts
    )
)


total_columns = sum(
    counts
)


status = (
    "RAW_BARCODE_WHITELIST_MATRIX_CELL_CALLING_REQUIRED"
    if all_whitelist
    else
    "REVIEW_REQUIRED"
)


result = {

    "dataset":
        "GSE234052",

    "matrix_columns":
        total_columns,

    "conditions":
        len(rows),

    "columns_per_condition":
        EXPECTED_WHITELIST
        if all_whitelist
        else "",

    "published_final_cells":
        23675,

    "published_mean_genes_per_cell":
        2020,

    "interpretation":
        status,

    "do_not_interpret_matrix_columns_as_cells":
        "YES",

    "next_action":
        "CELL_CALLING_QC_AND_RECLUSTERING",

    "source":
        "Buizza_et_al_Transl_Stroke_Res_2024_DOI_10.1007/s12975-023-01169-x",
}


with (
    OUT /
    "GSE234052_MATRIX_INTERPRETATION.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            result.keys()
        )
    )

    w.writeheader()
    w.writerow(result)


for k, v in result.items():

    print(
        f"{k}={v}"
    )


if not all_whitelist:

    raise SystemExit(
        "Expected whitelist pattern not detected"
    )

PY

}


# =============================================================================
# 03. GSE189432 ANNOTATED MYELOID TARGET EXPRESSION
# =============================================================================

analyse_gse189432() {

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

ANNOT = (
    ROOT /
    "GSE189432_annotations.csv.gz"
)


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


print(
    "ANNOTATED_CELLS=",
    len(ann)
)

print(
    "CLUSTERS=",
    sorted(
        ann["cluster"]
        .dropna()
        .astype(str)
        .unique()
    )
)


def read_lines_gz(path):

    with gzip.open(
        path,
        "rt",
        errors="replace"
    ) as f:

        return [
            x.rstrip("\n")
            for x in f
            if x.strip()
        ]


def normalize_barcode(x):

    x = str(x).strip()

    # Do not remove internal sample identifiers.
    # Only CellRanger terminal numeric suffix.
    return re.sub(
        r"-\d+$",
        "",
        x
    )


matrices = sorted(
    ROOT.rglob(
        "*_matrix.mtx.gz"
    )
)


expression_rows = []
match_rows = []


for matrix_path in matrices:

    prefix = matrix_path.name[
        :-len(
            "_matrix.mtx.gz"
        )
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

        print(
            "MISSING_TRIPLET",
            prefix
        )
        continue


    barcodes = read_lines_gz(
        barcode_path
    )


    feature_rows = []

    with gzip.open(
        feature_path,
        "rt",
        errors="replace"
    ) as f:

        for line in f:

            x = line.rstrip(
                "\n"
            ).split("\t")

            feature_rows.append(x)


    genes = []

    for x in feature_rows:

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


    gene_to_idx = {}

    for i, g in enumerate(
        genes
    ):

        if (
            g in TARGETS
            and g not in gene_to_idx
        ):
            gene_to_idx[g] = i


    a = ann[
        ann["sample"]
        .astype(str)
        == sample
    ].copy()


    raw_exact = {
        str(x): i
        for i, x in enumerate(
            barcodes
        )
    }


    raw_norm = {}

    duplicate_norm = set()

    for i, b in enumerate(
        barcodes
    ):

        n = normalize_barcode(b)

        if n in raw_norm:
            duplicate_norm.add(n)

        raw_norm[n] = i


    for x in duplicate_norm:
        raw_norm.pop(
            x,
            None
        )


    exact_matches = sum(
        str(x) in raw_exact

        for x in a["barcode"]
    )


    norm_matches = sum(
        normalize_barcode(x)
        in raw_norm

        for x in a["barcode"]
    )


    if len(a) > 0:

        exact_rate = (
            exact_matches /
            len(a)
        )

        norm_rate = (
            norm_matches /
            len(a)
        )

    else:

        exact_rate = 0.0
        norm_rate = 0.0


    if exact_rate >= 0.95:

        strategy = "EXACT"

        mapper = raw_exact

        a["matrix_col"] = [
            mapper.get(
                str(x),
                -1
            )
            for x in a[
                "barcode"
            ]
        ]


    elif norm_rate >= 0.95:

        strategy = (
            "STRIP_TERMINAL_CELLRANGER_SUFFIX"
        )

        mapper = raw_norm

        a["matrix_col"] = [
            mapper.get(
                normalize_barcode(x),
                -1
            )
            for x in a[
                "barcode"
            ]
        ]


    else:

        strategy = (
            "MATCH_FAILED"
        )

        a["matrix_col"] = -1


    a = a[
        a["matrix_col"] >= 0
    ].copy()


    match_rows.append({

        "sample":
            sample,

        "raw_cells":
            len(barcodes),

        "annotated_cells":
            len(
                ann[
                    ann["sample"]
                    .astype(str)
                    == sample
                ]
            ),

        "matched_cells":
            len(a),

        "exact_match_rate":
            exact_rate,

        "normalized_match_rate":
            norm_rate,

        "strategy":
            strategy,
    })


    print(
        sample,
        "raw=",
        len(barcodes),
        "annotated=",
        len(
            ann[
                ann["sample"]
                .astype(str)
                == sample
            ]
        ),
        "matched=",
        len(a),
        "strategy=",
        strategy,
    )


    if strategy == "MATCH_FAILED":

        continue


    print(
        "Loading",
        matrix_path.name
    )


    with gzip.open(
        matrix_path,
        "rb"
    ) as fh:

        X = mmread(
            fh
        ).tocsr()


    if X.shape[
        1
    ] != len(
        barcodes
    ):

        raise RuntimeError(
            f"{sample}: matrix/barcode mismatch"
        )


    total_umi = np.asarray(
        X.sum(
            axis=0
        )
    ).ravel()


    cols = a[
        "matrix_col"
    ].astype(int).to_numpy()


    clusters = (
        a["cluster"]
        .astype(str)
        .to_numpy()
    )


    unique_clusters = sorted(
        set(clusters)
    )


    for cluster in unique_clusters:

        sel = np.where(
            clusters
            == cluster
        )[0]

        matrix_cols = cols[
            sel
        ]

        lib_size = float(
            total_umi[
                matrix_cols
            ].sum()
        )


        for gene in TARGETS:

            if gene not in gene_to_idx:

                expression_rows.append({

                    "sample":
                        sample,

                    "cluster":
                        cluster,

                    "gene":
                        gene,

                    "n_cells":
                        len(
                            matrix_cols
                        ),

                    "sum_counts":
                        "",

                    "mean_counts":
                        "",

                    "detected_cells":
                        "",

                    "detection_fraction":
                        "",

                    "library_size":
                        lib_size,

                    "cpm":
                        "",

                    "gene_present":
                        False,
                })

                continue


            gi = gene_to_idx[
                gene
            ]

            vals = (
                X[
                    gi,
                    matrix_cols
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


            expression_rows.append({

                "sample":
                    sample,

                "cluster":
                    cluster,

                "gene":
                    gene,

                "n_cells":
                    len(
                        matrix_cols
                    ),

                "sum_counts":
                    s,

                "mean_counts":
                    float(
                        vals.mean()
                    )
                    if len(vals)
                    else 0,

                "detected_cells":
                    detected,

                "detection_fraction":
                    (
                        detected
                        /
                        len(vals)
                        if len(vals)
                        else 0
                    ),

                "library_size":
                    lib_size,

                "cpm":
                    (
                        s
                        /
                        lib_size
                        * 1e6
                        if lib_size > 0
                        else 0
                    ),

                "gene_present":
                    True,
            })


with (
    OUT /
    "GSE189432_BARCODE_MATCH.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            match_rows[0].keys()
        )
    )

    w.writeheader()
    w.writerows(
        match_rows
    )


with (
    OUT /
    "GSE189432_TARGET_CLUSTER_EXPRESSION.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            expression_rows[0].keys()
        )
    )

    w.writeheader()
    w.writerows(
        expression_rows
    )


# ------------------------------------------------------------------
# Cross-sample localization summary.
# Weighted by cells for descriptive localization only.
# ------------------------------------------------------------------

df = pd.DataFrame(
    expression_rows
)

df = df[
    df[
        "gene_present"
    ] == True
].copy()


summary = (
    df.groupby(
        [
            "cluster",
            "gene",
        ],
        dropna=False
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
        library_size=(
            "library_size",
            "sum"
        ),
        detected_cells=(
            "detected_cells",
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
    "GSE189432_TARGET_CLUSTER_SUMMARY.tsv",
    sep="\t",
    index=False
)


print()
print(
    "===== SH3PXD2A TOP CLUSTERS ====="
)


q = (
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
)


print(
    q.to_string(
        index=False
    )
)

PY

}


# =============================================================================
# 04. GSE225948 SAMPLE-AWARE TARGET PSEUDOBULK
# =============================================================================

analyse_gse225948() {

"$ENV/bin/python" - \
  "$DATA/GSE225948/processed" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import re
import sys

import numpy as np
import pandas as pd


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])


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


metadata_files = sorted(
    ROOT.glob(
        "*_metadata.csv.gz"
    )
)


design_rows = []
match_rows = []
pb_rows = []


def unique_value(
    df,
    col,
):

    if col not in df.columns:
        return ""

    vals = (
        df[col]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    return ";".join(
        sorted(vals)
    )


def parse_day(text):

    m = re.search(
        r"day\s*([0-9]+)",
        str(text),
        flags=re.I
    )

    return (
        m.group(1)
        if m
        else ""
    )


for meta_path in metadata_files:

    sample = meta_path.name.replace(
        "_metadata.csv.gz",
        ""
    )

    counts_path = (
        ROOT /
        f"{sample}_counts.csv.gz"
    )


    if not counts_path.exists():
        continue


    meta = pd.read_csv(
        meta_path,
        compression="gzip",
        index_col=0
    )

    meta.index = (
        meta.index
        .astype(str)
    )


    sample_description = unique_value(
        meta,
        "Sample_description"
    )

    treatment = unique_value(
        meta,
        "treatment"
    )

    replicate = unique_value(
        meta,
        "Replicate"
    )

    tissue = unique_value(
        meta,
        "tissue"
    )

    age = unique_value(
        meta,
        "age"
    )

    sex = unique_value(
        meta,
        "sex"
    )


    day = parse_day(
        sample_description
    )


    design_rows.append({

        "sample":
            sample,

        "n_cells":
            len(meta),

        "tissue":
            tissue,

        "age":
            age,

        "sex":
            sex,

        "treatment":
            treatment,

        "replicate":
            replicate,

        "day":
            day,

        "sample_description":
            sample_description,

        "sub_celltypes":
            meta[
                "sub.celltype"
            ].nunique()
            if "sub.celltype"
            in meta.columns
            else "",

        "parents":
            meta[
                "parent"
            ].nunique()
            if "parent"
            in meta.columns
            else "",
    })


    with gzip.open(
        counts_path,
        "rt",
        errors="replace"
    ) as fh:

        reader = csv.reader(
            fh
        )

        header = next(
            reader
        )

        count_cells = [
            str(x)
            for x in header[
                1:
            ]
        ]


        exact = sum(
            x in meta.index
            for x in count_cells
        )


        match_rate = (
            exact
            /
            len(count_cells)
            if count_cells
            else 0
        )


        match_rows.append({

            "sample":
                sample,

            "count_cells":
                len(
                    count_cells
                ),

            "metadata_cells":
                len(meta),

            "exact_matches":
                exact,

            "match_rate":
                match_rate,
        })


        print(
            sample,
            "cells=",
            len(
                count_cells
            ),
            "match=",
            match_rate,
            "treatment=",
            treatment,
            "day=",
            day,
            "rep=",
            replicate,
        )


        if match_rate < 0.95:

            print(
                "SKIP_LOW_CELLID_MATCH",
                sample
            )

            continue


        meta_aligned = meta.reindex(
            count_cells
        )


        target_vectors = {}


        for row in reader:

            if not row:
                continue

            gene = str(
                row[0]
            )

            if gene not in TARGETS:
                continue


            vals = np.fromiter(
                (
                    float(x)
                    if x not in {
                        "",
                        "NA",
                        "NaN",
                    }
                    else 0.0
                    for x in row[
                        1:
                    ]
                ),
                dtype=float,
                count=len(
                    count_cells
                ),
            )


            target_vectors[
                gene
            ] = vals


    if "nCount_RNA" not in meta_aligned.columns:

        raise RuntimeError(
            f"{sample}: nCount_RNA unavailable"
        )


    lib = pd.to_numeric(
        meta_aligned[
            "nCount_RNA"
        ],
        errors="coerce"
    ).fillna(
        0
    ).to_numpy(
        dtype=float
    )


    for level in [
        "parent",
        "sub.celltype",
    ]:

        if level not in meta_aligned.columns:
            continue


        groups = (
            meta_aligned[
                level
            ]
            .fillna(
                "NA"
            )
            .astype(str)
            .to_numpy()
        )


        for group in sorted(
            set(groups)
        ):

            idx = np.where(
                groups == group
            )[0]


            if len(idx) == 0:
                continue


            group_lib = float(
                lib[
                    idx
                ].sum()
            )


            for gene in TARGETS:

                vals = target_vectors.get(
                    gene
                )


                if vals is None:

                    pb_rows.append({

                        "sample":
                            sample,

                        "tissue":
                            tissue,

                        "age":
                            age,

                        "sex":
                            sex,

                        "treatment":
                            treatment,

                        "replicate":
                            replicate,

                        "day":
                            day,

                        "sample_description":
                            sample_description,

                        "celltype_level":
                            level,

                        "celltype":
                            group,

                        "gene":
                            gene,

                        "gene_present":
                            False,

                        "n_cells":
                            len(idx),

                        "library_size":
                            group_lib,

                        "target_sum":
                            "",

                        "target_mean":
                            "",

                        "detected_cells":
                            "",

                        "detection_fraction":
                            "",

                        "cpm":
                            "",
                    })

                    continue


                x = vals[
                    idx
                ]


                s = float(
                    x.sum()
                )


                detected = int(
                    np.sum(
                        x > 0
                    )
                )


                pb_rows.append({

                    "sample":
                        sample,

                    "tissue":
                        tissue,

                    "age":
                        age,

                    "sex":
                        sex,

                    "treatment":
                        treatment,

                    "replicate":
                        replicate,

                    "day":
                        day,

                    "sample_description":
                        sample_description,

                    "celltype_level":
                        level,

                    "celltype":
                        group,

                    "gene":
                        gene,

                    "gene_present":
                        True,

                    "n_cells":
                        len(idx),

                    "library_size":
                        group_lib,

                    "target_sum":
                        s,

                    "target_mean":
                        float(
                            x.mean()
                        ),

                    "detected_cells":
                        detected,

                    "detection_fraction":
                        (
                            detected
                            /
                            len(x)
                        ),

                    "cpm":
                        (
                            s
                            /
                            group_lib
                            * 1e6
                            if group_lib > 0
                            else 0
                        ),
                })


pd.DataFrame(
    design_rows
).to_csv(
    OUT /
    "GSE225948_SAMPLE_DESIGN.tsv",
    sep="\t",
    index=False
)


pd.DataFrame(
    match_rows
).to_csv(
    OUT /
    "GSE225948_CELLID_MATCH.tsv",
    sep="\t",
    index=False
)


pb = pd.DataFrame(
    pb_rows
)


pb.to_csv(
    OUT /
    "GSE225948_TARGET_PSEUDOBULK.tsv",
    sep="\t",
    index=False
)


# ------------------------------------------------------------------
# Replicate/design availability summary.
# ------------------------------------------------------------------

design = pd.DataFrame(
    design_rows
)


if not design.empty:

    summary = (
        design.groupby(
            [
                "tissue",
                "age",
                "day",
                "treatment",
            ],
            dropna=False
        )
        .agg(
            biological_samples=(
                "sample",
                "nunique"
            ),
            samples=(
                "sample",
                lambda x:
                    ";".join(
                        sorted(
                            set(x)
                        )
                    )
            ),
        )
        .reset_index()
    )


    summary.to_csv(
        OUT /
        "GSE225948_DESIGN_SUMMARY.tsv",
        sep="\t",
        index=False
    )


    print()
    print(
        "===== DESIGN SUMMARY ====="
    )

    print(
        summary.to_string(
            index=False
        )
    )


print()
print(
    "===== SH3PXD2A TOP PSEUDOBULK ROWS ====="
)


if not pb.empty:

    q = (
        pb[
            (
                pb["gene"]
                == "Sh3pxd2a"
            )
            &
            (
                pb[
                    "celltype_level"
                ]
                == "sub.celltype"
            )
        ]
        .sort_values(
            [
                "detection_fraction",
                "cpm",
            ],
            ascending=False
        )
        .head(40)
    )


    print(
        q[
            [
                "sample",
                "tissue",
                "age",
                "treatment",
                "day",
                "celltype",
                "n_cells",
                "detection_fraction",
                "cpm",
            ]
        ].to_string(
            index=False
        )
    )

PY

}


# =============================================================================
# 05. ACQUIRE AUTHOR-ANNOTATED HUMAN VASCULAR REFERENCE
# =============================================================================

download_human_annotated_reference() {

  D="$DATA/GSE256493"

  mkdir -p "$D"

  FN="GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz"

  URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256493/suppl/$FN"

  DEST="$D/$FN"

  echo "URL=$URL"
  echo "DEST=$DEST"


  if [[ ! -s "$DEST" ]]; then

    curl \
      -L \
      --fail \
      --retry 5 \
      --retry-delay 5 \
      --connect-timeout 30 \
      -C - \
      -o "$DEST.part" \
      "$URL"

    RC=$?

    if [[ $RC -ne 0 ]]; then
      echo "DOWNLOAD_FAIL=$RC"
      return 1
    fi

    mv \
      "$DEST.part" \
      "$DEST"

  else

    echo "ALREADY_PRESENT"
  fi


  echo
  echo "===== SIZE ====="

  ls -lh "$DEST"

  stat -c \
    'bytes=%s' \
    "$DEST"


  echo
  echo "===== GZIP TEST ====="

  gzip -t "$DEST"

  RC=$?

  if [[ $RC -ne 0 ]]; then
    echo "GZIP_TEST_FAIL"
    return 1
  fi


  sha256sum \
    "$DEST" \
    | tee \
      "$OUT/GSE256493_PRIMARY_RDS_SHA256.txt"


  echo
  echo "===== R ENV AUDIT ====="

  {
    echo "Rscript=$(command -v Rscript 2>/dev/null || true)"

    if command -v Rscript >/dev/null 2>&1; then

      Rscript - <<'RS'
cat(
  "R=",
  as.character(
    getRversion()
  ),
  "\n",
  sep=""
)

cat(
  "Seurat=",
  requireNamespace(
    "Seurat",
    quietly=TRUE
  ),
  "\n",
  sep=""
)

cat(
  "SeuratObject=",
  requireNamespace(
    "SeuratObject",
    quietly=TRUE
  ),
  "\n",
  sep=""
)
RS

    fi

    free -h

  } | tee \
      "$OUT/GSE256493_R_PARSE_PREFLIGHT.txt"


  echo "GSE256493_PRIMARY_RDS=PASS"
}


# =============================================================================
# 06. BUILD INTEGRATED TARGET EVIDENCE
# =============================================================================

build_target_master() {

"$ENV/bin/python" - \
  "$OUT" \
  "$DATA" <<'PY'

from pathlib import Path
import csv
import sys

import pandas as pd


OUT = Path(sys.argv[1])
DATA = Path(sys.argv[2])


g189 = (
    OUT /
    "GSE189432_TARGET_CLUSTER_SUMMARY.tsv"
)

g225 = (
    OUT /
    "GSE225948_TARGET_PSEUDOBULK.tsv"
)


rows = []


if g189.exists():

    x = pd.read_csv(
        g189,
        sep="\t"
    )

    for gene in [
        "Sh3pxd2a",
        "Fgf5",
        "Aldh2",
        "Col4a2",
        "Col4a1",
    ]:

        q = x[
            x["gene"] == gene
        ].copy()

        if q.empty:
            continue

        q = q.sort_values(
            [
                "detection_fraction",
                "cpm",
            ],
            ascending=False
        )

        top = q.iloc[0]

        rows.append({

            "human_gene":
                {
                    "Sh3pxd2a":
                        "SH3PXD2A",
                    "Fgf5":
                        "FGF5",
                    "Aldh2":
                        "ALDH2",
                    "Col4a2":
                        "COL4A2",
                    "Col4a1":
                        "COL4A1",
                }[
                    gene
                ],

            "dataset":
                "GSE189432",

            "species":
                "Mouse",

            "evidence_type":
                "STROKE_MYeloid_CELLTYPE_LOCALIZATION",

            "top_celltype":
                top["cluster"],

            "detection_fraction":
                top[
                    "detection_fraction"
                ],

            "cpm":
                top["cpm"],

            "interpretation":
                "DESCRIPTIVE_CELLTYPE_EVIDENCE",
        })


if g225.exists():

    x = pd.read_csv(
        g225,
        sep="\t"
    )


    x = x[
        (
            x[
                "celltype_level"
            ]
            == "sub.celltype"
        )
        &
        (
            x[
                "gene_present"
            ]
            .astype(str)
            .str.lower()
            == "true"
        )
    ].copy()


    for gene in [
        "Sh3pxd2a",
        "Fgf5",
        "Aldh2",
        "Col4a2",
        "Col4a1",
    ]:

        q = x[
            x["gene"]
            == gene
        ].copy()

        if q.empty:
            continue


        # Require >=20 cells to avoid tiny-category ranking.
        q = q[
            q[
                "n_cells"
            ] >= 20
        ]


        if q.empty:
            continue


        q = q.sort_values(
            [
                "detection_fraction",
                "cpm",
            ],
            ascending=False
        )


        top = q.iloc[0]


        rows.append({

            "human_gene":
                {
                    "Sh3pxd2a":
                        "SH3PXD2A",
                    "Fgf5":
                        "FGF5",
                    "Aldh2":
                        "ALDH2",
                    "Col4a2":
                        "COL4A2",
                    "Col4a1":
                        "COL4A1",
                }[
                    gene
                ],

            "dataset":
                "GSE225948",

            "species":
                "Mouse",

            "evidence_type":
                "STROKE_SAMPLE_AWARE_TARGET_EXPRESSION",

            "top_celltype":
                top[
                    "celltype"
                ],

            "detection_fraction":
                top[
                    "detection_fraction"
                ],

            "cpm":
                top["cpm"],

            "interpretation":
                "SAMPLE_AWARE_DESCRIPTIVE_PSEUDOBULK",
        })


# Human object only status; no unverified interpretation.
human = (
    DATA /
    "GSE256493" /
    "GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz"
)


rows.append({

    "human_gene":
        "FGF5;SH3PXD2A;COL4A2;COL4A1;ALDH2",

    "dataset":
        "GSE256493",

    "species":
        "Human",

    "evidence_type":
        "AUTHOR_ANNOTATED_VASCULAR_REFERENCE",

    "top_celltype":
        "",

    "detection_fraction":
        "",

    "cpm":
        "",

    "interpretation":
        (
            "RDS_ACQUIRED_PARSE_PENDING"
            if human.exists()
            else "RDS_NOT_AVAILABLE"
        ),
})


with (
    OUT /
    "CELLTYPE_TARGET_EVIDENCE_MASTER.tsv"
).open(
    "w",
    newline=""
) as f:

    fields = [
        "human_gene",
        "dataset",
        "species",
        "evidence_type",
        "top_celltype",
        "detection_fraction",
        "cpm",
        "interpretation",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(
        rows
    )


for x in rows:

    print(
        x
    )

PY

}


# =============================================================================
# 07. READINESS
# =============================================================================

build_readiness() {

"$ENV/bin/python" - \
  "$OUT" \
  "$DATA" <<'PY'

from pathlib import Path
import csv
import sys

import pandas as pd


OUT = Path(sys.argv[1])
DATA = Path(sys.argv[2])


# GSE189 barcode match
p189 = (
    OUT /
    "GSE189432_BARCODE_MATCH.tsv"
)

g189_ok = False

if p189.exists():

    x = pd.read_csv(
        p189,
        sep="\t"
    )

    if not x.empty:

        g189_ok = (
            x[
                "matched_cells"
            ].sum()
            > 30000
            and
            (
                x[
                    "matched_cells"
                ].sum()
                /
                x[
                    "annotated_cells"
                ].sum()
            )
            >= 0.95
        )


# GSE225 ID matching
p225 = (
    OUT /
    "GSE225948_CELLID_MATCH.tsv"
)

g225_ok = False

if p225.exists():

    x = pd.read_csv(
        p225,
        sep="\t"
    )

    if not x.empty:

        g225_ok = bool(
            (
                x[
                    "match_rate"
                ]
                >= 0.95
            ).all()
        )


human = (
    DATA /
    "GSE256493" /
    "GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz"
)


rows = [

    (
        "PHASE9F_C",
        "COMPUTATIONAL_COMPLETE"
    ),

    (
        "GSE189432_ANNOTATED_TARGET_EXPRESSION",
        (
            "COMPLETE"
            if g189_ok
            else "REVIEW_REQUIRED"
        )
    ),

    (
        "GSE189432_INFERENCE_LEVEL",
        "DESCRIPTIVE_CELLTYPE_DISEASE_CONTEXT"
    ),

    (
        "GSE225948_SAMPLE_AWARE_TARGET_PSEUDOBULK",
        (
            "COMPLETE"
            if g225_ok
            else "REVIEW_REQUIRED"
        )
    ),

    (
        "GSE225948_FULL_DGE",
        "NOT_YET_RUN_REQUIRE_DESIGN_CONFIRMATION"
    ),

    (
        "GSE234052_MATRIX",
        "RAW_BARCODE_WHITELIST_CELL_CALLING_REQUIRED"
    ),

    (
        "GSE234052_MATRIX_COLUMNS_ARE_CELLS",
        "NO"
    ),

    (
        "GSE234052_PUBLISHED_FINAL_CELLS",
        "23675"
    ),

    (
        "GSE256493_PRIMARY_HUMAN_REFERENCE",
        (
            "ACQUIRED_GZIP_VALIDATED"
            if human.exists()
            else "DOWNLOAD_FAILED"
        )
    ),

    (
        "HUMAN_RDS_PARSE",
        "PENDING_MEMORY_AND_SEURAT_PREFLIGHT"
    ),

    (
        "SH3PXD2A_BRANCH",
        (
            "MOUSE_CELLTYPE_EVIDENCE_READY"
            if (
                g189_ok
                and g225_ok
            )
            else "REVIEW"
        )
    ),

    (
        "COL4A2_COL4A1_BRANCH",
        "GSE234052_REQUIRES_CELL_CALLING_HUMAN_RDS_PENDING_PARSE"
    ),

    (
        "NEXT_STAGE",
        "PHASE9F_D_HUMAN_RDS_PARSE_AND_VALIDATED_PSEUDOBULK"
    ),
]


with (
    OUT /
    "PHASE9F_C_READINESS.tsv"
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
  "02_repair_GSE234052_status" \
  repair_gse234052_status

run_step \
  "03_GSE189432_target_expression" \
  analyse_gse189432

run_step \
  "04_GSE225948_target_pseudobulk" \
  analyse_gse225948

run_step \
  "05_download_GSE256493_human_reference" \
  download_human_annotated_reference

run_step \
  "06_build_target_master" \
  build_target_master

run_step \
  "07_build_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9F-C COMPLETE"
echo "================================================================================"


echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== GSE234052 CORRECTION ====="

cat \
  "$OUT/GSE234052_MATRIX_INTERPRETATION.tsv" \
  2>/dev/null \
  || true


echo
echo "===== GSE189432 BARCODE MATCH ====="

column -t -s $'\t' \
  "$OUT/GSE189432_BARCODE_MATCH.tsv" \
  2>/dev/null \
  || true


echo
echo "===== GSE189432 SH3PXD2A TOP CLUSTERS ====="

"$ENV/bin/python" - \
  "$OUT/GSE189432_TARGET_CLUSTER_SUMMARY.tsv" <<'PY'
import pandas as pd
import sys

p = sys.argv[1]

try:
    d = pd.read_csv(
        p,
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
        .head(15)
    )

    print(
        q.to_string(
            index=False
        )
    )

except Exception as e:
    print(
        "SUMMARY_ERROR",
        type(e).__name__,
        str(e)
    )
PY


echo
echo "===== GSE225948 DESIGN ====="

column -t -s $'\t' \
  "$OUT/GSE225948_DESIGN_SUMMARY.tsv" \
  2>/dev/null \
  || true


echo
echo "===== GSE225948 CELL ID MATCH ====="

column -t -s $'\t' \
  "$OUT/GSE225948_CELLID_MATCH.tsv" \
  2>/dev/null \
  || true


echo
echo "===== HUMAN RDS ====="

ls -lh \
  "$DATA/GSE256493/"*.rds.gz \
  2>/dev/null \
  || true

cat \
  "$OUT/GSE256493_R_PARSE_PREFLIGHT.txt" \
  2>/dev/null \
  || true


echo
echo "===== TARGET EVIDENCE MASTER ====="

column -t -s $'\t' \
  "$OUT/CELLTYPE_TARGET_EVIDENCE_MASTER.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/CELLTYPE_TARGET_EVIDENCE_MASTER.tsv"


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9F_C_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE9F_C_READINESS.tsv"


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"

