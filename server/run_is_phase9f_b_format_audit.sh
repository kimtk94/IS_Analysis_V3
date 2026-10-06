#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

P9FA="$ROOT/results/is/stage5_functional/phase9f_a_acquisition"

DATA="$ROOT/data/is/celltype"
OUT="$ROOT/results/is/stage5_functional/phase9f_b_format_audit"

ENV="$ROOT/envs/is_singlecell"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9f_b/$RUN_ID"
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

  echo "===== PHASE9F-A ====="

  if [[ ! -s "$P9FA/STEP_STATUS.tsv" ]]; then
    echo "PHASE9FA_STATUS_MISSING"
    return 1
  fi

  cat "$P9FA/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P9FA/STEP_STATUS.tsv"
  )"

  echo
  echo "PHASE9FA_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    return 1
  fi


  echo
  echo "===== ENV ====="

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

for x in [
    numpy,
    pandas,
    scipy,
    anndata,
    scanpy,
]:
    print(
        x.__name__,
        getattr(x, "__version__", "")
    )
PY

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. 10X TRIPLET FORMAT AUDIT
# =============================================================================

audit_10x_triplets() {

"$ENV/bin/python" - \
  "$DATA" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import re
import sys


DATA = Path(sys.argv[1])
OUT = Path(sys.argv[2])


datasets = [
    "GSE189432",
    "GSE256490",
]


TARGETS = {
    "GSE189432": [
        "Fgf5",
        "Aldh2",
        "Sh3pxd2a",
        "Col4a2",
        "Col4a1",
        "Calhm2",
        "Neurl1",
        "Ina",
    ],

    "GSE256490": [
        "FGF5",
        "ALDH2",
        "SH3PXD2A",
        "COL4A2",
        "COL4A1",
        "CALHM2",
        "NEURL1",
        "C4orf22",
        "INA",
    ],
}


def open_text(p):

    if p.name.endswith(".gz"):
        return gzip.open(
            p,
            "rt",
            errors="replace"
        )

    return p.open(
        "r",
        errors="replace"
    )


def prefix_from_name(name):

    patterns = [
        "_matrix.mtx.gz",
        "_matrix.mtx",
        "_barcodes.tsv.gz",
        "_barcodes.tsv",
        "_features.tsv.gz",
        "_features.tsv",
        "_genes.tsv.gz",
        "_genes.tsv",
    ]

    for x in patterns:

        if name.endswith(x):
            return name[:-len(x)]

    return None


def matrix_dims(p):

    with open_text(p) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith("%"):
                continue

            x = line.split()

            if len(x) >= 3:

                return (
                    int(x[0]),
                    int(x[1]),
                    int(x[2]),
                )

    raise RuntimeError(
        f"Cannot resolve MatrixMarket dims: {p}"
    )


def count_lines(p):

    n = 0

    with open_text(p) as f:

        for line in f:

            if line.strip():
                n += 1

    return n


def gene_positions(p, targets):

    found = {}

    with open_text(p) as f:

        for i, line in enumerate(
            f,
            start=1
        ):

            x = line.rstrip(
                "\n"
            ).split("\t")

            for token in x:

                if token in targets:

                    found[
                        token
                    ] = i

    return found


all_rows = []
gene_rows = []


for accession in datasets:

    root = DATA / accession

    files = [
        p
        for p in root.rglob("*")
        if p.is_file()
    ]

    groups = {}


    for p in files:

        prefix = prefix_from_name(
            p.name
        )

        if prefix is None:
            continue

        g = groups.setdefault(
            prefix,
            {}
        )


        if (
            "_matrix.mtx"
            in p.name
        ):
            g["matrix"] = p

        elif (
            "_barcodes.tsv"
            in p.name
        ):
            g["barcodes"] = p

        elif (
            "_features.tsv"
            in p.name
            or "_genes.tsv"
            in p.name
        ):
            g["features"] = p


    print(
        accession,
        "TRIPLETS=",
        len(groups)
    )


    for prefix, g in sorted(
        groups.items()
    ):

        matrix = g.get(
            "matrix"
        )

        features = g.get(
            "features"
        )

        barcodes = g.get(
            "barcodes"
        )


        complete = all([
            matrix,
            features,
            barcodes,
        ])


        row = {
            "accession":
                accession,

            "sample_prefix":
                prefix,

            "complete_triplet":
                complete,

            "matrix_path":
                str(matrix or ""),

            "features_path":
                str(features or ""),

            "barcodes_path":
                str(barcodes or ""),

            "matrix_genes":
                "",

            "matrix_cells":
                "",

            "matrix_nnz":
                "",

            "feature_rows":
                "",

            "barcode_rows":
                "",

            "dimension_status":
                "",
        }


        if not complete:

            row[
                "dimension_status"
            ] = "INCOMPLETE_TRIPLET"

            all_rows.append(row)
            continue


        nr, nc, nnz = matrix_dims(
            matrix
        )

        nf = count_lines(
            features
        )

        nb = count_lines(
            barcodes
        )


        row.update({

            "matrix_genes":
                nr,

            "matrix_cells":
                nc,

            "matrix_nnz":
                nnz,

            "feature_rows":
                nf,

            "barcode_rows":
                nb,

            "dimension_status":
                (
                    "PASS"
                    if (
                        nr == nf
                        and nc == nb
                    )
                    else "MISMATCH"
                ),
        })


        all_rows.append(row)


        gp = gene_positions(
            features,
            TARGETS[accession],
        )


        for gene in TARGETS[
            accession
        ]:

            gene_rows.append({

                "accession":
                    accession,

                "sample_prefix":
                    prefix,

                "gene":
                    gene,

                "present":
                    gene in gp,

                "matrix_row_1based":
                    gp.get(
                        gene,
                        ""
                    ),
            })


        print(
            accession,
            prefix,
            "genes=",
            nr,
            "cells=",
            nc,
            "nnz=",
            nnz,
            row[
                "dimension_status"
            ],
            "targets=",
            ",".join(
                sorted(
                    gp.keys()
                )
            ),
        )


with (
    OUT /
    "TENX_TRIPLET_AUDIT.tsv"
).open(
    "w",
    newline=""
) as f:

    fields = list(
        all_rows[0].keys()
    )

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(all_rows)


with (
    OUT /
    "TARGET_GENE_PRESENCE_10X.tsv"
).open(
    "w",
    newline=""
) as f:

    fields = [
        "accession",
        "sample_prefix",
        "gene",
        "present",
        "matrix_row_1based",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(gene_rows)


for accession in datasets:

    rr = [
        x for x in all_rows
        if x[
            "accession"
        ] == accession
    ]

    print(
        accession,
        "TOTAL=",
        len(rr),
        "PASS=",
        sum(
            x[
                "dimension_status"
            ] == "PASS"
            for x in rr
        ),
    )

PY

}


# =============================================================================
# 03. GSE189432 ANNOTATION AUDIT
# =============================================================================

audit_gse189432_annotations() {

"$ENV/bin/python" - \
  "$DATA/GSE189432" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import sys


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

p = (
    ROOT /
    "GSE189432_annotations.csv.gz"
)


if not p.exists():
    raise RuntimeError(
        "GSE189432 annotation missing"
    )


with gzip.open(
    p,
    "rt",
    errors="replace"
) as f:

    reader = csv.DictReader(f)

    fields = (
        reader.fieldnames
        or []
    )

    rows = []

    for r in reader:
        rows.append(r)


print(
    "ANNOTATION_ROWS=",
    len(rows)
)

print(
    "ANNOTATION_COLUMNS=",
    fields
)


def candidate_cols(keys):

    return [
        x for x in fields

        if any(
            k in x.lower()
            for k in keys
        )
    ]


barcode_cols = candidate_cols([
    "barcode",
    "cell",
])

celltype_cols = candidate_cols([
    "celltype",
    "cell_type",
    "cell type",
    "cluster",
    "annotation",
    "identity",
    "ident",
])

sample_cols = candidate_cols([
    "sample",
    "orig.ident",
    "gsm",
])

condition_cols = candidate_cols([
    "condition",
    "stroke",
    "time",
    "tissue",
    "group",
])


summary = [{
    "annotation_path":
        str(p),

    "n_rows":
        len(rows),

    "n_columns":
        len(fields),

    "columns":
        ";".join(fields),

    "barcode_candidates":
        ";".join(
            barcode_cols
        ),

    "celltype_candidates":
        ";".join(
            celltype_cols
        ),

    "sample_candidates":
        ";".join(
            sample_cols
        ),

    "condition_candidates":
        ";".join(
            condition_cols
        ),

    "status":
        (
            "PASS"
            if (
                rows
                and celltype_cols
            )
            else "REVIEW"
        ),
}]


with (
    OUT /
    "GSE189432_ANNOTATION_AUDIT.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            summary[0].keys()
        )
    )

    w.writeheader()
    w.writerows(summary)


# Value previews, not assumptions.
preview = []

for col in (
    celltype_cols
    + sample_cols
    + condition_cols
):

    values = []

    seen = set()

    for r in rows:

        v = str(
            r.get(
                col,
                ""
            )
        ).strip()

        if not v:
            continue

        if v in seen:
            continue

        seen.add(v)
        values.append(v)

        if len(values) >= 50:
            break

    preview.append({
        "column":
            col,

        "example_values":
            ";".join(values),
    })


with (
    OUT /
    "GSE189432_METADATA_VALUES.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "column",
            "example_values",
        ]
    )

    w.writeheader()
    w.writerows(preview)


for x in preview:
    print(
        x["column"],
        "=>",
        x["example_values"]
    )

PY

}


# =============================================================================
# 04. GSE225948 COUNTS/METADATA PAIR AUDIT
# =============================================================================

audit_gse225948() {

"$ENV/bin/python" - \
  "$DATA/GSE225948" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import sys


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

files = [
    p
    for p in ROOT.rglob("*.csv.gz")
]


groups = {}


for p in files:

    name = p.name

    if "_counts.csv.gz" in name:

        key = name.replace(
            "_counts.csv.gz",
            ""
        )

        groups.setdefault(
            key,
            {}
        )["counts"] = p


    elif "_metadata.csv.gz" in name:

        key = name.replace(
            "_metadata.csv.gz",
            ""
        )

        groups.setdefault(
            key,
            {}
        )["metadata"] = p


rows = []


def read_header(path):

    with gzip.open(
        path,
        "rt",
        errors="replace"
    ) as f:

        reader = csv.reader(f)

        try:
            header = next(reader)
        except StopIteration:
            return [], []

        first_rows = []

        for _ in range(3):

            try:
                first_rows.append(
                    next(reader)
                )
            except StopIteration:
                break

    return header, first_rows


for key, g in sorted(
    groups.items()
):

    c = g.get("counts")
    m = g.get("metadata")

    complete = bool(
        c and m
    )


    row = {
        "sample":
            key,

        "counts_path":
            str(c or ""),

        "metadata_path":
            str(m or ""),

        "pair_complete":
            complete,

        "counts_header_n":
            "",

        "counts_orientation_guess":
            "",

        "metadata_n_rows":
            "",

        "metadata_columns":
            "",

        "celltype_candidates":
            "",

        "condition_candidates":
            "",

        "status":
            "",
    }


    if not complete:

        row["status"] = (
            "INCOMPLETE_PAIR"
        )

        rows.append(row)
        continue


    ch, cr = read_header(c)

    row[
        "counts_header_n"
    ] = len(ch)


    # Typical submitted matrix:
    # first column gene symbol and remaining columns cells.
    first_data_first = (
        cr[0][0]
        if cr
        and cr[0]
        else ""
    )


    target_symbols = {
        "Fgf5",
        "Aldh2",
        "Sh3pxd2a",
        "Col4a2",
        "Col4a1",
        "Calhm2",
        "Neurl1",
        "Ina",
    }


    if any(
        x in target_symbols
        for x in [
            first_data_first,
            *[
                r[0]
                for r in cr
                if r
            ],
        ]
    ):

        orientation = (
            "GENES_BY_CELLS_LIKELY"
        )

    elif any(
        x in target_symbols
        for x in ch
    ):

        orientation = (
            "CELLS_BY_GENES_LIKELY"
        )

    elif len(ch) > 1000:

        orientation = (
            "GENES_BY_CELLS_LIKELY"
        )

    else:

        orientation = (
            "REVIEW"
        )


    row[
        "counts_orientation_guess"
    ] = orientation


    with gzip.open(
        m,
        "rt",
        errors="replace"
    ) as f:

        reader = csv.DictReader(f)

        mf = (
            reader.fieldnames
            or []
        )

        meta = list(reader)


    row[
        "metadata_n_rows"
    ] = len(meta)

    row[
        "metadata_columns"
    ] = ";".join(mf)


    cell_cols = [
        x for x in mf

        if any(
            z in x.lower()

            for z in [
                "celltype",
                "cell_type",
                "cell type",
                "cluster",
                "annotation",
                "ident",
            ]
        )
    ]


    cond_cols = [
        x for x in mf

        if any(
            z in x.lower()

            for z in [
                "stroke",
                "sham",
                "condition",
                "day",
                "age",
                "sex",
                "group",
                "sample",
            ]
        )
    ]


    row[
        "celltype_candidates"
    ] = ";".join(
        cell_cols
    )

    row[
        "condition_candidates"
    ] = ";".join(
        cond_cols
    )


    row["status"] = (
        "PASS"
        if (
            orientation != "REVIEW"
            and len(meta) > 0
        )
        else "REVIEW"
    )


    rows.append(row)


with (
    OUT /
    "GSE225948_PAIR_AUDIT.tsv"
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


print(
    "GSE225948_SAMPLE_PAIRS=",
    len(rows)
)

print(
    "PASS=",
    sum(
        x["status"] == "PASS"
        for x in rows
    )
)


for x in rows:

    print(
        x["sample"],
        x["status"],
        x[
            "counts_orientation_guess"
        ],
        "metadata_cells=",
        x[
            "metadata_n_rows"
        ],
        "celltype=",
        x[
            "celltype_candidates"
        ],
        "condition=",
        x[
            "condition_candidates"
        ],
    )

PY

}


# =============================================================================
# 05. GSE234052 AGGREGATED MATRIX + BARCODE SUFFIX AUDIT
# =============================================================================

audit_gse234052() {

"$ENV/bin/python" - \
  "$DATA/GSE234052" \
  "$OUT" <<'PY'

from pathlib import Path
from collections import Counter
import csv
import gzip
import re
import sys


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

barcodes = (
    ROOT /
    "GSE234052_barcodes.tsv.gz"
)

features = (
    ROOT /
    "GSE234052_features.tsv.gz"
)

matrix = (
    ROOT /
    "GSE234052_matrix.mtx.gz"
)


for p in [
    barcodes,
    features,
    matrix,
]:

    if not p.exists():
        raise RuntimeError(
            f"missing: {p}"
        )


# Matrix dimensions.
with gzip.open(
    matrix,
    "rt",
    errors="replace"
) as f:

    dims = None

    for line in f:

        line = line.strip()

        if not line:
            continue

        if line.startswith("%"):
            continue

        x = line.split()

        dims = tuple(
            map(
                int,
                x[:3]
            )
        )

        break


if dims is None:
    raise RuntimeError(
        "cannot read matrix dimensions"
    )


# Feature rows and target genes.
targets = {
    "Fgf5",
    "Aldh2",
    "Sh3pxd2a",
    "Col4a2",
    "Col4a1",
    "Calhm2",
    "Neurl1",
    "Ina",
}


target_positions = {}
feature_n = 0

with gzip.open(
    features,
    "rt",
    errors="replace"
) as f:

    for i, line in enumerate(
        f,
        start=1
    ):

        feature_n += 1

        x = line.rstrip(
            "\n"
        ).split("\t")

        for token in x:

            if token in targets:

                target_positions[
                    token
                ] = i


suffix_counts = Counter()
barcode_n = 0


with gzip.open(
    barcodes,
    "rt",
    errors="replace"
) as f:

    for line in f:

        b = line.strip()

        if not b:
            continue

        barcode_n += 1

        m = re.search(
            r"-(\d+)$",
            b
        )

        suffix = (
            m.group(1)
            if m
            else "UNRESOLVED"
        )

        suffix_counts[
            suffix
        ] += 1


mapping = {

    "1":
        "Contra_1h",

    "2":
        "IPSI_1h",

    "3":
        "Contra_12h",

    "4":
        "IPSI_12h",

    "5":
        "Contra_24h",

    "6":
        "IPSI_24h",
}


rows = []

for suffix in sorted(
    suffix_counts,
    key=lambda x:
        (
            x == "UNRESOLVED",
            x
        )
):

    rows.append({

        "suffix":
            suffix,

        "condition":
            mapping.get(
                suffix,
                "UNRESOLVED"
            ),

        "n_cells":
            suffix_counts[
                suffix
            ],
    })


with (
    OUT /
    "GSE234052_BARCODE_SUFFIX_AUDIT.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "suffix",
            "condition",
            "n_cells",
        ]
    )

    w.writeheader()
    w.writerows(rows)


gene_rows = []

for gene in sorted(
    targets
):

    gene_rows.append({

        "gene":
            gene,

        "present":
            gene
            in target_positions,

        "matrix_row_1based":
            target_positions.get(
                gene,
                ""
            ),
    })


with (
    OUT /
    "GSE234052_TARGET_GENE_PRESENCE.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "gene",
            "present",
            "matrix_row_1based",
        ]
    )

    w.writeheader()
    w.writerows(
        gene_rows
    )


print(
    "MATRIX_DIMENSIONS=",
    dims
)

print(
    "FEATURE_ROWS=",
    feature_n
)

print(
    "BARCODE_ROWS=",
    barcode_n
)

print(
    "DIMENSION_STATUS=",
    (
        "PASS"
        if (
            dims[0] == feature_n
            and dims[1] == barcode_n
        )
        else "MISMATCH"
    )
)


for x in rows:

    print(
        x["suffix"],
        x["condition"],
        x["n_cells"]
    )


print(
    "TARGETS=",
    target_positions
)

PY

}


# =============================================================================
# 06. GSE256490 SAMPLE MANIFEST
# =============================================================================

audit_gse256490_samples() {

"$ENV/bin/python" - \
  "$OUT/TENX_TRIPLET_AUDIT.tsv" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import re
import sys


INPUT = Path(sys.argv[1])
OUT = Path(sys.argv[2])


with INPUT.open() as f:

    rows = [
        x for x in csv.DictReader(
            f,
            delimiter="\t"
        )

        if x[
            "accession"
        ] == "GSE256490"
    ]


out = []


for x in rows:

    name = x[
        "sample_prefix"
    ]

    low = name.lower()


    if "fetal" in low:

        broad_group = (
            "FETAL"
        )

    elif any(
        z in low
        for z in [
            "tle",
            "tl_",
            "temporal",
        ]
    ):

        broad_group = (
            "ADULT_TEMPORAL_LOBE"
        )

    elif "avm" in low:

        broad_group = (
            "AVM"
        )

    elif "gbm" in low:

        broad_group = (
            "GBM"
        )

    elif "lgg" in low:

        broad_group = (
            "LGG"
        )

    elif "met" in low:

        broad_group = (
            "METASTASIS"
        )

    elif "men" in low:

        broad_group = (
            "MENINGIOMA"
        )

    else:

        broad_group = (
            "OTHER"
        )


    if "unsorted" in low:

        fraction = (
            "UNSORTED_ENDOTHELIAL_PERIVASCULAR"
        )

    elif (
        "sorted" in low
        or "_ec" in low
    ):

        fraction = (
            "SORTED_ENDOTHELIAL_ENRICHED"
        )

    else:

        fraction = (
            "UNKNOWN"
        )


    out.append({

        "sample_prefix":
            name,

        "broad_group_hint":
            broad_group,

        "fraction_hint":
            fraction,

        "n_cells":
            x[
                "matrix_cells"
            ],

        "dimension_status":
            x[
                "dimension_status"
            ],

        "use_for_target_analysis":
            (
                "CANDIDATE_REFERENCE_RAW"
                if (
                    broad_group
                    == "ADULT_TEMPORAL_LOBE"
                    and fraction
                    == "UNSORTED_ENDOTHELIAL_PERIVASCULAR"
                )
                else "SECONDARY"
            ),

        "annotation_status":
            "NO_AUTHOR_CELLTYPE_LABEL_IN_RAW_TRIPLET",
    })


with (
    OUT /
    "GSE256490_SAMPLE_AUDIT.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            out[0].keys()
        )
    )

    w.writeheader()
    w.writerows(out)


for x in out:

    if (
        x[
            "use_for_target_analysis"
        ]
        == "CANDIDATE_REFERENCE_RAW"
    ):

        print(
            "REFERENCE_CANDIDATE",
            x[
                "sample_prefix"
            ],
            "cells=",
            x[
                "n_cells"
            ],
        )

PY

}


# =============================================================================
# 07. GSE256493 SMALL MANIFEST FILES
# =============================================================================

prepare_author_annotation_manifest() {

  D="$DATA/GSE256493_manifest"

  mkdir -p "$D"

  BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256493/suppl"


  download_small() {

    local URL="$1"
    local DEST="$2"

    if [[ -s "$DEST" ]]; then
      echo "EXISTS $DEST"
      return 0
    fi

    curl \
      -L \
      --fail \
      --retry 3 \
      --connect-timeout 20 \
      -o "$DEST.part" \
      "$URL"

    RC=$?

    if [[ $RC -ne 0 ]]; then
      rm -f "$DEST.part"
      return 1
    fi

    mv \
      "$DEST.part" \
      "$DEST"

    return 0
  }


  download_small \
    "$BASE/GSE256493_readme_rds_descriptions.txt" \
    "$D/GSE256493_readme_rds_descriptions.txt"

  RC1=$?


  download_small \
    "$BASE/GSE256493_scRNAseq_seurat_objects.xlsx" \
    "$D/GSE256493_scRNAseq_seurat_objects.xlsx"

  RC2=$?


  echo
  echo "===== README ====="

  cat \
    "$D/GSE256493_readme_rds_descriptions.txt" \
    2>/dev/null \
    | head -200


  echo
  echo "===== TARGET AUTHOR OBJECT ====="

  cat > \
    "$OUT/GSE256493_ANNOTATED_OBJECT_PLAN.tsv" <<'EOF'
purpose	filename	approx_size	status
PRIMARY_HUMAN_VASCULAR_REFERENCE	GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds.gz	3.3GB	NOT_DOWNLOADED
SECONDARY_ENDOTHELIAL_REFERENCE	GSE256493_Adult_control_brain_temporal_lobe_sorted_endothelial_cells_seurat_object.rds.gz	3.3GB	NOT_DOWNLOADED
EOF

  cat \
    "$OUT/GSE256493_ANNOTATED_OBJECT_PLAN.tsv"


  if [[ $RC1 -ne 0 || $RC2 -ne 0 ]]; then

    echo "MANIFEST_DOWNLOAD_PARTIAL"

    return 1

  fi

  echo "AUTHOR_MANIFEST=PASS"
}


# =============================================================================
# 08. BUILD READINESS
# =============================================================================

build_readiness() {

"$ENV/bin/python" - \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


def read(name):

    p = OUT / name

    if not p.exists():
        return []

    with p.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


tenx = read(
    "TENX_TRIPLET_AUDIT.tsv"
)

g189 = read(
    "GSE189432_ANNOTATION_AUDIT.tsv"
)

g225 = read(
    "GSE225948_PAIR_AUDIT.tsv"
)

g234 = read(
    "GSE234052_BARCODE_SUFFIX_AUDIT.tsv"
)


def triplet_status(acc):

    rows = [
        x for x in tenx

        if x[
            "accession"
        ] == acc
    ]

    if not rows:
        return "MISSING"

    if all(
        x[
            "dimension_status"
        ] == "PASS"

        for x in rows
    ):
        return (
            f"PASS_{len(rows)}_TRIPLETS"
        )

    return "REVIEW"


status189 = (
    "PASS"
    if (
        g189
        and g189[0][
            "status"
        ] == "PASS"
    )
    else "REVIEW"
)


status225 = (
    "PASS"
    if (
        g225
        and all(
            x["status"] == "PASS"
            for x in g225
        )
    )
    else "REVIEW"
)


suffixes = {
    x["suffix"]
    for x in g234
}


status234 = (
    "PASS"
    if suffixes
    >= {
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
    }
    else "REVIEW"
)


rows = [

    (
        "PHASE9F_B_FORMAT_AUDIT",
        "COMPLETE"
    ),

    (
        "GSE189432_10X",
        triplet_status(
            "GSE189432"
        )
    ),

    (
        "GSE189432_ANNOTATION",
        status189
    ),

    (
        "GSE225948_COUNTS_METADATA",
        status225
    ),

    (
        "GSE234052_AGGREGATED_MATRIX",
        status234
    ),

    (
        "GSE256490_RAW_10X",
        triplet_status(
            "GSE256490"
        )
    ),

    (
        "GSE256490_CELLTYPE_LABELS",
        "AUTHOR_ANNOTATED_SEURAT_OBJECT_REQUIRED"
    ),

    (
        "HUMAN_VASCULAR_REFERENCE_NEXT",
        "ACQUIRE_GSE256493_ADULT_CONTROL_UNSORTED_RDS"
    ),

    (
        "SH3PXD2A_MOUSE_NEXT",
        (
            "READY_FOR_TARGET_EXPRESSION"
            if (
                status189 == "PASS"
                or status225 == "PASS"
            )
            else "PENDING_FORMAT_REVIEW"
        )
    ),

    (
        "COL4A2_COL4A1_MOUSE_NEXT",
        (
            "READY_FOR_TARGET_EXPRESSION"
            if status234 == "PASS"
            else "PENDING_FORMAT_REVIEW"
        )
    ),

    (
        "NEXT_STAGE",
        "PHASE9F_C_TARGET_EXPRESSION_AND_HUMAN_ANNOTATED_REFERENCE"
    ),
]


with (
    OUT /
    "PHASE9F_B_READINESS.tsv"
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
  "02_tenx_triplet_audit" \
  audit_10x_triplets

run_step \
  "03_gse189432_annotation_audit" \
  audit_gse189432_annotations

run_step \
  "04_gse225948_pair_audit" \
  audit_gse225948

run_step \
  "05_gse234052_suffix_audit" \
  audit_gse234052

run_step \
  "06_gse256490_sample_audit" \
  audit_gse256490_samples

run_step \
  "07_author_annotation_manifest" \
  prepare_author_annotation_manifest

run_step \
  "08_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9F-B COMPLETE"
echo "================================================================================"


echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== 10X SUMMARY ====="

python3 - \
  "$OUT/TENX_TRIPLET_AUDIT.tsv" <<'PY'

import csv
import collections
import sys

p = sys.argv[1]

d = collections.defaultdict(
    lambda: [
        0,
        0,
        0,
    ]
)

with open(p) as f:

    for x in csv.DictReader(
        f,
        delimiter="\t"
    ):

        a = x["accession"]

        d[a][0] += 1

        d[a][1] += (
            x["dimension_status"]
            == "PASS"
        )

        try:
            d[a][2] += int(
                x["matrix_cells"]
            )
        except Exception:
            pass


for a, v in sorted(
    d.items()
):

    print(
        a,
        "triplets=",
        v[0],
        "pass=",
        v[1],
        "cells_total=",
        v[2],
    )

PY


echo
echo "===== GSE189432 ANNOTATION ====="

cat \
  "$OUT/GSE189432_ANNOTATION_AUDIT.tsv" \
  2>/dev/null || true


echo
echo "===== GSE225948 ====="

column -t -s $'\t' \
  "$OUT/GSE225948_PAIR_AUDIT.tsv" \
  2>/dev/null \
  | head -40 \
  || true


echo
echo "===== GSE234052 CONDITIONS ====="

cat \
  "$OUT/GSE234052_BARCODE_SUFFIX_AUDIT.tsv" \
  2>/dev/null || true


echo
echo "===== HUMAN ANNOTATED OBJECT PLAN ====="

cat \
  "$OUT/GSE256493_ANNOTATED_OBJECT_PLAN.tsv" \
  2>/dev/null || true


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9F_B_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE9F_B_READINESS.tsv"


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
