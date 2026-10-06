#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_phase9a/$RUN_ID"

OUTROOT="$ROOT/results/is/stage5_functional/phase9"
QCROOT="$ROOT/results/is/stage4_cross_eas/phase8_qc"

STATUS="$LOGROOT/STATUS.tsv"
READINESS="$ROOT/results/is/stage0_registry/IS_PHASE9A_STATUS.tsv"

mkdir -p \
    "$LOGROOT" \
    "$OUTROOT" \
    "$QCROOT/kriging_fixed"

cd "$REPO" || exit 1

printf "step\tstatus\tmessage\tlog\n" > "$STATUS"


record() {
    printf "%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "$4" >> "$STATUS"
}


run_step() {

    STEP="$1"
    shift

    LOG="$LOGROOT/${STEP}.log"

    echo
    echo "===================================================="
    echo "STEP=$STEP"
    echo "START=$(date)"
    echo "LOG=$LOG"
    echo "===================================================="

    (
        "$@"
    ) > >(tee "$LOG") 2>&1

    RC=$?

    if [ "$RC" -eq 0 ]; then

        echo "[PASS] $STEP"

        record \
            "$STEP" \
            "PASS" \
            "completed" \
            "$LOG"

    else

        echo "[FAIL] $STEP rc=$RC"

        record \
            "$STEP" \
            "FAIL" \
            "rc=$RC" \
            "$LOG"

    fi

    return 0
}


# ============================================================
# STEP 01
# BUILD BBJ CREDIBLE-SET MASTER
# ============================================================

build_credible_set_master() {

python3 - <<'PY'
import csv
from pathlib import Path

ROOT = Path("/srv/is-analysis")

BBJ = (
    ROOT /
    "results/is/stage3_finemap/japan/bbj"
)

OUT = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "BBJ_CREDIBLE_SET_MASTER.tsv"
)

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004",
]

rows = []
fail = 0

for locus in loci:

    pipfile = (
        BBJ /
        "susie_v3" /
        locus /
        "PIP.tsv"
    )

    csfile = (
        BBJ /
        "susie_v3" /
        locus /
        "CREDIBLE_SETS.tsv"
    )

    if not pipfile.exists():
        print("MISSING", pipfile)
        fail += 1
        continue

    if not csfile.exists():
        print("MISSING", csfile)
        fail += 1
        continue

    pips = {}

    with pipfile.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):
            pips[r["variant_id"]] = r

    with csfile.open() as f:

        csrows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    for r in csrows:

        vid = r["variant_id"]

        parts = vid.split(":")

        if len(parts) != 4:
            print(
                "BAD_VARIANT_ID",
                locus,
                vid
            )
            continue

        chrom, pos, ref, alt = parts

        p = pips.get(
            vid,
            {}
        )

        rows.append({
            "locus":
                locus,

            "credible_set":
                r.get(
                    "credible_set",
                    ""
                ),

            "variant_id":
                vid,

            "chr":
                chrom,

            "pos_grch37":
                pos,

            "ref":
                ref,

            "alt":
                alt,

            "bbj_pip":
                r.get(
                    "pip",
                    p.get(
                        "pip",
                        ""
                    )
                ),

            "bbj_beta":
                p.get(
                    "beta",
                    ""
                ),

            "bbj_se":
                p.get(
                    "se",
                    ""
                ),

            "bbj_p":
                p.get(
                    "p",
                    ""
                )
        })


if fail:
    raise SystemExit(
        f"INPUT_FAILURES={fail}"
    )

if not rows:
    raise SystemExit(
        "NO_CREDIBLE_VARIANTS"
    )


fields = list(
    rows[0].keys()
)


with OUT.open(
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(
        rows
    )


print(
    "N_CREDIBLE_VARIANTS=",
    len(rows)
)

print()


expected = {
    "BBJ_IS_L001": 3,
    "BBJ_IS_L002": 14,
    "BBJ_IS_L003": 5,
    "BBJ_IS_L004": 2,
}


for locus in loci:

    x = [
        r
        for r in rows
        if r["locus"] == locus
    ]

    print(
        locus,
        "N=",
        len(x),
        "EXPECTED=",
        expected[locus],
        "PASS=",
        len(x) == expected[locus]
    )

    for r in sorted(
        x,
        key=lambda z:
            -float(
                z["bbj_pip"]
            )
    ):

        print(
            " ",
            r["variant_id"],
            "PIP=",
            r["bbj_pip"]
        )


bad = [
    locus
    for locus, n in expected.items()
    if sum(
        r["locus"] == locus
        for r in rows
    ) != n
]

print(
    "\nOUT=",
    OUT
)

raise SystemExit(
    0 if not bad else 2
)
PY

    return $?
}


# ============================================================
# STEP 02
# GRCh37 LOCUS -> GENE UNIVERSE
# ============================================================

annotate_genes_grch37() {

    INP="$OUTROOT/BBJ_CREDIBLE_SET_MASTER.tsv"

    if [ ! -s "$INP" ]; then
        echo "CREDIBLE_SET_MASTER_MISSING"
        return 21
    fi


python3 - <<'PY'
import csv
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path("/srv/is-analysis")

INP = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "BBJ_CREDIBLE_SET_MASTER.tsv"
)

OUT = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "GRCH37_LOCUS_GENE_UNIVERSE.tsv"
)

PROTEIN_OUT = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "GRCH37_PROTEIN_CODING_GENES.tsv"
)


locus_bounds = {}


with INP.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        locus = r["locus"]
        chrom = r["chr"]
        pos = int(
            r["pos_grch37"]
        )

        if locus not in locus_bounds:

            locus_bounds[locus] = {
                "chr":
                    chrom,

                "min":
                    pos,

                "max":
                    pos
            }

        else:

            locus_bounds[locus]["min"] = min(
                locus_bounds[locus]["min"],
                pos
            )

            locus_bounds[locus]["max"] = max(
                locus_bounds[locus]["max"],
                pos
            )


rows = []


for locus in sorted(
    locus_bounds
):

    x = locus_bounds[locus]

    chrom = x["chr"]

    start = max(
        1,
        x["min"] - 500000
    )

    end = (
        x["max"] + 500000
    )

    region = (
        f"{chrom}:{start}-{end}"
    )

    url = (
        "https://grch37.rest.ensembl.org/"
        f"overlap/region/human/{region}"
        "?feature=gene"
    )

    print()
    print(
        "QUERY",
        locus,
        region
    )

    req = urllib.request.Request(
        url,
        headers={
            "Content-Type":
                "application/json"
        }
    )

    data = None

    for attempt in range(
        1,
        4
    ):

        try:

            with urllib.request.urlopen(
                req,
                timeout=90
            ) as response:

                data = json.loads(
                    response.read()
                )

            break

        except Exception as e:

            print(
                "ATTEMPT",
                attempt,
                "ERROR",
                repr(e)
            )

            time.sleep(
                2 * attempt
            )


    if data is None:

        print(
            "FAILED",
            locus
        )

        continue


    print(
        "GENES=",
        len(data)
    )


    for g in data:

        rows.append({
            "locus":
                locus,

            "query_region":
                region,

            "ensembl_gene_id":
                g.get(
                    "id",
                    ""
                ),

            "gene_symbol":
                g.get(
                    "external_name",
                    ""
                ),

            "biotype":
                g.get(
                    "biotype",
                    ""
                ),

            "chr":
                g.get(
                    "seq_region_name",
                    ""
                ),

            "start_grch37":
                g.get(
                    "start",
                    ""
                ),

            "end_grch37":
                g.get(
                    "end",
                    ""
                ),

            "strand":
                g.get(
                    "strand",
                    ""
                )
        })


    time.sleep(
        0.25
    )


if not rows:

    raise SystemExit(
        "NO_GENES_RETURNED"
    )


fields = list(
    rows[0].keys()
)


with OUT.open(
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
        rows
    )


protein_rows = [
    r
    for r in rows
    if r["biotype"] == "protein_coding"
]


with PROTEIN_OUT.open(
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
        protein_rows
    )


print()
print(
    "===== LOCUS GENE COUNTS ====="
)


for locus in sorted(
    locus_bounds
):

    all_genes = [
        r
        for r in rows
        if r["locus"] == locus
    ]

    protein = [
        r
        for r in all_genes
        if r["biotype"]
        == "protein_coding"
    ]

    print()
    print(
        locus,
        "ALL=",
        len(all_genes),
        "PROTEIN_CODING=",
        len(protein)
    )

    for r in protein:

        print(
            " ",
            r["gene_symbol"],
            r["ensembl_gene_id"],
            f'{r["chr"]}:{r["start_grch37"]}-{r["end_grch37"]}'
        )


print()
print(
    "OUT=",
    OUT
)

print(
    "PROTEIN_OUT=",
    PROTEIN_OUT
)
PY

    return $?
}


# ============================================================
# STEP 03
# LOCAL FUNCTIONAL RESOURCE INVENTORY
# ============================================================

inventory_resources() {

python3 - <<'PY'
from pathlib import Path
import csv
import os

ROOTS = [
    Path("/srv/is-analysis/data"),
    Path("/srv/is-analysis/results"),
]

OUT = Path(
    "/srv/is-analysis/results/is/"
    "stage5_functional/phase9/"
    "LOCAL_FUNCTIONAL_RESOURCE_INVENTORY.tsv"
)

SUMMARY = Path(
    "/srv/is-analysis/results/is/"
    "stage5_functional/phase9/"
    "LOCAL_FUNCTIONAL_RESOURCE_SUMMARY.tsv"
)


keywords = [
    "eqtl",
    "gtex",
    "brain",
    "blood",
    "artery",
    "vascular",
    "endothelial",
    "pericyte",
    "smooth_muscle",
    "singlecell",
    "single_cell",
    "scrna",
    "snrna",
    "atac",
    "chromatin",
    "pqtl",
    "ukbppp",
    "decode",
    "fenland",
    "coloc",
    "stroke",
    "ischemic",
]


rows = []
seen = set()


for root in ROOTS:

    if not root.exists():
        continue

    for dirpath, dirnames, filenames in os.walk(
        root
    ):

        # Skip hidden/system-like trees.
        dirnames[:] = [
            d
            for d in dirnames
            if not d.startswith(".")
        ]

        for filename in filenames:

            p = Path(
                dirpath
            ) / filename

            full = str(
                p
            )

            low = full.lower()

            hits = [
                k
                for k in keywords
                if k in low
            ]

            if not hits:
                continue

            if full in seen:
                continue

            seen.add(
                full
            )

            try:
                size = p.stat().st_size
            except Exception:
                size = -1


            rows.append({
                "path":
                    full,

                "filename":
                    filename,

                "size_bytes":
                    size,

                "size_mb":
                    (
                        round(
                            size /
                            1024 /
                            1024,
                            3
                        )
                        if size >= 0
                        else ""
                    ),

                "keywords":
                    ",".join(
                        hits
                    )
            })


rows.sort(
    key=lambda x:
        (
            x["keywords"],
            x["path"]
        )
)


fields = [
    "path",
    "filename",
    "size_bytes",
    "size_mb",
    "keywords"
]


with OUT.open(
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
        rows
    )


counts = {}

for r in rows:

    for k in r[
        "keywords"
    ].split(","):

        counts[k] = (
            counts.get(
                k,
                0
            )
            + 1
        )


summary_rows = [
    {
        "keyword":
            k,

        "n_files":
            v
    }
    for k, v in sorted(
        counts.items()
    )
]


with SUMMARY.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "keyword",
            "n_files"
        ]
    )

    w.writeheader()
    w.writerows(
        summary_rows
    )


print(
    "RESOURCE_FILES=",
    len(rows)
)


print()
print(
    "===== RESOURCE SUMMARY ====="
)


for r in summary_rows:

    print(
        r["keyword"],
        r["n_files"]
    )


print()
print(
    "===== LARGE / LIKELY USEFUL FILES ====="
)


large = sorted(
    rows,
    key=lambda x:
        x["size_bytes"],
    reverse=True
)


for r in large[:100]:

    print(
        r["size_mb"],
        "MB",
        r["keywords"],
        r["path"]
    )


print()
print(
    "OUT=",
    OUT
)

print(
    "SUMMARY=",
    SUMMARY
)
PY

    return $?
}


# ============================================================
# STEP 04
# FIXED BBJ KRIGING / ALLELE-SWITCH AUDIT
# ============================================================

kriging_audit() {

Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

suppressPackageStartupMessages(
    library(susieR)
)


ROOT <- paste0(
    "/srv/is-analysis/results/is/",
    "stage3_finemap/japan/bbj"
)

OUT <- paste0(
    "/srv/is-analysis/results/is/",
    "stage4_cross_eas/phase8_qc/",
    "kriging_fixed"
)

dir.create(
    OUT,
    recursive=TRUE,
    showWarnings=FALSE
)


loci <- c(
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
)


N_EFF <- round(
    4 /
    (
        1/22664 +
        1/152022
    )
)


read_ld <- function(
    path,
    n
) {

    expected <- n*n*8

    actual <- file.info(
        path
    )$size


    if (
        is.na(actual)
        || actual != expected
    ) {

        stop(
            paste0(
                "LD size mismatch ",
                actual,
                " != ",
                expected
            )
        )
    }


    con <- file(
        path,
        "rb"
    )

    on.exit(
        close(con)
    )


    x <- readBin(
        con,
        what="numeric",
        n=n*n,
        size=8,
        endian="little"
    )


    matrix(
        x,
        nrow=n,
        ncol=n,
        byrow=TRUE
    )
}


summary_list <- list()


for (locus in loci) {

    cat(
        "\n========================================\n",
        locus,
        "\n========================================\n",
        sep=""
    )


    result <- tryCatch({

        ss <- fread(
            file.path(
                ROOT,
                "susie_inputs_v3",
                paste0(
                    locus,
                    ".tsv"
                )
            )
        )


        vars <- readLines(
            file.path(
                ROOT,
                "ld_v3",
                paste0(
                    locus,
                    ".unphased.vcor1.bin.vars"
                )
            )
        )


        if (
            nrow(ss) != length(vars)
        ) {

            stop(
                "summary/LD dimension mismatch"
            )
        }


        if (
            !all(
                ss$variant_id == vars
            )
        ) {

            stop(
                "summary/LD variant order mismatch"
            )
        }


        R <- read_ld(
            file.path(
                ROOT,
                "ld_v3",
                paste0(
                    locus,
                    ".unphased.vcor1.bin"
                )
            ),
            length(vars)
        )


        R <- (
            R +
            t(R)
        ) / 2

        diag(R) <- 1


        z <- (
            as.numeric(
                ss$beta
            ) /
            as.numeric(
                ss$se
            )
        )


        s <- estimate_s_rss(
            z=z,
            R=R,
            n=N_EFF,
            method="null-mle"
        )


        kr <- kriging_rss(
            z=z,
            R=R,
            n=N_EFF,
            s=s
        )


        cat(
            "RETURN_NAMES=",
            paste(
                names(kr),
                collapse=","
            ),
            "\n",
            sep=""
        )


        n_switch <- NA_integer_
        table_name <- ""


        if (
            !is.null(
                kr$conditional_dist
            )
        ) {

            tab <- as.data.table(
                kr$conditional_dist
            )

            table_name <- (
                "conditional_dist"
            )

        } else if (
            !is.null(
                kr$table
            )
        ) {

            tab <- as.data.table(
                kr$table
            )

            table_name <- (
                "table"
            )

        } else {

            tab <- data.table()

        }


        if (
            nrow(tab) > 0
        ) {

            if (
                nrow(tab)
                != length(vars)
            ) {

                stop(
                    paste0(
                        "kriging rows=",
                        nrow(tab),
                        " variants=",
                        length(vars)
                    )
                )
            }


            tab[
                ,
                variant_id := vars
            ]


            fwrite(
                tab,
                file.path(
                    OUT,
                    paste0(
                        locus,
                        ".",
                        table_name,
                        ".tsv"
                    )
                ),
                sep="\t"
            )


            cat(
                "TABLE=",
                table_name,
                "\n",
                sep=""
            )

            cat(
                "COLUMNS=",
                paste(
                    names(tab),
                    collapse=","
                ),
                "\n",
                sep=""
            )


            lr_col <- grep(
                "^logLR$|log.*LR|log.*lik",
                names(tab),
                ignore.case=TRUE,
                value=TRUE
            )


            z_col <- grep(
                "^z$|zscore|z_score",
                names(tab),
                ignore.case=TRUE,
                value=TRUE
            )


            if (
                length(lr_col) >= 1
            ) {

                lr <- suppressWarnings(
                    as.numeric(
                        tab[[lr_col[1]]]
                    )
                )


                if (
                    length(z_col) >= 1
                ) {

                    z_for_switch <- suppressWarnings(
                        as.numeric(
                            tab[[z_col[1]]]
                        )
                    )

                } else {

                    z_for_switch <- z
                }


                suspect <- tab[
                    is.finite(lr)
                    &
                    lr > 2
                    &
                    abs(
                        z_for_switch
                    ) > 2
                ]


                n_switch <- nrow(
                    suspect
                )


                fwrite(
                    suspect,
                    file.path(
                        OUT,
                        paste0(
                            locus,
                            ".possible_switch.tsv"
                        )
                    ),
                    sep="\t"
                )
            }

        }


        cat(
            "s=",
            s,
            "\n",
            sep=""
        )

        cat(
            "possible_switch_n=",
            n_switch,
            "\n",
            sep=""
        )


        data.table(
            locus=
                locus,

            status=
                "PASS",

            n_variants=
                length(vars),

            n_eff=
                N_EFF,

            s=
                s,

            kriging_table=
                table_name,

            possible_switch_n=
                n_switch,

            error=
                ""
        )

    }, error=function(e) {

        cat(
            "ERROR=",
            conditionMessage(e),
            "\n",
            sep=""
        )


        data.table(
            locus=
                locus,

            status=
                "FAIL",

            n_variants=
                NA_integer_,

            n_eff=
                N_EFF,

            s=
                NA_real_,

            kriging_table=
                "",

            possible_switch_n=
                NA_integer_,

            error=
                conditionMessage(e)
        )
    })


    summary_list[
        [length(summary_list) + 1]
    ] <- result
}


summary <- rbindlist(
    summary_list,
    fill=TRUE
)


fwrite(
    summary,
    file.path(
        OUT,
        "BBJ_KRIGING_FIXED_SUMMARY.tsv"
    ),
    sep="\t"
)


cat(
    "\n===== KRIGING SUMMARY =====\n"
)

print(
    summary
)

RS

    return $?
}


# ============================================================
# STEP 05
# BUILD PHASE 9A LOCUS SUMMARY
# ============================================================

build_phase9_summary() {

python3 - <<'PY'
import csv
from pathlib import Path

ROOT = Path(
    "/srv/is-analysis"
)

CS = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "BBJ_CREDIBLE_SET_MASTER.tsv"
)

GENE = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "GRCH37_PROTEIN_CODING_GENES.tsv"
)

KR = (
    ROOT /
    "results/is/stage4_cross_eas/"
    "phase8_qc/kriging_fixed/"
    "BBJ_KRIGING_FIXED_SUMMARY.tsv"
)

OUT = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "PHASE9A_LOCUS_SUMMARY.tsv"
)


loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004",
]


cs = {}

with CS.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        cs.setdefault(
            r["locus"],
            []
        ).append(
            r
        )


genes = {}

with GENE.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        genes.setdefault(
            r["locus"],
            []
        ).append(
            r
        )


kr = {}

if KR.exists():

    with KR.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):

            kr[
                r["locus"]
            ] = r


rows = []


for locus in loci:

    variants = cs.get(
        locus,
        []
    )

    protein = genes.get(
        locus,
        []
    )


    top = max(
        variants,
        key=lambda x:
            float(
                x["bbj_pip"]
            )
    )


    gene_symbols = sorted(
        {
            x["gene_symbol"]
            for x in protein
            if x[
                "gene_symbol"
            ]
        }
    )


    k = kr.get(
        locus,
        {}
    )


    rows.append({
        "locus":
            locus,

        "n_bbj_cs_variants":
            len(
                variants
            ),

        "bbj_top_variant":
            top[
                "variant_id"
            ],

        "bbj_top_pip":
            top[
                "bbj_pip"
            ],

        "n_protein_coding_genes":
            len(
                gene_symbols
            ),

        "protein_coding_genes":
            ";".join(
                gene_symbols
            ),

        "ld_s":
            k.get(
                "s",
                ""
            ),

        "possible_switch_n":
            k.get(
                "possible_switch_n",
                ""
            ),

        "kriging_status":
            k.get(
                "status",
                ""
            )
    })


fields = list(
    rows[0].keys()
)


with OUT.open(
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
        rows
    )


print(
    "\t".join(
        fields
    )
)

for r in rows:

    print(
        "\t".join(
            str(r[x])
            for x in fields
        )
    )


print(
    "\nOUT=",
    OUT
)
PY

    return $?
}


# ============================================================
# STEP 06
# FINAL READINESS
# ============================================================

final_status() {

    CS="$OUTROOT/BBJ_CREDIBLE_SET_MASTER.tsv"
    GENE="$OUTROOT/GRCH37_LOCUS_GENE_UNIVERSE.tsv"
    PCG="$OUTROOT/GRCH37_PROTEIN_CODING_GENES.tsv"
    INV="$OUTROOT/LOCAL_FUNCTIONAL_RESOURCE_INVENTORY.tsv"

    KR="$QCROOT/kriging_fixed/BBJ_KRIGING_FIXED_SUMMARY.tsv"

    SUMMARY="$OUTROOT/PHASE9A_LOCUS_SUMMARY.tsv"


    CS_N=0
    LOCUS_N=0
    KR_PASS=0


    if [ -s "$CS" ]; then

        CS_N="$(
            awk '
                NR>1 {n++}
                END {print n+0}
            ' "$CS"
        )"

        LOCUS_N="$(
            awk -F '\t' '
                NR>1 {
                    x[$1]=1
                }
                END {
                    print length(x)
                }
            ' "$CS"
        )"
    fi


    if [ -s "$KR" ]; then

        KR_PASS="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {
                    n++
                }
                END {
                    print n+0
                }
            ' "$KR"
        )"
    fi


    printf \
        "component\tstatus\tblocker\n" \
        > "$READINESS"


    if [ "$CS_N" -eq 24 ] && \
       [ "$LOCUS_N" -eq 4 ]
    then

        printf \
            "BBJ_CREDIBLE_SET_MASTER\tREADY\t\n" \
            >> "$READINESS"

    else

        printf \
            "BBJ_CREDIBLE_SET_MASTER\tCHECK\tvariants=%s loci=%s\n" \
            "$CS_N" \
            "$LOCUS_N" \
            >> "$READINESS"
    fi


    if [ -s "$GENE" ] && \
       [ -s "$PCG" ]
    then

        printf \
            "GRCH37_GENE_UNIVERSE\tREADY\t\n" \
            >> "$READINESS"

    else

        printf \
            "GRCH37_GENE_UNIVERSE\tPARTIAL\tannotation missing\n" \
            >> "$READINESS"
    fi


    if [ -s "$INV" ]; then

        printf \
            "LOCAL_FUNCTIONAL_INVENTORY\tREADY\t\n" \
            >> "$READINESS"

    else

        printf \
            "LOCAL_FUNCTIONAL_INVENTORY\tPARTIAL\tinventory missing\n" \
            >> "$READINESS"
    fi


    if [ "$KR_PASS" -eq 4 ]; then

        printf \
            "BBJ_KRIGING_AUDIT\tREADY\t\n" \
            >> "$READINESS"

    else

        printf \
            "BBJ_KRIGING_AUDIT\tPARTIAL\t%s/4\n" \
            "$KR_PASS" \
            >> "$READINESS"
    fi


    if [ -s "$SUMMARY" ]; then

        printf \
            "PHASE9A_LOCUS_SUMMARY\tREADY\t\n" \
            >> "$READINESS"

    else

        printf \
            "PHASE9A_LOCUS_SUMMARY\tBLOCKED\tmissing\n" \
            >> "$READINESS"
    fi


    printf \
        "GIGASTROKE_NAWARE\tREADY\t\n" \
        >> "$READINESS"

    printf \
        "GIGASTROKE_LD_DIAGNOSTIC\tREADY\t\n" \
        >> "$READINESS"

    printf \
        "GIGASTROKE_ANCESTRY\tUNVERIFIED\tmetadata ancestry blank\n" \
        >> "$READINESS"

    printf \
        "TPMI\tWAITING\trate_limit\n" \
        >> "$READINESS"

    printf \
        "CKB\tWAITING\tdecryption_key\n" \
        >> "$READINESS"


    echo
    echo "===== READINESS ====="

    column -t -s $'\t' \
        "$READINESS"

    return 0
}


# ============================================================
# RUN
# ============================================================

echo "===================================================="
echo "IS MASTER PHASE 9A"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="


run_step \
    "01_credible_set_master" \
    build_credible_set_master


run_step \
    "02_grch37_gene_annotation" \
    annotate_genes_grch37


run_step \
    "03_local_resource_inventory" \
    inventory_resources


run_step \
    "04_kriging_audit" \
    kriging_audit


run_step \
    "05_phase9_locus_summary" \
    build_phase9_summary


run_step \
    "06_final_status" \
    final_status


# ============================================================
# FINAL PRINT
# ============================================================

echo
echo "===================================================="
echo "FINAL STEP STATUS"
echo "===================================================="

column -t -s $'\t' \
    "$STATUS"


echo
echo "===== BBJ CREDIBLE SET MASTER ====="

column -t -s $'\t' \
    "$OUTROOT/BBJ_CREDIBLE_SET_MASTER.tsv" \
    2>/dev/null


echo
echo "===== PROTEIN-CODING GENES ====="

column -t -s $'\t' \
    "$OUTROOT/GRCH37_PROTEIN_CODING_GENES.tsv" \
    2>/dev/null


echo
echo "===== RESOURCE SUMMARY ====="

column -t -s $'\t' \
    "$OUTROOT/LOCAL_FUNCTIONAL_RESOURCE_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== KRIGING SUMMARY ====="

column -t -s $'\t' \
    "$QCROOT/kriging_fixed/BBJ_KRIGING_FIXED_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== PHASE9A LOCUS SUMMARY ====="

column -t -s $'\t' \
    "$OUTROOT/PHASE9A_LOCUS_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== FINAL READINESS ====="

column -t -s $'\t' \
    "$READINESS" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="

