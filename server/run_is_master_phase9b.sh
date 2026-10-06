#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_phase9b/$RUN_ID"

P9="$ROOT/results/is/stage5_functional/phase9"
QC="$ROOT/results/is/stage4_cross_eas/phase8_qc"
OUT="$P9/phase9b"

STATUS="$LOGROOT/STATUS.tsv"
READINESS="$ROOT/results/is/stage0_registry/IS_PHASE9B_STATUS.tsv"

mkdir -p \
  "$LOGROOT" \
  "$OUT" \
  "$QC/kriging_fixed"

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
        record "$STEP" "PASS" "completed" "$LOG"
    else
        echo "[FAIL] $STEP rc=$RC"
        record "$STEP" "FAIL" "rc=$RC" "$LOG"
    fi

    return 0
}


# ============================================================
# 01. REPAIR BBJ KRIGING AUDIT
# ============================================================

repair_kriging() {

Rscript - <<'RS'

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(susieR))

ROOT <- "/srv/is-analysis/results/is/stage3_finemap/japan/bbj"

OUT <- paste0(
    "/srv/is-analysis/results/is/",
    "stage4_cross_eas/phase8_qc/kriging_fixed"
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
    4 / (
        1/22664 +
        1/152022
    )
)


read_ld <- function(path, n) {

    expected <- n*n*8
    actual <- file.info(path)$size

    if (
        is.na(actual) ||
        actual != expected
    ) {
        stop(
            paste0(
                "LD_BYTE_MISMATCH ",
                actual,
                " != ",
                expected
            )
        )
    }

    con <- file(path, "rb")
    on.exit(close(con))

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
        "\n============================\n",
        locus,
        "\n============================\n",
        sep=""
    )

    result <- tryCatch({

        ss <- fread(
            file.path(
                ROOT,
                "susie_inputs_v3",
                paste0(locus, ".tsv")
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

        if (nrow(ss) != length(vars)) {
            stop("dimension mismatch")
        }

        if (!all(ss$variant_id == vars)) {
            stop("variant ordering mismatch")
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

        # Diagnostics were already performed before this;
        # this only guards floating-point asymmetry.
        R <- (R + t(R)) / 2
        diag(R) <- 1

        z <- (
            as.numeric(ss$beta) /
            as.numeric(ss$se)
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

        table_name <- ""
        tab <- NULL

        if (!is.null(kr$conditional_dist)) {

            table_name <- "conditional_dist"

            tab <- as.data.table(
                kr$conditional_dist
            )

        } else if (!is.null(kr$table)) {

            table_name <- "table"

            tab <- as.data.table(
                kr$table
            )
        }

        n_switch <- NA_integer_


        if (
            !is.null(tab) &&
            nrow(tab) > 0
        ) {

            cat(
                "KRIGING_ROWS=",
                nrow(tab),
                "\n",
                sep=""
            )

            cat(
                "KRIGING_COLUMNS=",
                paste(
                    names(tab),
                    collapse=","
                ),
                "\n",
                sep=""
            )

            if (nrow(tab) == length(vars)) {
                tab[, variant_id := vars]
            }

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


            lr_col <- grep(
                "logLR|log.*LR|log.*likelihood",
                names(tab),
                ignore.case=TRUE,
                value=TRUE
            )


            if (length(lr_col) >= 1) {

                lr <- suppressWarnings(
                    as.numeric(
                        tab[[lr_col[1]]]
                    )
                )

                # Use original GWAS z because its ordering
                # exactly matches LD vars.
                if (length(lr) == length(z)) {

                    switch_idx <- which(
                        is.finite(lr) &
                        lr > 2 &
                        abs(z) > 2
                    )

                    n_switch <- length(
                        switch_idx
                    )

                    suspect <- data.table(
                        variant_id=
                            vars[switch_idx],

                        z=
                            z[switch_idx],

                        logLR=
                            lr[switch_idx]
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
            locus=locus,
            status="PASS",
            n_variants=length(vars),
            n_eff=N_EFF,
            s=s,
            kriging_object=table_name,
            possible_switch_n=n_switch,
            error=""
        )

    }, error=function(e) {

        cat(
            "ERROR=",
            conditionMessage(e),
            "\n",
            sep=""
        )

        data.table(
            locus=locus,
            status="FAIL",
            n_variants=NA_integer_,
            n_eff=N_EFF,
            s=NA_real_,
            kriging_object="",
            possible_switch_n=NA_integer_,
            error=conditionMessage(e)
        )
    })


    summary_list[[length(summary_list) + 1]] <- result
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

print(summary)

RS

    return $?
}


# ============================================================
# 02. VARIANT -> GENE DISTANCE / OVERLAP
# ============================================================

map_variant_gene_distance() {

python3 - <<'PY'
import csv
from pathlib import Path

ROOT = Path("/srv/is-analysis")

CS = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "BBJ_CREDIBLE_SET_MASTER.tsv"
)

GENES = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "GRCH37_PROTEIN_CODING_GENES.tsv"
)

OUT = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "phase9b/"
    "BBJ_CS_VARIANT_GENE_DISTANCE.tsv"
)


variants = []

with CS.open() as f:

    variants = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


genes = []

with GENES.open() as f:

    genes = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


rows = []


for v in variants:

    pos = int(
        v["pos_grch37"]
    )

    locus = v["locus"]

    gs = [
        g for g in genes
        if g["locus"] == locus
    ]

    for g in gs:

        start = int(
            g["start_grch37"]
        )

        end = int(
            g["end_grch37"]
        )


        if start <= pos <= end:

            relation = "INTRAGENIC"
            distance = 0

        elif pos < start:

            relation = "UPSTREAM_COORDINATE"
            distance = start - pos

        else:

            relation = "DOWNSTREAM_COORDINATE"
            distance = pos - end


        rows.append({
            "locus":
                locus,

            "variant_id":
                v["variant_id"],

            "bbj_pip":
                v["bbj_pip"],

            "gene_symbol":
                g["gene_symbol"],

            "ensembl_gene_id":
                g["ensembl_gene_id"],

            "gene_start":
                start,

            "gene_end":
                end,

            "relation":
                relation,

            "distance_bp":
                distance
        })


rows.sort(
    key=lambda x: (
        x["locus"],
        -float(x["bbj_pip"]),
        int(x["distance_bp"])
    )
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
    w.writerows(rows)


print(
    "N_VARIANT_GENE_PAIRS=",
    len(rows)
)


for locus in sorted(
    set(
        x["locus"]
        for x in rows
    )
):

    print()
    print(
        "=====",
        locus,
        "====="
    )

    locus_rows = [
        x for x in rows
        if x["locus"] == locus
    ]

    top_variant = max(
        {
            x["variant_id"]:
                float(x["bbj_pip"])
            for x in locus_rows
        },
        key=lambda x: {
            y["variant_id"]:
                float(y["bbj_pip"])
            for y in locus_rows
        }[x]
    )

    top_rows = [
        x
        for x in locus_rows
        if x["variant_id"] == top_variant
    ]

    top_rows.sort(
        key=lambda x:
            int(x["distance_bp"])
    )

    print(
        "TOP_VARIANT=",
        top_variant
    )

    for x in top_rows[:5]:

        print(
            x["gene_symbol"],
            x["relation"],
            "DIST=",
            x["distance_bp"]
        )


print(
    "\nOUT=",
    OUT
)
PY

    return $?
}


# ============================================================
# 03. BUILD POSITIONAL CANDIDATE TABLE
# ============================================================

build_positional_candidates() {

python3 - <<'PY'
import csv
from pathlib import Path

ROOT = Path("/srv/is-analysis")

DIST = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "phase9b/"
    "BBJ_CS_VARIANT_GENE_DISTANCE.tsv"
)

CS = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "BBJ_CREDIBLE_SET_MASTER.tsv"
)

OUT = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "phase9b/"
    "POSITIONAL_GENE_CANDIDATES.tsv"
)


with DIST.open() as f:

    drows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


with CS.open() as f:

    csrows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


top_by_locus = {}

for r in csrows:

    locus = r["locus"]

    if (
        locus not in top_by_locus
        or float(r["bbj_pip"]) >
           float(
               top_by_locus[locus][
                   "bbj_pip"
               ]
           )
    ):
        top_by_locus[locus] = r


rows = []


for locus, top in top_by_locus.items():

    candidates = [
        x
        for x in drows
        if x["locus"] == locus
        and x["variant_id"] ==
            top["variant_id"]
    ]

    candidates.sort(
        key=lambda x:
            int(x["distance_bp"])
    )

    for rank, x in enumerate(
        candidates,
        start=1
    ):

        if rank > 5:
            break

        rows.append({
            "locus":
                locus,

            "bbj_top_variant":
                top["variant_id"],

            "bbj_top_pip":
                top["bbj_pip"],

            "positional_rank":
                rank,

            "gene_symbol":
                x["gene_symbol"],

            "relation":
                x["relation"],

            "distance_bp":
                x["distance_bp"],

            "evidence_level":
                (
                    "POSITIONAL_ONLY"
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
    w.writerows(rows)


print(
    "\t".join(fields)
)

for r in rows:

    print(
        "\t".join(
            str(r[x])
            for x in fields
        )
    )


print(
    "\nIMPORTANT:",
    "POSITIONAL_ONLY is not causal-gene evidence."
)
PY

    return $?
}


# ============================================================
# 04. AUDIT LOCAL eQTL / BRAIN / VASCULAR RESOURCES
# ============================================================

audit_functional_resources() {

python3 - <<'PY'
import csv
from pathlib import Path

INV = Path(
    "/srv/is-analysis/results/is/"
    "stage5_functional/phase9/"
    "LOCAL_FUNCTIONAL_RESOURCE_INVENTORY.tsv"
)

OUT = Path(
    "/srv/is-analysis/results/is/"
    "stage5_functional/phase9/"
    "phase9b/"
    "FUNCTIONAL_RESOURCE_READINESS.tsv"
)


with INV.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


categories = {
    "brain_eqtl": [],
    "vascular_eqtl": [],
    "general_eqtl": [],
    "stroke_functional": [],
    "pqtl": [],
}


for r in rows:

    p = r["path"].lower()

    if (
        "eqtl" in p and
        "brain" in p
    ):
        categories[
            "brain_eqtl"
        ].append(r)

    if (
        "eqtl" in p and
        (
            "vascular" in p
            or "artery" in p
            or "endothelial" in p
        )
    ):
        categories[
            "vascular_eqtl"
        ].append(r)

    if "eqtl" in p:
        categories[
            "general_eqtl"
        ].append(r)

    if (
        "stroke" in p
        or "/is/" in p
    ):
        categories[
            "stroke_functional"
        ].append(r)

    if (
        "pqtl" in p
        or "ukbppp" in p
    ):
        categories[
            "pqtl"
        ].append(r)


outrows = []


for category, xs in categories.items():

    outrows.append({
        "category":
            category,

        "n_files":
            len(xs),

        "ready":
            (
                "YES"
                if len(xs) > 0
                else "NO"
            ),

        "example":
            (
                xs[0]["path"]
                if xs
                else ""
            )
    })


with OUT.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "category",
            "n_files",
            "ready",
            "example"
        ]
    )

    w.writeheader()
    w.writerows(
        outrows
    )


for r in outrows:

    print(
        r["category"],
        "N=",
        r["n_files"],
        "READY=",
        r["ready"]
    )

    if r["example"]:
        print(
            " EXAMPLE=",
            r["example"]
        )
PY

    return $?
}


# ============================================================
# 05. BUILD PHASE 9B INTEGRATED LOCUS TABLE
# ============================================================

build_integrated_table() {

python3 - <<'PY'
import csv
from pathlib import Path

ROOT = Path("/srv/is-analysis")

CS = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "BBJ_CREDIBLE_SET_MASTER.tsv"
)

POS = (
    ROOT /
    "results/is/stage5_functional/phase9/"
    "phase9b/"
    "POSITIONAL_GENE_CANDIDATES.tsv"
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
    "phase9b/"
    "PHASE9B_LOCUS_MASTER.tsv"
)


with CS.open() as f:

    cs = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


with POS.open() as f:

    pos = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


kr = {}

if KR.exists():

    with KR.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):
            kr[r["locus"]] = r


loci = sorted(
    set(
        x["locus"]
        for x in cs
    )
)


rows = []


for locus in loci:

    c = [
        x
        for x in cs
        if x["locus"] == locus
    ]

    top = max(
        c,
        key=lambda x:
            float(
                x["bbj_pip"]
            )
    )

    p = [
        x
        for x in pos
        if x["locus"] == locus
    ]

    p.sort(
        key=lambda x:
            int(
                x["positional_rank"]
            )
    )


    rows.append({
        "locus":
            locus,

        "n_cs_variants":
            len(c),

        "top_variant":
            top["variant_id"],

        "top_pip":
            top["bbj_pip"],

        "position_candidate_1":
            (
                p[0]["gene_symbol"]
                if len(p) > 0
                else ""
            ),

        "position_candidate_1_relation":
            (
                p[0]["relation"]
                if len(p) > 0
                else ""
            ),

        "position_candidate_1_distance":
            (
                p[0]["distance_bp"]
                if len(p) > 0
                else ""
            ),

        "position_candidate_2":
            (
                p[1]["gene_symbol"]
                if len(p) > 1
                else ""
            ),

        "position_candidate_3":
            (
                p[2]["gene_symbol"]
                if len(p) > 2
                else ""
            ),

        "ld_s":
            kr.get(
                locus,
                {}
            ).get(
                "s",
                ""
            ),

        "possible_switch_n":
            kr.get(
                locus,
                {}
            ).get(
                "possible_switch_n",
                ""
            ),

        "causal_gene_status":
            "UNRESOLVED_NEEDS_MOLECULAR_COLOC"
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
    "\t".join(fields)
)

for r in rows:

    print(
        "\t".join(
            str(r[x])
            for x in fields
        )
    )
PY

    return $?
}


# ============================================================
# 06. FINAL STATUS
# ============================================================

final_status() {

    KR="$QC/kriging_fixed/BBJ_KRIGING_FIXED_SUMMARY.tsv"
    DIST="$OUT/BBJ_CS_VARIANT_GENE_DISTANCE.tsv"
    POS="$OUT/POSITIONAL_GENE_CANDIDATES.tsv"
    RES="$OUT/FUNCTIONAL_RESOURCE_READINESS.tsv"
    MASTER="$OUT/PHASE9B_LOCUS_MASTER.tsv"

    KR_PASS=0

    if [ -s "$KR" ]; then

        KR_PASS="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {n++}
                END {print n+0}
            ' "$KR"
        )"
    fi


    printf \
        "component\tstatus\tblocker\n" \
        > "$READINESS"


    if [ "$KR_PASS" -eq 4 ]; then
        printf "BBJ_KRIGING_AUDIT\tREADY\t\n" >> "$READINESS"
    else
        printf "BBJ_KRIGING_AUDIT\tPARTIAL\t%s/4\n" "$KR_PASS" >> "$READINESS"
    fi


    [ -s "$DIST" ] \
        && printf "VARIANT_GENE_DISTANCE\tREADY\t\n" >> "$READINESS" \
        || printf "VARIANT_GENE_DISTANCE\tBLOCKED\tmissing\n" >> "$READINESS"


    [ -s "$POS" ] \
        && printf "POSITIONAL_CANDIDATES\tREADY\t\n" >> "$READINESS" \
        || printf "POSITIONAL_CANDIDATES\tBLOCKED\tmissing\n" >> "$READINESS"


    [ -s "$RES" ] \
        && printf "FUNCTIONAL_RESOURCE_AUDIT\tREADY\t\n" >> "$READINESS" \
        || printf "FUNCTIONAL_RESOURCE_AUDIT\tBLOCKED\tmissing\n" >> "$READINESS"


    [ -s "$MASTER" ] \
        && printf "PHASE9B_LOCUS_MASTER\tREADY\t\n" >> "$READINESS" \
        || printf "PHASE9B_LOCUS_MASTER\tBLOCKED\tmissing\n" >> "$READINESS"


    printf \
        "MOLECULAR_COLOC\tPENDING\tbrain/vascular molecular QTL required\n" \
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


    column -t -s $'\t' \
        "$READINESS"

    return 0
}


# ============================================================
# RUN
# ============================================================

echo "===================================================="
echo "IS MASTER PHASE 9B"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="


run_step \
    "01_repair_kriging" \
    repair_kriging


run_step \
    "02_variant_gene_distance" \
    map_variant_gene_distance


run_step \
    "03_positional_candidates" \
    build_positional_candidates


run_step \
    "04_functional_resource_audit" \
    audit_functional_resources


run_step \
    "05_integrated_locus_table" \
    build_integrated_table


run_step \
    "06_final_status" \
    final_status


echo
echo "===================================================="
echo "FINAL STEP STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"


echo
echo "===== KRIGING ====="

column -t -s $'\t' \
    "$QC/kriging_fixed/BBJ_KRIGING_FIXED_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== POSITIONAL CANDIDATES ====="

column -t -s $'\t' \
    "$OUT/POSITIONAL_GENE_CANDIDATES.tsv" \
    2>/dev/null


echo
echo "===== FUNCTIONAL RESOURCE READINESS ====="

column -t -s $'\t' \
    "$OUT/FUNCTIONAL_RESOURCE_READINESS.tsv" \
    2>/dev/null


echo
echo "===== PHASE9B LOCUS MASTER ====="

column -t -s $'\t' \
    "$OUT/PHASE9B_LOCUS_MASTER.tsv" \
    2>/dev/null


echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$READINESS" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="
