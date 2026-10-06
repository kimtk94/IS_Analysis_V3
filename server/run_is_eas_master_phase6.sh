#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase6/$RUN_ID"
STATUS="$LOGROOT/STATUS.tsv"

mkdir -p "$LOGROOT"
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

    RC=${PIPESTATUS[0]}

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
# 01. INPUT / SAMPLE SIZE AUDIT
# ============================================================

input_audit() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    echo "===== BBJ SAMPLE SIZE ====="

    CASES=22664
    CONTROLS=152022
    TOTAL=$((CASES + CONTROLS))

    python3 - <<'PY'
cases = 22664
controls = 152022

n_total = cases + controls

n_eff = 4.0 / (
    1.0/cases +
    1.0/controls
)

print("CASES =", cases)
print("CONTROLS =", controls)
print("N_TOTAL =", n_total)
print("N_EFF =", n_eff)
print("N_EFF_ROUNDED =", round(n_eff))
PY

    echo
    echo "===== SUSIE INPUTS ====="

    FAIL=0

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        SS="$BASE/susie_inputs_v3/${LOCUS}.tsv"
        LD="$BASE/ld_v3/${LOCUS}.unphased.vcor1.bin"
        VARS="$BASE/ld_v3/${LOCUS}.unphased.vcor1.bin.vars"

        echo
        echo "$LOCUS"

        for F in "$SS" "$LD" "$VARS"
        do

            if [ -s "$F" ]; then
                echo "OK $F"
            else
                echo "MISSING $F"
                FAIL=$((FAIL + 1))
            fi

        done

    done

    echo
    echo "===== GIGASTROKE ====="

    GDIR="$ROOT/data/is/processed/gigastroke/eas"

    for SPEC in \
        "GCST90104544 AS" \
        "GCST90104545 AIS" \
        "GCST90104546 CES" \
        "GCST90104547 LAS" \
        "GCST90104548 SVS"
    do

        ACC="$(echo "$SPEC" | awk '{print $1}')"
        P="$(echo "$SPEC" | awk '{print $2}')"

        F="$GDIR/${ACC}_${P}_GRCh37.canonical.tsv.gz"

        if [ -s "$F" ]; then
            echo "OK $P $F"
        else
            echo "MISSING $P $F"
            FAIL=$((FAIL + 1))
        fi

    done

    if [ "$FAIL" -gt 0 ]; then
        return 11
    fi

    return 0
}


# ============================================================
# 02. BBJ SUSIE N_EFF SENSITIVITY
# ============================================================

run_bbj_neff_susie() {

    Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

suppressPackageStartupMessages(
    library(susieR)
)

root <- paste0(
    "/srv/is-analysis/results/is/",
    "stage3_finemap/japan/bbj"
)

loci <- c(
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
)

N_CASE <- 22664
N_CTRL <- 152022

N_TOTAL <- N_CASE + N_CTRL

N_EFF <- round(
    4 /
    (
        1/N_CASE +
        1/N_CTRL
    )
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


read_ld <- function(path, n) {

    expected <- n*n*8

    actual <- file.info(
        path
    )$size

    if (
        is.na(actual) ||
        actual != expected
    ) {
        stop(
            paste0(
                "LD byte mismatch: ",
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


master <- list()


for (locus in loci) {

    cat(
        "\n========================================\n",
        locus,
        "\n========================================\n",
        sep=""
    )

    result <- tryCatch({

        ssfile <- file.path(
            root,
            "susie_inputs_v3",
            paste0(
                locus,
                ".tsv"
            )
        )

        varsfile <- file.path(
            root,
            "ld_v3",
            paste0(
                locus,
                ".unphased.vcor1.bin.vars"
            )
        )

        ldfile <- file.path(
            root,
            "ld_v3",
            paste0(
                locus,
                ".unphased.vcor1.bin"
            )
        )

        if (
            !file.exists(ssfile) ||
            !file.exists(varsfile) ||
            !file.exists(ldfile)
        ) {
            stop(
                "input missing"
            )
        }

        ss <- fread(
            ssfile
        )

        vars <- readLines(
            varsfile
        )

        nvar <- length(
            vars
        )

        if (
            nrow(ss) != nvar
        ) {
            stop(
                "dimension mismatch"
            )
        }

        if (
            !all(
                ss$variant_id == vars
            )
        ) {
            stop(
                "variant order mismatch"
            )
        }

        R <- read_ld(
            ldfile,
            nvar
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

        fit <- susie_rss(
            z=z,
            R=R,
            n=N_EFF,
            L=10,
            estimate_residual_variance=FALSE,
            coverage=0.95,
            min_abs_corr=0.5,
            max_iter=1000,
            verbose=FALSE
        )

        outdir <- file.path(
            root,
            "susie_neff",
            locus
        )

        dir.create(
            outdir,
            recursive=TRUE,
            showWarnings=FALSE
        )

        saveRDS(
            fit,
            file.path(
                outdir,
                "susie_fit.rds"
            )
        )

        pip <- data.table(
            variant_id=vars,
            beta=ss$beta,
            se=ss$se,
            z=z,
            p=ss$p,
            pip=fit$pip
        )

        setorder(
            pip,
            -pip
        )

        fwrite(
            pip,
            file.path(
                outdir,
                "PIP.tsv"
            ),
            sep="\t"
        )

        cs <- fit$sets$cs

        csout <- data.table()

        if (
            !is.null(cs) &&
            length(cs) > 0
        ) {

            for (
                i in seq_along(cs)
            ) {

                idx <- cs[[i]]

                csout <- rbind(
                    csout,
                    data.table(
                        credible_set=i,
                        variant_id=vars[idx],
                        pip=fit$pip[idx]
                    )
                )
            }

        }

        fwrite(
            csout,
            file.path(
                outdir,
                "CREDIBLE_SETS.tsv"
            ),
            sep="\t"
        )

        top <- which.max(
            fit$pip
        )

        ncs <- if (
            is.null(cs)
        ) {
            0
        } else {
            length(cs)
        }

        cat(
            "converged=",
            fit$converged,
            "\n",
            sep=""
        )

        cat(
            "n_cs=",
            ncs,
            "\n",
            sep=""
        )

        cat(
            "top_variant=",
            vars[top],
            "\n",
            sep=""
        )

        cat(
            "top_pip=",
            fit$pip[top],
            "\n",
            sep=""
        )

        list(
            locus=locus,
            status="PASS",
            converged=fit$converged,
            n_eff=N_EFF,
            nvar=nvar,
            n_cs=ncs,
            top_variant=vars[top],
            top_pip=fit$pip[top],
            error=""
        )

    }, error=function(e) {

        cat(
            "ERROR=",
            conditionMessage(e),
            "\n",
            sep=""
        )

        list(
            locus=locus,
            status="FAIL",
            converged=FALSE,
            n_eff=N_EFF,
            nvar=NA_integer_,
            n_cs=NA_integer_,
            top_variant=NA_character_,
            top_pip=NA_real_,
            error=conditionMessage(e)
        )
    })

    master[[length(master) + 1]] <- result
}


summary <- rbindlist(
    lapply(
        master,
        as.data.table
    ),
    fill=TRUE
)

fwrite(
    summary,
    file.path(
        root,
        "SUSIE_NEFF_MASTER_SUMMARY.tsv"
    ),
    sep="\t"
)

cat(
    "\n===== N_EFF SUSIE =====\n"
)

print(
    summary
)

RS

    return $?
}


# ============================================================
# 03. COMPARE N_TOTAL vs N_EFF
# ============================================================

compare_susie_models() {

    python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is/"
    "stage3_finemap/japan/bbj"
)

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004",
]

rows = []

def read_pip(path):
    out = {}
    if not path.exists():
        return out

    with path.open() as f:
        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):
            try:
                out[r["variant_id"]] = float(
                    r["pip"]
                )
            except Exception:
                pass

    return out


def read_cs(path):
    out = set()

    if not path.exists() or path.stat().st_size == 0:
        return out

    with path.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for r in reader:

            v = r.get(
                "variant_id",
                ""
            )

            if v:
                out.add(v)

    return out


for locus in loci:

    total_pip = read_pip(
        root /
        "susie_v3" /
        locus /
        "PIP.tsv"
    )

    eff_pip = read_pip(
        root /
        "susie_neff" /
        locus /
        "PIP.tsv"
    )

    total_cs = read_cs(
        root /
        "susie_v3" /
        locus /
        "CREDIBLE_SETS.tsv"
    )

    eff_cs = read_cs(
        root /
        "susie_neff" /
        locus /
        "CREDIBLE_SETS.tsv"
    )

    top_total = (
        max(
            total_pip,
            key=total_pip.get
        )
        if total_pip
        else ""
    )

    top_eff = (
        max(
            eff_pip,
            key=eff_pip.get
        )
        if eff_pip
        else ""
    )

    union = (
        total_cs |
        eff_cs
    )

    intersection = (
        total_cs &
        eff_cs
    )

    jaccard = (
        len(intersection) /
        len(union)
        if union
        else 1.0
    )

    common = (
        set(total_pip) &
        set(eff_pip)
    )

    max_abs_pip_diff = (
        max(
            abs(
                total_pip[v] -
                eff_pip[v]
            )
            for v in common
        )
        if common
        else None
    )

    rows.append({
        "locus":
            locus,

        "top_total":
            top_total,

        "top_total_pip":
            total_pip.get(
                top_total,
                ""
            ),

        "top_neff":
            top_eff,

        "top_neff_pip":
            eff_pip.get(
                top_eff,
                ""
            ),

        "same_top_variant":
            "YES"
            if top_total == top_eff
            else "NO",

        "cs_total_size":
            len(total_cs),

        "cs_neff_size":
            len(eff_cs),

        "cs_intersection":
            len(intersection),

        "cs_union":
            len(union),

        "cs_jaccard":
            jaccard,

        "max_abs_pip_diff":
            max_abs_pip_diff
    })


out = (
    root /
    "SUSIE_NTOTAL_VS_NEFF_COMPARISON.tsv"
)

fields = list(
    rows[0].keys()
)

with out.open(
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
            str(
                r[x]
            )
            for x in fields
        )
    )
PY

    return $?
}


# ============================================================
# 04. EXTRACT GIGASTROKE 5 PHENOTYPES × 4 BBJ REGIONS
# ============================================================

extract_gigastroke_regions() {

    python3 - <<'PY'
import csv
import gzip
from pathlib import Path

root = Path(
    "/srv/is-analysis"
)

indir = (
    root /
    "data/is/processed/gigastroke/eas"
)

region_file = (
    root /
    "results/is/stage3_finemap/"
    "japan/bbj/"
    "BBJ_IS_FINEMAP_REGIONS.tsv"
)

outdir = (
    root /
    "results/is/stage4_cross_eas/"
    "gigastroke_regions"
)

outdir.mkdir(
    parents=True,
    exist_ok=True
)


with region_file.open() as f:

    all_regions = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

regions = [
    r
    for r in all_regions
    if r["locus_id"] in {
        "BBJ_IS_L001",
        "BBJ_IS_L002",
        "BBJ_IS_L003",
        "BBJ_IS_L004",
    }
]


specs = {
    "AS":
        "GCST90104544_AS_GRCh37.canonical.tsv.gz",

    "AIS":
        "GCST90104545_AIS_GRCh37.canonical.tsv.gz",

    "CES":
        "GCST90104546_CES_GRCh37.canonical.tsv.gz",

    "LAS":
        "GCST90104547_LAS_GRCh37.canonical.tsv.gz",

    "SVS":
        "GCST90104548_SVS_GRCh37.canonical.tsv.gz",
}


for phenotype, filename in specs.items():

    inp = (
        indir /
        filename
    )

    print()
    print(
        "PHENOTYPE",
        phenotype,
        inp
    )

    handles = {}
    writers = {}
    counts = {}

    fields = [
        "dataset",
        "phenotype",
        "ancestry",
        "build",
        "chr",
        "pos",
        "effect_allele",
        "other_allele",
        "beta",
        "se",
        "p",
        "eaf",
        "or",
        "variant_id"
    ]


    for r in regions:

        locus = r["locus_id"]

        out = (
            outdir /
            f"{phenotype}.{locus}.tsv"
        )

        h = out.open(
            "w",
            newline=""
        )

        w = csv.DictWriter(
            h,
            delimiter="\t",
            fieldnames=fields
        )

        w.writeheader()

        handles[locus] = h
        writers[locus] = w
        counts[locus] = 0


    with gzip.open(
        inp,
        "rt"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            chrom = str(
                row["chr"]
            )

            try:
                pos = int(
                    row["pos"]
                )
            except Exception:
                continue


            for r in regions:

                if chrom != str(
                    r["chr"]
                ):
                    continue

                start = int(
                    r["region_start"]
                )

                end = int(
                    r["region_end"]
                )

                if (
                    start <= pos <= end
                ):

                    writers[
                        r["locus_id"]
                    ].writerow(
                        {
                            k:
                                row.get(
                                    k,
                                    ""
                                )
                            for k in fields
                        }
                    )

                    counts[
                        r["locus_id"]
                    ] += 1


    for h in handles.values():
        h.close()


    for locus, n in counts.items():

        print(
            phenotype,
            locus,
            "N=",
            n
        )
PY

    return $?
}


# ============================================================
# 05. HARMONIZE GIGASTROKE TO BBJ ALLELE ORIENTATION
# ============================================================

harmonize_gigastroke() {

    python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is"
)

bbj_root = (
    root /
    "stage3_finemap/japan/bbj"
)

giga_root = (
    root /
    "stage4_cross_eas/gigastroke_regions"
)

outdir = (
    root /
    "stage4_cross_eas/harmonized"
)

outdir.mkdir(
    parents=True,
    exist_ok=True
)


phenotypes = [
    "AS",
    "AIS",
    "CES",
    "LAS",
    "SVS"
]

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
]


COMP = {
    "A":"T",
    "T":"A",
    "C":"G",
    "G":"C"
}


def comp(a):

    return "".join(
        COMP.get(
            x,
            x
        )
        for x in a.upper()
    )


def palindromic(a, b):

    pair = {
        a.upper(),
        b.upper()
    }

    return pair in [
        {"A","T"},
        {"C","G"}
    ]


def harmonize(
    ref,
    alt,
    bbj_eaf,
    ea,
    oa,
    giga_eaf,
    beta
):

    ref = ref.upper()
    alt = alt.upper()

    ea = ea.upper()
    oa = oa.upper()

    bbj_eaf = float(
        bbj_eaf
    )

    giga_eaf = float(
        giga_eaf
    )

    beta = float(
        beta
    )


    # Exact orientation:
    # BBJ beta is ALT/effect allele.
    if (
        ea == alt and
        oa == ref
    ):

        return (
            "MATCH",
            beta,
            giga_eaf
        )


    if (
        ea == ref and
        oa == alt
    ):

        return (
            "SWAP",
            -beta,
            1-giga_eaf
        )


    # Complement strand
    cea = comp(ea)
    coa = comp(oa)


    if (
        cea == alt and
        coa == ref
    ):

        # Palindromic variants require
        # EAF support.
        if palindromic(
            ref,
            alt
        ):

            if abs(
                giga_eaf -
                bbj_eaf
            ) > 0.15:

                return (
                    "AMBIGUOUS_PALINDROME",
                    None,
                    None
                )

        return (
            "COMPLEMENT",
            beta,
            giga_eaf
        )


    if (
        cea == ref and
        coa == alt
    ):

        aligned_eaf = (
            1-giga_eaf
        )

        if palindromic(
            ref,
            alt
        ):

            if abs(
                aligned_eaf -
                bbj_eaf
            ) > 0.15:

                return (
                    "AMBIGUOUS_PALINDROME",
                    None,
                    None
                )

        return (
            "COMPLEMENT_SWAP",
            -beta,
            aligned_eaf
        )


    return (
        "ALLELE_MISMATCH",
        None,
        None
    )


for locus in loci:

    bbj_file = (
        bbj_root /
        "intersection" /
        f"{locus}.summary.tsv"
    )

    bbj = {}

    with bbj_file.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):

            key = (
                r["chr"],
                int(
                    r["pos"]
                )
            )

            bbj[key] = r


    for pheno in phenotypes:

        giga_file = (
            giga_root /
            f"{pheno}.{locus}.tsv"
        )

        outfile = (
            outdir /
            f"{pheno}.{locus}.harmonized.tsv"
        )

        output = []

        if not giga_file.exists():

            print(
                pheno,
                locus,
                "INPUT_MISSING"
            )

            continue


        with giga_file.open() as f:

            reader = csv.DictReader(
                f,
                delimiter="\t"
            )

            for g in reader:

                try:

                    key = (
                        g["chr"],
                        int(
                            g["pos"]
                        )
                    )

                except Exception:

                    continue


                if key not in bbj:
                    continue


                b = bbj[key]


                try:

                    status, abeta, aeaf = harmonize(
                        b["bbj_ref"],
                        b["bbj_alt"],
                        b["eaf"],
                        g["effect_allele"],
                        g["other_allele"],
                        g["eaf"],
                        g["beta"]
                    )

                except Exception:

                    status = "NUMERIC_ERROR"
                    abeta = None
                    aeaf = None


                row = {
                    "phenotype":
                        pheno,

                    "locus":
                        locus,

                    "chr":
                        b["chr"],

                    "pos":
                        b["pos"],

                    "variant_id":
                        b["variant_id"],

                    "bbj_ref":
                        b["bbj_ref"],

                    "bbj_alt":
                        b["bbj_alt"],

                    "bbj_beta":
                        b["beta"],

                    "bbj_p":
                        b["p"],

                    "bbj_eaf":
                        b["eaf"],

                    "giga_effect_allele":
                        g["effect_allele"],

                    "giga_other_allele":
                        g["other_allele"],

                    "giga_beta_raw":
                        g["beta"],

                    "giga_beta_aligned":
                        ""
                        if abeta is None
                        else abeta,

                    "giga_p":
                        g["p"],

                    "giga_eaf_raw":
                        g["eaf"],

                    "giga_eaf_aligned":
                        ""
                        if aeaf is None
                        else aeaf,

                    "harmonization":
                        status
                }

                output.append(
                    row
                )


        fields = [
            "phenotype",
            "locus",
            "chr",
            "pos",
            "variant_id",
            "bbj_ref",
            "bbj_alt",
            "bbj_beta",
            "bbj_p",
            "bbj_eaf",
            "giga_effect_allele",
            "giga_other_allele",
            "giga_beta_raw",
            "giga_beta_aligned",
            "giga_p",
            "giga_eaf_raw",
            "giga_eaf_aligned",
            "harmonization"
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
                output
            )


        from collections import Counter

        c = Counter(
            x["harmonization"]
            for x in output
        )

        print(
            pheno,
            locus,
            "N=",
            len(output),
            dict(c)
        )

PY

    return $?
}


# ============================================================
# 06. BUILD 20-COMPARISON BENCHMARK SUMMARY
# ============================================================

build_cross_eas_summary() {

    python3 - <<'PY'
import csv
import math
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is"
)

bbj_root = (
    root /
    "stage3_finemap/japan/bbj"
)

harm_root = (
    root /
    "stage4_cross_eas/harmonized"
)

outroot = (
    root /
    "stage4_cross_eas"
)

outroot.mkdir(
    parents=True,
    exist_ok=True
)


phenotypes = [
    "AS",
    "AIS",
    "CES",
    "LAS",
    "SVS"
]

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
]


# Load BBJ total-N SuSiE PIP and credible sets
pip = {}
credible = {}


for locus in loci:

    pip_file = (
        bbj_root /
        "susie_v3" /
        locus /
        "PIP.tsv"
    )

    cs_file = (
        bbj_root /
        "susie_v3" /
        locus /
        "CREDIBLE_SETS.tsv"
    )

    pip[locus] = {}

    with pip_file.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):

            pip[locus][
                r["variant_id"]
            ] = float(
                r["pip"]
            )


    credible[locus] = set()

    if (
        cs_file.exists() and
        cs_file.stat().st_size > 0
    ):

        with cs_file.open() as f:

            for r in csv.DictReader(
                f,
                delimiter="\t"
            ):

                v = r.get(
                    "variant_id",
                    ""
                )

                if v:
                    credible[
                        locus
                    ].add(
                        v
                    )


summary = []
details = []


valid_harmonization = {
    "MATCH",
    "SWAP",
    "COMPLEMENT",
    "COMPLEMENT_SWAP"
}


for pheno in phenotypes:

    for locus in loci:

        f = (
            harm_root /
            f"{pheno}.{locus}.harmonized.tsv"
        )

        rows = []

        if f.exists():

            with f.open() as fh:

                rows = list(
                    csv.DictReader(
                        fh,
                        delimiter="\t"
                    )
                )


        valid = [
            r for r in rows
            if r[
                "harmonization"
            ] in valid_harmonization
        ]


        for r in valid:

            r["bbj_pip"] = pip[
                locus
            ].get(
                r["variant_id"],
                0.0
            )

            r["in_bbj_cs"] = (
                "YES"
                if r[
                    "variant_id"
                ] in credible[
                    locus
                ]
                else "NO"
            )

            try:

                bb = float(
                    r["bbj_beta"]
                )

                gb = float(
                    r[
                        "giga_beta_aligned"
                    ]
                )

                r[
                    "direction_concordant"
                ] = (
                    "YES"
                    if (
                        bb > 0
                        and gb > 0
                    ) or (
                        bb < 0
                        and gb < 0
                    )
                    else "NO"
                )

            except Exception:

                r[
                    "direction_concordant"
                ] = "NA"


            details.append(
                r
            )


        lead = None

        if valid:

            lead = min(
                valid,
                key=lambda x:
                    float(
                        x["giga_p"]
                    )
            )


        bbj_top = None

        if pip[locus]:

            bbj_top = max(
                pip[locus],
                key=pip[locus].get
            )


        bbj_top_row = next(
            (
                r
                for r in valid
                if r[
                    "variant_id"
                ] == bbj_top
            ),
            None
        )


        cs_rows = [
            r for r in valid
            if r[
                "variant_id"
            ] in credible[
                locus
            ]
        ]


        pip01_rows = [
            r for r in valid
            if pip[
                locus
            ].get(
                r["variant_id"],
                0
            ) >= 0.1
        ]


        concord_cs = [
            r for r in cs_rows
            if r[
                "direction_concordant"
            ] == "YES"
        ]


        concord_pip01 = [
            r for r in pip01_rows
            if r[
                "direction_concordant"
            ] == "YES"
        ]


        best_cs_p = (
            min(
                float(
                    r["giga_p"]
                )
                for r in cs_rows
            )
            if cs_rows
            else None
        )


        summary.append({
            "phenotype":
                pheno,

            "locus":
                locus,

            "analysis_role":
                "EAS_CONSORTIUM_BENCHMARK",

            "n_position_overlap":
                len(rows),

            "n_harmonized":
                len(valid),

            "harmonization_rate":
                (
                    len(valid) /
                    len(rows)
                    if rows
                    else 0
                ),

            "giga_locus_lead":
                (
                    lead[
                        "variant_id"
                    ]
                    if lead
                    else ""
                ),

            "giga_locus_lead_p":
                (
                    lead[
                        "giga_p"
                    ]
                    if lead
                    else ""
                ),

            "bbj_top_variant":
                bbj_top or "",

            "bbj_top_pip":
                (
                    pip[locus].get(
                        bbj_top,
                        ""
                    )
                    if bbj_top
                    else ""
                ),

            "bbj_top_present_in_giga":
                (
                    "YES"
                    if bbj_top_row
                    else "NO"
                ),

            "bbj_top_giga_p":
                (
                    bbj_top_row[
                        "giga_p"
                    ]
                    if bbj_top_row
                    else ""
                ),

            "bbj_top_giga_beta_aligned":
                (
                    bbj_top_row[
                        "giga_beta_aligned"
                    ]
                    if bbj_top_row
                    else ""
                ),

            "bbj_top_direction_concordant":
                (
                    bbj_top_row[
                        "direction_concordant"
                    ]
                    if bbj_top_row
                    else ""
                ),

            "bbj_cs_size":
                len(
                    credible[
                        locus
                    ]
                ),

            "bbj_cs_present_in_giga":
                len(
                    cs_rows
                ),

            "bbj_cs_best_giga_p":
                (
                    best_cs_p
                    if best_cs_p is not None
                    else ""
                ),

            "bbj_cs_direction_concordant_n":
                len(
                    concord_cs
                ),

            "bbj_cs_direction_concordance":
                (
                    len(
                        concord_cs
                    ) /
                    len(
                        cs_rows
                    )
                    if cs_rows
                    else ""
                ),

            "pip_ge_0.1_present_n":
                len(
                    pip01_rows
                ),

            "pip_ge_0.1_direction_concordant_n":
                len(
                    concord_pip01
                ),

            "pip_ge_0.1_direction_concordance":
                (
                    len(
                        concord_pip01
                    ) /
                    len(
                        pip01_rows
                    )
                    if pip01_rows
                    else ""
                )
        })


summary_file = (
    outroot /
    "GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv"
)

fields = list(
    summary[0].keys()
)

with summary_file.open(
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
        summary
    )


detail_file = (
    outroot /
    "GIGASTROKE_BBJ_HARMONIZED_VARIANTS.tsv"
)

if details:

    dfields = list(
        details[0].keys()
    )

    with detail_file.open(
        "w",
        newline=""
    ) as f:

        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=dfields
        )

        w.writeheader()
        w.writerows(
            details
        )


print(
    "N_COMPARISONS=",
    len(summary)
)

print()
print(
    "\t".join(
        fields
    )
)

for r in summary:

    print(
        "\t".join(
            str(
                r[x]
            )
            for x in fields
        )
    )
PY

    return $?
}


# ============================================================
# 07. CREATE LOCUS × PHENOTYPE MATRIX
# ============================================================

build_matrix() {

    python3 - <<'PY'
import csv
from pathlib import Path

inp = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/"
    "GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv"
)

out = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/"
    "GIGASTROKE_BBJ_TOP_P_MATRIX.tsv"
)

with inp.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


phenotypes = [
    "AS",
    "AIS",
    "CES",
    "LAS",
    "SVS"
]

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
]


index = {
    (
        r["locus"],
        r["phenotype"]
    ):
        r
    for r in rows
}


with out.open(
    "w",
    newline=""
) as f:

    fields = [
        "locus"
    ] + phenotypes

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    w.writeheader()

    for locus in loci:

        row = {
            "locus":
                locus
        }

        for p in phenotypes:

            x = index.get(
                (
                    locus,
                    p
                ),
                {}
            )

            row[p] = x.get(
                "bbj_top_giga_p",
                ""
            )

        w.writerow(
            row
        )


print(
    out
)

print(
    out.read_text()
)
PY

    return $?
}


# ============================================================
# 08. FINAL STATUS
# ============================================================

final_status() {

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE6_STATUS.tsv"

    NEFF="$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_NEFF_MASTER_SUMMARY.tsv"

    COMP="$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_NTOTAL_VS_NEFF_COMPARISON.tsv"

    BENCH="$ROOT/results/is/stage4_cross_eas/GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv"

    NEFF_PASS=0

    if [ -s "$NEFF" ]; then

        NEFF_PASS="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {n++}
                END {print n+0}
            ' "$NEFF"
        )"

    fi

    BENCH_N=0

    if [ -s "$BENCH" ]; then

        BENCH_N="$(
            awk '
                END {
                    print NR > 0 ? NR-1 : 0
                }
            ' "$BENCH"
        )"

    fi


    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"


    if [ "$NEFF_PASS" -eq 4 ]; then

        printf \
            "BBJ_NEFF_SUSIE\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "BBJ_NEFF_SUSIE\tPARTIAL\t%s/4\n" \
            "$NEFF_PASS" \
            >> "$OUT"

    fi


    if [ -s "$COMP" ]; then

        printf \
            "BBJ_N_SENSITIVITY\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "BBJ_N_SENSITIVITY\tBLOCKED\tcomparison_missing\n" \
            >> "$OUT"

    fi


    if [ "$BENCH_N" -eq 20 ]; then

        printf \
            "GIGASTROKE_20_REGION_BENCHMARK\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "GIGASTROKE_20_REGION_BENCHMARK\tPARTIAL\t%s/20\n" \
            "$BENCH_N" \
            >> "$OUT"

    fi


    printf \
        "BBJ_PRIMARY_FINEMAP\tREADY\t\n" \
        >> "$OUT"

    printf \
        "TPMI\tWAITING\trate_limit\n" \
        >> "$OUT"

    printf \
        "CKB\tWAITING\tdecryption_key\n" \
        >> "$OUT"


    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 6"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="


run_step \
    "01_input_audit" \
    input_audit


run_step \
    "02_bbj_neff_susie" \
    run_bbj_neff_susie


run_step \
    "03_compare_susie_models" \
    compare_susie_models


run_step \
    "04_extract_gigastroke_regions" \
    extract_gigastroke_regions


run_step \
    "05_harmonize_gigastroke" \
    harmonize_gigastroke


run_step \
    "06_build_cross_eas_summary" \
    build_cross_eas_summary


run_step \
    "07_build_matrix" \
    build_matrix


run_step \
    "08_final_status" \
    final_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"


echo
echo "===== N_TOTAL VS N_EFF ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_NTOTAL_VS_NEFF_COMPARISON.tsv" \
    2>/dev/null


echo
echo "===== GIGASTROKE 20-COMPARISON ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== BBJ TOP VARIANT P MATRIX ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/GIGASTROKE_BBJ_TOP_P_MATRIX.tsv" \
    2>/dev/null


echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE6_STATUS.tsv" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="

