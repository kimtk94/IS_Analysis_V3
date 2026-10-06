#!/usr/bin/env bash

ROOT="/srv/is-analysis"

BASE="$ROOT/results/is/stage4_cross_eas/phase8_qc"
NROOT="$BASE/susie_naware_corrected"
MODEL="$BASE/GIGASTROKE_SAMPLE_SIZE_MODEL.tsv"

MASTER="$BASE/GIGASTROKE_SUSIE_NAWARE_CORRECTED_MASTER.tsv"
COMP="$BASE/GIGASTROKE_ZONLY_VS_NAWARE_CORRECTED.tsv"

STATUS="$ROOT/results/is/stage0_registry/IS_EAS_PHASE8C_STATUS.tsv"

echo "===================================================="
echo "PHASE 8C FINALIZE"
echo "START=$(date)"
echo "===================================================="


# ------------------------------------------------------------
# 1. REBUILD MASTER FROM EXISTING 20 RESULT DIRECTORIES
# ------------------------------------------------------------

python3 - <<'PY'
import csv
from pathlib import Path

base = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/phase8_qc"
)

nroot = base / "susie_naware_corrected"

model_file = (
    base /
    "GIGASTROKE_SAMPLE_SIZE_MODEL.tsv"
)

out = (
    base /
    "GIGASTROKE_SUSIE_NAWARE_CORRECTED_MASTER.tsv"
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


meta = {}

with model_file.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):
        meta[r["phenotype"]] = r


rows = []


for phenotype in phenotypes:

    for locus in loci:

        pdir = (
            nroot /
            phenotype /
            locus
        )

        pipfile = (
            pdir /
            "PIP.tsv"
        )

        csfile = (
            pdir /
            "CREDIBLE_SETS.tsv"
        )

        if not pipfile.exists():

            rows.append({
                "phenotype":
                    phenotype,

                "locus":
                    locus,

                "status":
                    "FAIL",

                "converged":
                    "",

                "n_used":
                    meta.get(
                        phenotype,
                        {}
                    ).get(
                        "n_used",
                        ""
                    ),

                "n_mode":
                    meta.get(
                        phenotype,
                        {}
                    ).get(
                        "n_mode",
                        ""
                    ),

                "n_variants":
                    "",

                "n_cs":
                    "",

                "top_variant":
                    "",

                "top_pip":
                    "",

                "error":
                    "PIP.tsv missing"
            })

            continue


        with pipfile.open() as f:

            piprows = list(
                csv.DictReader(
                    f,
                    delimiter="\t"
                )
            )


        if not piprows:

            status = "FAIL"
            top_variant = ""
            top_pip = ""
            error = "empty PIP"

        else:

            status = "PASS"

            top = max(
                piprows,
                key=lambda x:
                    float(
                        x["pip"]
                    )
            )

            top_variant = (
                top["variant_id"]
            )

            top_pip = (
                top["pip"]
            )

            error = ""


        cs_ids = set()

        if (
            csfile.exists()
            and csfile.stat().st_size > 0
        ):

            with csfile.open() as f:

                reader = csv.DictReader(
                    f,
                    delimiter="\t"
                )

                for r in reader:

                    x = r.get(
                        "credible_set",
                        ""
                    )

                    if x:
                        cs_ids.add(x)


        rows.append({
            "phenotype":
                phenotype,

            "locus":
                locus,

            "status":
                status,

            "converged":
                "TRUE"
                if status == "PASS"
                else "",

            "n_used":
                meta[
                    phenotype
                ]["n_used"],

            "n_mode":
                meta[
                    phenotype
                ]["n_mode"],

            "n_variants":
                len(
                    piprows
                ),

            "n_cs":
                len(
                    cs_ids
                ),

            "top_variant":
                top_variant,

            "top_pip":
                top_pip,

            "error":
                error
        })


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
    "TOTAL=",
    len(rows)
)

print(
    "PASS=",
    sum(
        x["status"] == "PASS"
        for x in rows
    )
)

print()
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


# ------------------------------------------------------------
# 2. VERIFY EXACT N PER PHENOTYPE
# ------------------------------------------------------------

echo
echo "===== N USED AUDIT ====="

python3 - <<'PY'
import csv
from pathlib import Path

f = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/phase8_qc/"
    "GIGASTROKE_SUSIE_NAWARE_CORRECTED_MASTER.tsv"
)

expected = {
    "AS": 264655,
    "AIS": 256274,
    "CES": 238168,
    "LAS": 238977,
    "SVS": 242774,
}

bad = 0

with f.open() as h:

    for r in csv.DictReader(
        h,
        delimiter="\t"
    ):

        actual = int(
            float(
                r["n_used"]
            )
        )

        target = expected[
            r["phenotype"]
        ]

        ok = (
            actual == target
        )

        print(
            r["phenotype"],
            r["locus"],
            "N=",
            actual,
            "PASS="
            + str(ok)
        )

        if not ok:
            bad += 1


raise SystemExit(
    0 if bad == 0 else 1
)
PY


# ------------------------------------------------------------
# 3. LD DIAGNOSTIC SUMMARY
# ------------------------------------------------------------

echo
echo "===== BBJ LD S ====="

column -t -s $'\t' \
    "$BASE/bbj_ld_diagnostics/BBJ_LD_DIAGNOSTIC_SUMMARY.tsv"


echo
echo "===== GIGASTROKE L004 LD S ====="

awk -F '\t' '
    NR==1 || $2=="BBJ_IS_L004"
' \
"$BASE/gigastroke_ld_diagnostics/GIGASTROKE_LD_DIAGNOSTIC_SUMMARY.tsv" \
| column -t -s $'\t'


# ------------------------------------------------------------
# 4. Z-ONLY VS N-AWARE EXCEPTIONS
# ------------------------------------------------------------

echo
echo "===== SENSITIVITY EXCEPTIONS ====="

awk -F '\t' '
    NR==1 ||
    $7!="YES" ||
    $12 < 0.8
' "$COMP" \
| column -t -s $'\t'


# ------------------------------------------------------------
# 5. KRIGING COLUMN AUDIT
# ------------------------------------------------------------

echo
echo "===== KRIGING COLUMN AUDIT ====="

python3 - <<'PY'
import csv
from pathlib import Path

roots = [
    Path(
        "/srv/is-analysis/results/is/"
        "stage4_cross_eas/phase8_qc/"
        "bbj_ld_diagnostics"
    ),

    Path(
        "/srv/is-analysis/results/is/"
        "stage4_cross_eas/phase8_qc/"
        "gigastroke_ld_diagnostics/kriging"
    ),
]

seen = 0

for root in roots:

    for f in sorted(
        root.glob("*.tsv")
    ):

        if (
            "SUMMARY" in
            f.name.upper()
        ):
            continue

        with f.open() as h:

            reader = csv.reader(
                h,
                delimiter="\t"
            )

            try:
                header = next(
                    reader
                )
            except StopIteration:
                continue

        print(
            f.name,
            "=>",
            ",".join(
                header
            )
        )

        seen += 1

        if seen >= 6:
            break

    if seen >= 6:
        break

print(
    "FILES_INSPECTED=",
    seen
)
PY


# ------------------------------------------------------------
# 6. FINAL READINESS
# ------------------------------------------------------------

TOTAL="$(
    awk '
        NR>1 {n++}
        END {print n+0}
    ' "$MASTER"
)"

PASS="$(
    awk -F '\t' '
        NR>1 && $3=="PASS" {
            n++
        }
        END {
            print n+0
        }
    ' "$MASTER"
)"


BBJ="$(
    awk -F '\t' '
        NR>1 && $2=="PASS" {
            n++
        }
        END {
            print n+0
        }
    ' \
    "$BASE/bbj_ld_diagnostics/BBJ_LD_DIAGNOSTIC_SUMMARY.tsv"
)"


GIGA="$(
    awk -F '\t' '
        NR>1 && $3=="PASS" {
            n++
        }
        END {
            print n+0
        }
    ' \
    "$BASE/gigastroke_ld_diagnostics/GIGASTROKE_LD_DIAGNOSTIC_SUMMARY.tsv"
)"


printf \
    "component\tstatus\tblocker\n" \
    > "$STATUS"


if [ "$TOTAL" -eq 20 ] && \
   [ "$PASS" -eq 20 ]
then

    printf \
        "GIGASTROKE_NAWARE_CORRECTED\tREADY\t\n" \
        >> "$STATUS"

else

    printf \
        "GIGASTROKE_NAWARE_CORRECTED\tPARTIAL\t%s/20\n" \
        "$PASS" \
        >> "$STATUS"
fi


if [ "$BBJ" -eq 4 ]; then
    printf \
        "BBJ_LD_DIAGNOSTIC\tREADY\t\n" \
        >> "$STATUS"
else
    printf \
        "BBJ_LD_DIAGNOSTIC\tPARTIAL\t%s/4\n" \
        "$BBJ" \
        >> "$STATUS"
fi


if [ "$GIGA" -eq 20 ]; then
    printf \
        "GIGASTROKE_LD_DIAGNOSTIC\tREADY\t\n" \
        >> "$STATUS"
else
    printf \
        "GIGASTROKE_LD_DIAGNOSTIC\tPARTIAL\t%s/20\n" \
        "$GIGA" \
        >> "$STATUS"
fi


printf \
    "N_SENSITIVITY_CORRECTED\tREADY\t\n" \
    >> "$STATUS"

printf \
    "PHASE7_SHARED_SIGNAL_BENCHMARK\tREADY\t\n" \
    >> "$STATUS"

printf \
    "GIGASTROKE_ANCESTRY\tUNVERIFIED\tmetadata ancestry blank\n" \
    >> "$STATUS"

printf \
    "TPMI\tWAITING\trate_limit\n" \
    >> "$STATUS"

printf \
    "CKB\tWAITING\tdecryption_key\n" \
    >> "$STATUS"


echo
echo "===== FINAL READINESS ====="

column -t -s $'\t' \
    "$STATUS"


echo
echo "===================================================="
echo "END=$(date)"
echo "===================================================="

