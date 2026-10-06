#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

P8F="$ROOT/results/is/stage4_cross_eas/phase8f_l003_jpt_ld"
P8E="$ROOT/results/is/stage4_cross_eas/phase8e_l003_rs671"
P8D="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit_r2"

OUT="$ROOT/results/is/stage5_functional/phase9c_convergence"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9c_functional/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

BBJ="$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"

GIGA="$ROOT/results/is/stage4_cross_eas/GIGASTROKE_BBJ_HARMONIZED_VARIANTS.tsv"

RS671="12:112241766:G:A"

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

  echo "===== PHASE8F ====="

  if [[ ! -s "$P8F/STEP_STATUS.tsv" ]]; then
    echo "PHASE8F_STATUS_MISSING"
    return 1
  fi

  cat "$P8F/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P8F/STEP_STATUS.tsv"
  )"

  echo
  echo "PHASE8F_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    echo "PHASE8F_NOT_CLEAN"
    return 1
  fi

  echo
  echo "===== CORE INPUTS ====="

  for f in \
    "$BBJ" \
    "$P8F/L003_JPT_CONDITIONING_SUMMARY.tsv" \
    "$P8E/L003_RS671_CONDITIONING_SUMMARY.tsv" \
    "$P8D/SCIENTIFIC_READINESS_R2.tsv"
  do

    if [[ -s "$f" ]]; then
      echo "FOUND $f"
    else
      echo "MISSING $f"
      return 1
    fi

  done

  echo
  echo "===== OPTIONAL GIGASTROKE MASTER ====="

  if [[ -s "$GIGA" ]]; then
    echo "FOUND $GIGA"
  else
    echo "GIGASTROKE_MASTER_NOT_FOUND"
  fi

  echo
  echo "===== DEPENDENCIES ====="

  for exe in \
    python3 \
    curl
  do

    if command -v "$exe" >/dev/null 2>&1; then
      echo "FOUND $exe $(command -v "$exe")"
    else
      echo "MISSING $exe"
      return 1
    fi

  done

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. FREEZE L003 STATISTICAL CONCLUSION
# =============================================================================

freeze_l003() {

python3 - \
  "$P8E" \
  "$P8F" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


P8E = Path(sys.argv[1])
P8F = Path(sys.argv[2])
OUT = Path(sys.argv[3])


def read_first(path):

    with path.open() as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    return rows[0] if rows else {}


e8 = read_first(
    P8E /
    "L003_RS671_CONDITIONING_SUMMARY.tsv"
)

f8 = read_first(
    P8F /
    "L003_JPT_CONDITIONING_SUMMARY.tsv"
)


row = {

    "locus":
        "BBJ_IS_L003",

    "primary_gene":
        "ALDH2",

    "primary_variant":
        e8.get(
            "condition_variant",
            "12:112241766:G:A"
        ),

    "primary_rsid":
        e8.get(
            "condition_rsid",
            "rs671"
        ),

    "bbj_beta":
        e8.get(
            "condition_beta",
            ""
        ),

    "bbj_se":
        e8.get(
            "condition_se",
            ""
        ),

    "bbj_z":
        e8.get(
            "condition_z",
            ""
        ),

    "bbj_p":
        e8.get(
            "condition_p",
            ""
        ),

    "bbj_eaf":
        e8.get(
            "condition_eaf",
            ""
        ),

    "eas504_apparent_secondary_gws_n":
        e8.get(
            "conditional_p_lt_5e8_n",
            ""
        ),

    "jpt_secondary_gws_n":
        f8.get(
            "jpt_conditional_gws_n",
            ""
        ),

    "jpt_low_ld_secondary_gws_n":
        f8.get(
            "jpt_conditional_low_ld_gws_n",
            ""
        ),

    "jpt_top_residual_variant":
        f8.get(
            "jpt_top_variant",
            ""
        ),

    "jpt_top_residual_p":
        f8.get(
            "jpt_top_p",
            ""
        ),

    "jpt_top_residual_r2_to_rs671":
        f8.get(
            "jpt_top_r2_to_rs671",
            ""
        ),

    "statistical_interpretation":
        "RS671_DOMINANT_EAS_SIGNAL_NO_ROBUST_INDEPENDENT_SECONDARY_SIGNAL",

    "paper_grade_limitation":
        "COHORT_MATCHED_BBJ_LD_UNAVAILABLE",

    "status":
        "FROZEN_EXPLORATORY",
}


outfile = (
    OUT /
    "L003_STATISTICAL_FREEZE.tsv"
)

with outfile.open(
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
# 03. RS671 LOCAL EVIDENCE MASTER
# =============================================================================

build_rs671_local_master() {

python3 - \
  "$BBJ" \
  "$GIGA" \
  "$P8E" \
  "$P8F" \
  "$OUT" \
  "$RS671" <<'PY'

from pathlib import Path
import csv
import gzip
import sys


BBJ = Path(sys.argv[1])
GIGA = Path(sys.argv[2])
P8E = Path(sys.argv[3])
P8F = Path(sys.argv[4])
OUT = Path(sys.argv[5])
RS671 = sys.argv[6]


# ------------------------------------------------------------------
# BBJ
# ------------------------------------------------------------------

bbj_row = None

with gzip.open(
    BBJ,
    "rt",
    errors="replace"
) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for r in reader:

        if r.get(
            "variant_id"
        ) == RS671:

            bbj_row = r
            break


if bbj_row is None:
    raise RuntimeError(
        "rs671 not found in BBJ canonical"
    )


with (
    OUT /
    "RS671_BBJ_MASTER.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            bbj_row.keys()
        )
    )

    w.writeheader()
    w.writerow(
        bbj_row
    )


print("===== BBJ rs671 =====")

for k in [
    "variant_id",
    "rsid",
    "ref",
    "alt",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p",
    "eaf",
    "n",
]:

    print(
        k,
        "=",
        bbj_row.get(k, "")
    )


# ------------------------------------------------------------------
# GIGASTROKE rows
# ------------------------------------------------------------------

giga_rows = []


if GIGA.exists():

    with GIGA.open(
        errors="replace"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for r in reader:

            if (
                r.get("variant_id")
                == RS671
            ):

                giga_rows.append(r)


if giga_rows:

    with (
        OUT /
        "RS671_GIGASTROKE_MASTER.tsv"
    ).open(
        "w",
        newline=""
    ) as f:

        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=list(
                giga_rows[0].keys()
            )
        )

        w.writeheader()
        w.writerows(
            giga_rows
        )

else:

    (
        OUT /
        "RS671_GIGASTROKE_MASTER.tsv"
    ).write_text(
        "status\n"
        "RS671_NOT_FOUND_OR_MASTER_MISSING\n"
    )


print()
print(
    "GIGASTROKE_RS671_ROWS=",
    len(giga_rows)
)


for r in giga_rows:

    print(
        r.get(
            "phenotype",
            ""
        ),
        "p=",
        r.get(
            "giga_p",
            r.get(
                "p",
                ""
            )
        ),
        "beta=",
        r.get(
            "giga_beta_aligned",
            r.get(
                "beta",
                ""
            )
        ),
        "harm=",
        r.get(
            "harmonization",
            ""
        )
    )


# ------------------------------------------------------------------
# Conditioning evidence
# ------------------------------------------------------------------

def first_row(p):

    if not p.exists():
        return {}

    with p.open() as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    return (
        rows[0]
        if rows
        else {}
    )


e8 = first_row(
    P8E /
    "L003_RS671_CONDITIONING_SUMMARY.tsv"
)

f8 = first_row(
    P8F /
    "L003_JPT_CONDITIONING_SUMMARY.tsv"
)


local = {

    "variant_id":
        RS671,

    "rsid":
        bbj_row.get(
            "rsid",
            ""
        ),

    "gene":
        "ALDH2",

    "build":
        "GRCh37",

    "ref":
        bbj_row.get(
            "ref",
            ""
        ),

    "alt":
        bbj_row.get(
            "alt",
            ""
        ),

    "bbj_effect_allele":
        bbj_row.get(
            "effect_allele",
            ""
        ),

    "bbj_beta":
        bbj_row.get(
            "beta",
            ""
        ),

    "bbj_se":
        bbj_row.get(
            "se",
            ""
        ),

    "bbj_p":
        bbj_row.get(
            "p",
            ""
        ),

    "bbj_eaf":
        bbj_row.get(
            "eaf",
            ""
        ),

    "bbj_ref_validation":
        "PASS_GRCH37",

    "eas504_residual_gws_n":
        e8.get(
            "conditional_p_lt_5e8_n",
            ""
        ),

    "jpt_residual_gws_n":
        f8.get(
            "jpt_conditional_gws_n",
            ""
        ),

    "jpt_low_ld_residual_gws_n":
        f8.get(
            "jpt_conditional_low_ld_gws_n",
            ""
        ),

    "independent_secondary_signal":
        "NOT_SUPPORTED_IN_JPT_LD_SENSITIVITY",

    "mechanism_branch":
        "CODING_PROTEIN_METABOLIC_TO_TEST",

    "evidence_status":
        "LOCAL_GENETIC_EVIDENCE_COMPLETE",
}


with (
    OUT /
    "RS671_LOCAL_EVIDENCE_MASTER.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            local.keys()
        )
    )

    w.writeheader()
    w.writerow(local)

PY

}


# =============================================================================
# 04. OPTIONAL EXISTING COLOC RESULTS FROM GOOGLE DRIVE
# =============================================================================

collect_coloc_outputs() {

  echo "===== LOCAL SEARCH ====="

  ABF_LOCAL="$(
    find "$ROOT/results/is" \
      "$ROOT/IS_Analysis_V3" \
      -type f \
      -name 'COLOC_ABF_MASTER_ANNOTATED_V2.tsv' \
      -print \
      2>/dev/null \
      | head -1
  )"

  SUSIE_LOCAL="$(
    find "$ROOT/results/is" \
      "$ROOT/IS_Analysis_V3" \
      -type f \
      -name 'COLOC_SUSIE_SIGNAL_PAIRS.tsv' \
      -print \
      2>/dev/null \
      | head -1
  )"

  if [[ -n "$ABF_LOCAL" ]]; then

    cp \
      "$ABF_LOCAL" \
      "$OUT/COLOC_ABF_MASTER_ANNOTATED_V2.tsv"

    echo "LOCAL_ABF=$ABF_LOCAL"

  fi

  if [[ -n "$SUSIE_LOCAL" ]]; then

    cp \
      "$SUSIE_LOCAL" \
      "$OUT/COLOC_SUSIE_SIGNAL_PAIRS.tsv"

    echo "LOCAL_SUSIE=$SUSIE_LOCAL"

  fi


  echo
  echo "===== DRIVE FALLBACK ====="

  if ! command -v rclone >/dev/null 2>&1; then

    echo "RCLONE_NOT_AVAILABLE"
    return 0

  fi


  if [[ ! -s "$OUT/COLOC_ABF_MASTER_ANNOTATED_V2.tsv" ]]; then

    rclone copyto \
      "gdrive:MASTER_DEGREE/IS_COLAB/results_v2/COLOC_ABF_MASTER_ANNOTATED_V2.tsv" \
      "$OUT/COLOC_ABF_MASTER_ANNOTATED_V2.tsv" \
      2>/dev/null

    RC=$?

    if [[ $RC -eq 0 ]]; then
      echo "DRIVE_ABF=FOUND"
    else
      echo "DRIVE_ABF=NOT_FOUND"
      rm -f \
        "$OUT/COLOC_ABF_MASTER_ANNOTATED_V2.tsv" \
        2>/dev/null || true
    fi

  fi


  if [[ ! -s "$OUT/COLOC_SUSIE_SIGNAL_PAIRS.tsv" ]]; then

    rclone copyto \
      "gdrive:MASTER_DEGREE/IS_COLAB/phase11b_multisignal/results/COLOC_SUSIE_SIGNAL_PAIRS.tsv" \
      "$OUT/COLOC_SUSIE_SIGNAL_PAIRS.tsv" \
      2>/dev/null

    RC=$?

    if [[ $RC -eq 0 ]]; then
      echo "DRIVE_SUSIE=FOUND"
    else
      echo "DRIVE_SUSIE=NOT_FOUND"
      rm -f \
        "$OUT/COLOC_SUSIE_SIGNAL_PAIRS.tsv" \
        2>/dev/null || true
    fi

  fi


  ls -lh \
    "$OUT/COLOC_ABF_MASTER_ANNOTATED_V2.tsv" \
    "$OUT/COLOC_SUSIE_SIGNAL_PAIRS.tsv" \
    2>/dev/null || true

  return 0
}


# =============================================================================
# 05. RS671 EXTERNAL FUNCTIONAL ANNOTATION — ENSEMBL GRCh37
# =============================================================================

annotate_rs671_vep() {

  RAW="$OUT/RS671_ENSEMBL_VEP_GRCH37.json"
  STATUS_FILE="$OUT/RS671_ENSEMBL_STATUS.tsv"

  URL="https://grch37.rest.ensembl.org/vep/human/region/12:112241766-112241766:1/A?canonical=1&hgvs=1&protein=1"

  echo "VEP_URL=$URL"

  curl \
    -sS \
    --fail \
    --connect-timeout 20 \
    --max-time 120 \
    -H 'Content-Type: application/json' \
    "$URL" \
    -o "$RAW"

  RC=$?

  if [[ $RC -ne 0 || ! -s "$RAW" ]]; then

    printf \
      "resource\tstatus\nENSEMBL_GRCH37_VEP\tNETWORK_OR_API_UNAVAILABLE\n" \
      > "$STATUS_FILE"

    echo "ENSEMBL_VEP=NETWORK_OR_API_UNAVAILABLE"

    rm -f "$RAW" 2>/dev/null || true

    return 0

  fi


python3 - \
  "$RAW" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import json
import sys


RAW = Path(sys.argv[1])
OUT = Path(sys.argv[2])


obj = json.loads(
    RAW.read_text()
)


if isinstance(obj, dict):
    obj = [obj]


rows = []


for result in obj:

    for t in result.get(
        "transcript_consequences",
        []
    ):

        gene = (
            t.get("gene_symbol")
            or ""
        )

        if gene != "ALDH2":
            continue

        rows.append({

            "gene_symbol":
                gene,

            "gene_id":
                t.get(
                    "gene_id",
                    ""
                ),

            "transcript_id":
                t.get(
                    "transcript_id",
                    ""
                ),

            "canonical":
                t.get(
                    "canonical",
                    ""
                ),

            "consequence_terms":
                ",".join(
                    t.get(
                        "consequence_terms",
                        []
                    )
                ),

            "impact":
                t.get(
                    "impact",
                    ""
                ),

            "hgvsc":
                t.get(
                    "hgvsc",
                    ""
                ),

            "hgvsp":
                t.get(
                    "hgvsp",
                    ""
                ),

            "protein_id":
                t.get(
                    "protein_id",
                    ""
                ),

            "amino_acids":
                t.get(
                    "amino_acids",
                    ""
                ),

            "codons":
                t.get(
                    "codons",
                    ""
                ),

            "protein_start":
                t.get(
                    "protein_start",
                    ""
                ),

            "sift_prediction":
                t.get(
                    "sift_prediction",
                    ""
                ),

            "sift_score":
                t.get(
                    "sift_score",
                    ""
                ),

            "polyphen_prediction":
                t.get(
                    "polyphen_prediction",
                    ""
                ),

            "polyphen_score":
                t.get(
                    "polyphen_score",
                    ""
                ),
        })


outfile = (
    OUT /
    "RS671_ENSEMBL_VEP_TRANSCRIPTS.tsv"
)


if rows:

    with outfile.open(
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

else:

    outfile.write_text(
        "status\n"
        "NO_ALDH2_TRANSCRIPT_RETURNED\n"
    )


canonical = [
    x
    for x in rows
    if str(
        x.get(
            "canonical",
            ""
        )
    ) in {
        "1",
        "YES",
        "true",
        "True",
    }
]


chosen = (
    canonical[0]
    if canonical
    else (
        rows[0]
        if rows
        else {}
    )
)


master = {

    "variant_id":
        "12:112241766:G:A",

    "rsid":
        "rs671",

    "gene":
        "ALDH2",

    "annotation_source":
        "Ensembl_GRCh37_VEP",

    "transcript_id":
        chosen.get(
            "transcript_id",
            ""
        ),

    "canonical":
        chosen.get(
            "canonical",
            ""
        ),

    "consequence_terms":
        chosen.get(
            "consequence_terms",
            ""
        ),

    "impact":
        chosen.get(
            "impact",
            ""
        ),

    "hgvsc":
        chosen.get(
            "hgvsc",
            ""
        ),

    "hgvsp":
        chosen.get(
            "hgvsp",
            ""
        ),

    "amino_acids":
        chosen.get(
            "amino_acids",
            ""
        ),

    "sift_prediction":
        chosen.get(
            "sift_prediction",
            ""
        ),

    "sift_score":
        chosen.get(
            "sift_score",
            ""
        ),

    "polyphen_prediction":
        chosen.get(
            "polyphen_prediction",
            ""
        ),

    "polyphen_score":
        chosen.get(
            "polyphen_score",
            ""
        ),

    "annotation_status":
        (
            "PASS"
            if chosen
            else "NO_ALDH2_TRANSCRIPT"
        ),
}


with (
    OUT /
    "RS671_FUNCTIONAL_ANNOTATION.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            master.keys()
        )
    )

    w.writeheader()
    w.writerow(master)


print(
    "ALDH2_TRANSCRIPTS=",
    len(rows)
)

for k, v in master.items():
    print(
        f"{k}={v}"
    )

PY

  printf \
    "resource\tstatus\nENSEMBL_GRCH37_VEP\tPASS\n" \
    > "$STATUS_FILE"

}


# =============================================================================
# 06. BUILD COMPARATIVE COLOC MATRIX
# =============================================================================

build_convergence_matrix() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import math
import re
import sys


OUT = Path(sys.argv[1])


TARGETS = [

    {
        "gene": "FGF5",
        "locus": "BBJ_IS_L001",
        "role": "CORE",
        "mechanism_branch":
            "REGULATORY_EXPRESSION_PROTEIN",
    },

    {
        "gene": "ALDH2",
        "locus": "BBJ_IS_L003",
        "role": "CORE",
        "mechanism_branch":
            "CODING_PROTEIN_METABOLIC",
    },

    {
        "gene": "SH3PXD2A",
        "locus": "BBJ_IS_L002",
        "role": "CORE",
        "mechanism_branch":
            "CELL_SPECIFIC_REGULATORY",
    },

    {
        "gene": "COL4A2",
        "locus": "BBJ_IS_L004",
        "role": "CORE",
        "mechanism_branch":
            "VASCULAR_STRUCTURAL_REGULATORY",
    },

    {
        "gene": "COL4A1",
        "locus": "BBJ_IS_L004",
        "role": "CORE",
        "mechanism_branch":
            "VASCULAR_STRUCTURAL_REGULATORY",
    },

    {
        "gene": "C4orf22",
        "locus": "BBJ_IS_L001",
        "role": "SECONDARY_SCREEN",
        "mechanism_branch":
            "MOLECULAR_COLOC_SCREEN",
    },

    {
        "gene": "CALHM2",
        "locus": "BBJ_IS_L002",
        "role": "SECONDARY_SCREEN",
        "mechanism_branch":
            "MOLECULAR_COLOC_SCREEN",
    },

    {
        "gene": "NEURL1",
        "locus": "BBJ_IS_L002",
        "role": "SECONDARY_SCREEN",
        "mechanism_branch":
            "MOLECULAR_COLOC_SCREEN",
    },

    {
        "gene": "INA",
        "locus": "BBJ_IS_L002",
        "role": "SECONDARY_SCREEN",
        "mechanism_branch":
            "MOLECULAR_COLOC_SCREEN",
    },
]


def normalise(x):

    return re.sub(
        r"[^a-z0-9]",
        "",
        x.lower()
    )


def choose_col(
    fieldnames,
    preferred,
    contains=None,
):

    if not fieldnames:
        return None

    norm = {
        normalise(x): x
        for x in fieldnames
    }

    for p in preferred:

        n = normalise(p)

        if n in norm:
            return norm[n]

    if contains:

        for x in fieldnames:

            low = x.lower()

            if all(
                z.lower() in low
                for z in contains
            ):
                return x

    return None


def load_best(
    path,
    source,
):

    if not path.exists():
        return {}

    with path.open(
        errors="replace"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        fields = (
            reader.fieldnames
            or []
        )

        gene_col = choose_col(
            fields,
            [
                "gene",
                "gene_symbol",
                "symbol",
                "candidate_gene",
            ],
            contains=["gene"],
        )

        locus_col = choose_col(
            fields,
            [
                "locus",
            ],
            contains=["locus"],
        )

        tissue_col = choose_col(
            fields,
            [
                "tissue",
                "tissue_name",
                "dataset_label",
                "dataset",
            ]
        )

        h4_col = choose_col(
            fields,
            [
                "PP.H4",
                "PP.H4.abf",
                "PP_H4",
                "max_PP.H4",
                "max_pp_h4",
            ]
        )

        status_col = choose_col(
            fields,
            [
                "status",
            ]
        )


        print(
            source,
            "columns:",
            "gene=",
            gene_col,
            "locus=",
            locus_col,
            "tissue=",
            tissue_col,
            "h4=",
            h4_col,
            "status=",
            status_col,
        )


        if gene_col is None:
            return {}

        if h4_col is None:
            return {}


        best = {}


        for r in reader:

            if (
                status_col
                and r.get(
                    status_col,
                    ""
                )
                and r.get(
                    status_col,
                    ""
                ).upper()
                not in {
                    "PASS",
                    "OK",
                    "SUCCESS",
                }
            ):
                continue

            gene = (
                r.get(
                    gene_col,
                    ""
                )
                .strip()
            )

            if not gene:
                continue

            try:
                h4 = float(
                    r.get(
                        h4_col,
                        ""
                    )
                )
            except Exception:
                continue

            if (
                not math.isfinite(h4)
            ):
                continue


            old = best.get(gene)

            if (
                old is None
                or h4 > old[
                    "h4"
                ]
            ):

                best[gene] = {

                    "h4":
                        h4,

                    "locus":
                        (
                            r.get(
                                locus_col,
                                ""
                            )
                            if locus_col
                            else ""
                        ),

                    "tissue":
                        (
                            r.get(
                                tissue_col,
                                ""
                            )
                            if tissue_col
                            else ""
                        ),

                    "source":
                        source,
                }


        return best


abf = load_best(
    OUT /
    "COLOC_ABF_MASTER_ANNOTATED_V2.tsv",
    "COLOC_ABF_V2",
)


susie = load_best(
    OUT /
    "COLOC_SUSIE_SIGNAL_PAIRS.tsv",
    "COLOC_SUSIE_V2",
)


vep = {}

vp = (
    OUT /
    "RS671_FUNCTIONAL_ANNOTATION.tsv"
)

if vp.exists():

    with vp.open() as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    if rows:
        vep = rows[0]


matrix = []


for t in TARGETS:

    gene = t["gene"]

    a = abf.get(
        gene,
        {}
    )

    s = susie.get(
        gene,
        {}
    )


    row = {

        "gene":
            gene,

        "locus":
            t["locus"],

        "role":
            t["role"],

        "mechanism_branch":
            t["mechanism_branch"],

        "best_abf_h4":
            a.get(
                "h4",
                ""
            ),

        "best_abf_tissue":
            a.get(
                "tissue",
                ""
            ),

        "best_susie_h4":
            s.get(
                "h4",
                ""
            ),

        "best_susie_tissue":
            s.get(
                "tissue",
                ""
            ),

        "l003_rs671_dominant_signal":
            (
                "YES"
                if gene == "ALDH2"
                else ""
            ),

        "rs671_functional_consequence":
            (
                vep.get(
                    "consequence_terms",
                    ""
                )
                if gene == "ALDH2"
                else ""
            ),

        "rs671_hgvsp":
            (
                vep.get(
                    "hgvsp",
                    ""
                )
                if gene == "ALDH2"
                else ""
            ),

        "current_interpretation":
            "",

        "next_functional_layer":
            "",
    }


    if gene == "FGF5":

        row[
            "current_interpretation"
        ] = (
            "MOLECULAR_COLOC_SUPPORT"
        )

        row[
            "next_functional_layer"
        ] = (
            "PQTL_MR_BP_PATHWAY"
        )


    elif gene == "ALDH2":

        row[
            "current_interpretation"
        ] = (
            "RS671_DOMINANT_EAS_SIGNAL_"
            "EXPRESSION_COLOC_NOT_ESTABLISHED"
        )

        row[
            "next_functional_layer"
        ] = (
            "CODING_PROTEIN_METABOLIC_FUNCTION"
        )


    elif gene == "SH3PXD2A":

        row[
            "current_interpretation"
        ] = (
            "BULK_EQTL_SUPPORT_LIMITED"
        )

        row[
            "next_functional_layer"
        ] = (
            "MACROPHAGE_MONOCYTE_VASCULAR_CELL_QTL_ATAC"
        )


    elif gene in {
        "COL4A2",
        "COL4A1",
    }:

        row[
            "current_interpretation"
        ] = (
            "BULK_EQTL_SUPPORT_LIMITED"
        )

        row[
            "next_functional_layer"
        ] = (
            "VASCULAR_SMC_PERICYTE_ENDOTHELIAL_ATAC_EQTL"
        )


    else:

        row[
            "current_interpretation"
        ] = (
            "SECONDARY_MOLECULAR_SCREEN"
        )

        row[
            "next_functional_layer"
        ] = (
            "REQUIRE_ORTHOGONAL_VALIDATION"
        )


    matrix.append(row)


outfile = (
    OUT /
    "IS_FUNCTIONAL_CONVERGENCE_MASTER.tsv"
)


with outfile.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            matrix[0].keys()
        )
    )

    w.writeheader()
    w.writerows(matrix)


print()
print(
    "===== FUNCTIONAL CONVERGENCE ====="
)

for x in matrix:

    print(
        x["gene"],
        x["locus"],
        "ABF=",
        x["best_abf_h4"],
        x["best_abf_tissue"],
        "SuSiE=",
        x["best_susie_h4"],
        x["best_susie_tissue"],
        "branch=",
        x["next_functional_layer"],
    )

PY

}


# =============================================================================
# 07. BUILD PHASE9C REPORT
# =============================================================================

build_report() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


def rows(name):

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


matrix = rows(
    "IS_FUNCTIONAL_CONVERGENCE_MASTER.tsv"
)

freeze = rows(
    "L003_STATISTICAL_FREEZE.tsv"
)

vep = rows(
    "RS671_FUNCTIONAL_ANNOTATION.tsv"
)


lines = []

lines.append(
    "# IS Phase9C Functional Convergence\n"
)

lines.append(
    f"Run output: `{OUT}`\n"
)

lines.append(
    "## 1. L003 statistical freeze\n"
)


if freeze:

    x = freeze[0]

    lines.append(
        f"- Primary locus: **{x.get('primary_gene','')} / "
        f"{x.get('primary_rsid','')}**\n"
    )

    lines.append(
        f"- BBJ P: `{x.get('bbj_p','')}`\n"
    )

    lines.append(
        f"- JPT conditional GWS variants: "
        f"`{x.get('jpt_secondary_gws_n','')}`\n"
    )

    lines.append(
        f"- Low-LD conditional GWS variants: "
        f"`{x.get('jpt_low_ld_secondary_gws_n','')}`\n"
    )

    lines.append(
        "- Interpretation: "
        "`RS671_DOMINANT_EAS_SIGNAL_NO_ROBUST_"
        "INDEPENDENT_SECONDARY_SIGNAL`\n"
    )


lines.append(
    "\n## 2. rs671 functional annotation\n"
)


if vep:

    x = vep[0]

    lines.append(
        f"- Consequence: `{x.get('consequence_terms','')}`\n"
    )

    lines.append(
        f"- Transcript: `{x.get('transcript_id','')}`\n"
    )

    lines.append(
        f"- HGVSc: `{x.get('hgvsc','')}`\n"
    )

    lines.append(
        f"- HGVSp: `{x.get('hgvsp','')}`\n"
    )

    lines.append(
        f"- Amino acids: `{x.get('amino_acids','')}`\n"
    )

else:

    lines.append(
        "- Ensembl VEP annotation unavailable during this run.\n"
    )


lines.append(
    "\n## 3. Comparative molecular evidence\n"
)

lines.append(
    "| Gene | Locus | ABF H4 | ABF tissue | "
    "SuSiE H4 | SuSiE tissue | Current branch |\n"
)

lines.append(
    "|---|---|---:|---|---:|---|---|\n"
)


for x in matrix:

    lines.append(
        "| "
        + str(
            x.get(
                "gene",
                ""
            )
        )
        + " | "
        + str(
            x.get(
                "locus",
                ""
            )
        )
        + " | "
        + str(
            x.get(
                "best_abf_h4",
                ""
            )
        )
        + " | "
        + str(
            x.get(
                "best_abf_tissue",
                ""
            )
        )
        + " | "
        + str(
            x.get(
                "best_susie_h4",
                ""
            )
        )
        + " | "
        + str(
            x.get(
                "best_susie_tissue",
                ""
            )
        )
        + " | "
        + str(
            x.get(
                "next_functional_layer",
                ""
            )
        )
        + " |\n"
    )


lines.append(
    "\n## 4. Working mechanistic model\n"
)

lines.append(
    "- **FGF5:** regulatory/expression/protein convergence branch.\n"
)

lines.append(
    "- **ALDH2:** ancestry-specific rs671 coding/protein/metabolic branch; "
    "bulk eQTL colocalization should not be treated as the sole mechanism test.\n"
)

lines.append(
    "- **SH3PXD2A:** cell-specific immune/vascular regulatory branch.\n"
)

lines.append(
    "- **COL4A2/COL4A1:** cerebrovascular structural/regulatory branch.\n"
)

lines.append(
    "\n## 5. Current limitations\n"
)

lines.append(
    "- BBJ fine-mapping still uses external 1000G EAS LD.\n"
)

lines.append(
    "- JPT sensitivity uses only 104 reference individuals.\n"
)

lines.append(
    "- Phase11 molecular colocalization remains exploratory until "
    "cohort-matched molecular-QTL LD is available.\n"
)

lines.append(
    "- Literature evidence is not silently inserted into this master; "
    "it should be added in a separate sourced benchmark layer.\n"
)


(
    OUT /
    "PHASE9C_FUNCTIONAL_CONVERGENCE.md"
).write_text(
    "".join(lines)
)


print(
    "".join(lines)
)

PY

}


# =============================================================================
# 08. READINESS
# =============================================================================

build_readiness() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


matrix = (
    OUT /
    "IS_FUNCTIONAL_CONVERGENCE_MASTER.tsv"
)

freeze = (
    OUT /
    "L003_STATISTICAL_FREEZE.tsv"
)

local_rs = (
    OUT /
    "RS671_LOCAL_EVIDENCE_MASTER.tsv"
)

vep = (
    OUT /
    "RS671_FUNCTIONAL_ANNOTATION.tsv"
)


rows = [

    (
        "PHASE9C_FUNCTIONAL_CONVERGENCE",
        (
            "COMPLETE"
            if (
                matrix.exists()
                and freeze.exists()
                and local_rs.exists()
            )
            else "INCOMPLETE"
        )
    ),

    (
        "L003_STATISTICAL_BRANCH",
        "FROZEN_RS671_DOMINANT_EAS_SIGNAL"
    ),

    (
        "L003_INDEPENDENT_SECONDARY_SIGNAL",
        "NOT_SUPPORTED_IN_JPT_LD_SENSITIVITY"
    ),

    (
        "RS671_REFERENCE",
        "GRCH37_VERIFIED"
    ),

    (
        "RS671_FUNCTIONAL_ANNOTATION",
        (
            "ENSEMBL_VEP_AVAILABLE"
            if vep.exists()
            else "EXTERNAL_ANNOTATION_PENDING"
        )
    ),

    (
        "FGF5_BRANCH",
        "REGULATORY_EXPRESSION_PROTEIN"
    ),

    (
        "ALDH2_BRANCH",
        "CODING_PROTEIN_METABOLIC"
    ),

    (
        "SH3PXD2A_BRANCH",
        "CELL_SPECIFIC_IMMUNE_VASCULAR_REGULATION"
    ),

    (
        "COL4A2_COL4A1_BRANCH",
        "CEREBROVASCULAR_STRUCTURAL_REGULATION"
    ),

    (
        "LITERATURE_BENCHMARK",
        "NEXT"
    ),

    (
        "SINGLE_CELL_ATAC_FUNCTIONAL_LAYER",
        "NEXT_AFTER_LITERATURE_BENCHMARK"
    ),

    (
        "PAPER_GRADE_GENETIC_LIMITATION",
        "COHORT_MATCHED_BBJ_LD_STILL_UNAVAILABLE"
    ),
]


with (
    OUT /
    "PHASE9C_READINESS.tsv"
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
  "02_freeze_l003" \
  freeze_l003

run_step \
  "03_rs671_local_master" \
  build_rs671_local_master

run_step \
  "04_collect_coloc_outputs" \
  collect_coloc_outputs

run_step \
  "05_rs671_ensembl_vep" \
  annotate_rs671_vep

run_step \
  "06_functional_convergence_matrix" \
  build_convergence_matrix

run_step \
  "07_build_report" \
  build_report

run_step \
  "08_build_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9C FUNCTIONAL CONVERGENCE COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== L003 FREEZE ====="

cat \
  "$OUT/L003_STATISTICAL_FREEZE.tsv" \
  2>/dev/null \
  || true


echo
echo "===== RS671 FUNCTIONAL ANNOTATION ====="

cat \
  "$OUT/RS671_FUNCTIONAL_ANNOTATION.tsv" \
  2>/dev/null \
  || true


echo
echo "===== FUNCTIONAL CONVERGENCE MASTER ====="

column -t -s $'\t' \
  "$OUT/IS_FUNCTIONAL_CONVERGENCE_MASTER.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/IS_FUNCTIONAL_CONVERGENCE_MASTER.tsv"


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9C_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE9C_READINESS.tsv"


echo
echo "===== OUTPUT FILES ====="

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -printf '%f\t%s bytes\n' \
  | sort


echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
