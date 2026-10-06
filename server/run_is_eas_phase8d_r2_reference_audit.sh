#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

R1="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit_r1"

OUT="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit_r2"
REFDIR="$ROOT/reference/grch37_1000g"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_eas_phase8d_r2/$RUN_ID"
STATUS="$OUT/STEP_STATUS_R2.tsv"

mkdir -p "$OUT" "$REFDIR" "$LOGDIR"

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

  echo "===== R1 READINESS ====="

  if [[ ! -s "$R1/SCIENTIFIC_READINESS_R1.tsv" ]]; then
    echo "R1_READINESS_MISSING"
    return 1
  fi

  cat "$R1/SCIENTIFIC_READINESS_R1.tsv"

  echo
  echo "===== R1 STEP STATUS ====="

  if [[ ! -s "$R1/STEP_STATUS_R1.tsv" ]]; then
    echo "R1_STEP_STATUS_MISSING"
    return 1
  fi

  cat "$R1/STEP_STATUS_R1.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$R1/STEP_STATUS_R1.tsv"
  )"

  echo
  echo "R1_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    echo "R1_NOT_CLEAN"
    return 1
  fi

  echo
  echo "===== DEPENDENCIES ====="

  for exe in \
    python3 \
    curl \
    gzip \
    samtools \
    sha256sum
  do

    if command -v "$exe" >/dev/null 2>&1; then
      echo "FOUND $exe $(command -v "$exe")"
    else
      echo "MISSING $exe"
      return 1
    fi

  done

  echo
  echo "===== DISK ====="

  df -h "$ROOT"

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. GRCh37 / 1000G REFERENCE
# =============================================================================

prepare_reference() {

  REF_GZ="$REFDIR/human_g1k_v37.fasta.gz"
  REF="$REFDIR/human_g1k_v37.fasta"
  FAI="$REF.fai"

  URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

  echo "REFERENCE_SOURCE_URL=$URL"

  if [[ ! -s "$REF" ]]; then

    if [[ ! -s "$REF_GZ" ]]; then

      echo "===== DOWNLOAD GRCh37 REFERENCE ====="

      curl \
        -L \
        --fail \
        --retry 3 \
        --retry-delay 5 \
        -o "$REF_GZ.part" \
        "$URL"

      RC=$?

      if [[ $RC -ne 0 ]]; then
        echo "REFERENCE_DOWNLOAD_FAIL=$RC"
        return 1
      fi

      mv "$REF_GZ.part" "$REF_GZ"

    else

      echo "REFERENCE_GZ_ALREADY_EXISTS"

    fi

    echo
    echo "===== DECOMPRESS ====="

    gzip -dc "$REF_GZ" > "$REF.part"

    RC=$?

    if [[ $RC -ne 0 || ! -s "$REF.part" ]]; then
      echo "REFERENCE_DECOMPRESS_FAIL=$RC"
      return 1
    fi

    mv "$REF.part" "$REF"

  else

    echo "REFERENCE_FASTA_ALREADY_EXISTS"

  fi

  echo
  echo "===== SHA256 ====="

  sha256sum "$REF" \
    | tee "$OUT/GRCH37_REFERENCE_SHA256.txt"

  echo
  echo "===== FAI ====="

  if [[ ! -s "$FAI" ]]; then
    samtools faidx "$REF"
  fi

  if [[ ! -s "$FAI" ]]; then
    echo "FAI_CREATION_FAIL"
    return 1
  fi

  head -30 "$FAI"

  echo
  echo "===== CONTIG LENGTH CHECK ====="

  python3 - "$FAI" "$OUT" <<'PY'
from pathlib import Path
import csv
import sys

fai = Path(sys.argv[1])
out = Path(sys.argv[2])

expected = {
    "1": 249250621,
    "4": 191154276,
    "10": 135534747,
    "12": 133851895,
    "13": 115169878,
    "19": 59128983,
}

got = {}

with fai.open() as f:
    for line in f:
        x = line.rstrip("\n").split("\t")
        if len(x) < 2:
            continue
        got[x[0]] = int(x[1])

rows = []

for chrom, exp in expected.items():

    candidates = [
        chrom,
        "chr" + chrom,
    ]

    found = next(
        (
            c for c in candidates
            if c in got
        ),
        None
    )

    actual = (
        got.get(found)
        if found
        else None
    )

    status = (
        "PASS"
        if actual == exp
        else "FAIL"
    )

    rows.append({
        "chromosome": chrom,
        "resolved_contig": found or "",
        "expected_grch37_length": exp,
        "actual_length": actual if actual is not None else "",
        "status": status,
    })

with (
    out /
    "GRCH37_REFERENCE_CONTIG_QC.tsv"
).open("w", newline="") as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            rows[0].keys()
        ),
    )

    w.writeheader()
    w.writerows(rows)

for r in rows:
    print(
        r["chromosome"],
        r["resolved_contig"],
        r["expected_grch37_length"],
        r["actual_length"],
        r["status"]
    )

if any(
    r["status"] != "PASS"
    for r in rows
):
    raise SystemExit(
        "GRCh37 contig fingerprint failed"
    )

print("GRCH37_CONTIG_FINGERPRINT=PASS")
PY

}


# =============================================================================
# 03. BBJ PRIMARY-LOCUS REF AUDIT
# =============================================================================

audit_bbj_ref() {

python3 - \
  "$REFDIR/human_g1k_v37.fasta" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import subprocess
import sys
import re


FASTA = Path(sys.argv[1])
OUT = Path(sys.argv[2])

BBJ = Path(
    "/srv/is-analysis/data/is/processed/japan/bbj/"
    "BBJ_IS_GRCh37.canonical.tsv.gz"
)

if not BBJ.exists():
    raise SystemExit(
        f"BBJ canonical missing: {BBJ}"
    )


LOCI = {
    "BBJ_IS_L001": ("4", 80681072, 81684341),
    "BBJ_IS_L002": ("10", 104839048, 105960401),
    "BBJ_IS_L003": ("12", 110823939, 113617897),
    "BBJ_IS_L004": ("13", 110454237, 111459643),
}


# ------------------------------------------------------------------
# FAI → resolve 1 vs chr1
# ------------------------------------------------------------------

contigs = {}

with Path(
    str(FASTA) + ".fai"
).open() as f:

    for line in f:

        x = line.rstrip("\n").split("\t")

        if len(x) >= 2:
            contigs[x[0]] = int(x[1])


def resolve_contig(chrom):

    if chrom in contigs:
        return chrom

    x = "chr" + chrom

    if x in contigs:
        return x

    raise RuntimeError(
        f"contig unavailable: {chrom}"
    )


# ------------------------------------------------------------------
# Fetch four whole sequences once
# ------------------------------------------------------------------

seqs = {}

for locus, (
    chrom,
    start,
    end,
) in LOCI.items():

    contig = resolve_contig(chrom)

    region = (
        f"{contig}:{start}-{end}"
    )

    p = subprocess.run(
        [
            "samtools",
            "faidx",
            str(FASTA),
            region,
        ],
        capture_output=True,
        text=True,
    )

    if p.returncode != 0:
        raise RuntimeError(
            p.stderr
        )

    lines = p.stdout.splitlines()

    seq = "".join(
        x.strip()
        for x in lines[1:]
    ).upper()

    expected_len = end - start + 1

    if len(seq) != expected_len:
        raise RuntimeError(
            f"{locus}: sequence length "
            f"{len(seq)} != {expected_len}"
        )

    seqs[locus] = seq


# ------------------------------------------------------------------
# Header detection
# ------------------------------------------------------------------

with gzip.open(
    BBJ,
    "rt",
    errors="replace",
) as f:

    reader = csv.reader(
        f,
        delimiter="\t",
    )

    header = next(reader)


lower = {
    x.lower(): x
    for x in header
}


def choose(names):

    for n in names:

        if n in lower:
            return lower[n]

    return None


chr_col = choose([
    "chr",
    "chrom",
    "chromosome",
])

pos_col = choose([
    "pos",
    "position",
    "base_pair_location",
    "bp",
])

ref_col = choose([
    "ref",
    "reference",
    "reference_allele",
    "bbj_ref",
])

alt_col = choose([
    "alt",
    "alternate",
    "alternate_allele",
    "bbj_alt",
])

print("BBJ_HEADER=", header)
print(
    "RESOLVED=",
    chr_col,
    pos_col,
    ref_col,
    alt_col,
)

if None in [
    chr_col,
    pos_col,
    ref_col,
    alt_col,
]:

    raise RuntimeError(
        "Cannot resolve BBJ canonical columns"
    )


# ------------------------------------------------------------------
# Streaming audit
# ------------------------------------------------------------------

stats = {
    locus: {
        "rows": 0,
        "ref_match": 0,
        "ref_mismatch": 0,
        "non_acgt_ref": 0,
        "out_of_bounds": 0,
    }
    for locus in LOCI
}

mismatches = []


def normalize_chr(x):

    x = str(x).strip()

    if x.lower().startswith("chr"):
        x = x[3:]

    return x


with gzip.open(
    BBJ,
    "rt",
    errors="replace",
) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t",
    )

    for r in reader:

        chrom = normalize_chr(
            r[chr_col]
        )

        try:
            pos = int(
                float(
                    r[pos_col]
                )
            )
        except Exception:
            continue

        locus = None

        for loc, (
            c,
            lo,
            hi,
        ) in LOCI.items():

            if (
                chrom == c
                and lo <= pos <= hi
            ):
                locus = loc
                break

        if locus is None:
            continue

        ref = (
            str(
                r[ref_col]
            )
            .strip()
            .upper()
        )

        alt = (
            str(
                r[alt_col]
            )
            .strip()
            .upper()
        )

        st = stats[locus]
        st["rows"] += 1

        c, lo, hi = LOCI[locus]

        offset = pos - lo

        if (
            offset < 0
            or (
                offset
                + len(ref)
                > len(
                    seqs[locus]
                )
            )
        ):

            st["out_of_bounds"] += 1
            continue

        if not re.fullmatch(
            r"[ACGT]+",
            ref,
        ):

            st["non_acgt_ref"] += 1
            continue

        reference = seqs[locus][
            offset:
            offset + len(ref)
        ]

        if ref == reference:

            st["ref_match"] += 1

        else:

            st["ref_mismatch"] += 1

            if len(mismatches) < 5000:

                mismatches.append({
                    "locus": locus,
                    "chr": chrom,
                    "pos": pos,
                    "bbj_ref": ref,
                    "bbj_alt": alt,
                    "reference": reference,
                })


rows = []

for locus, s in stats.items():

    informative = (
        s["ref_match"]
        + s["ref_mismatch"]
    )

    rate = (
        s["ref_match"] / informative
        if informative
        else float("nan")
    )

    status = (
        "PASS"
        if (
            informative > 0
            and s["ref_mismatch"] == 0
        )
        else "REVIEW"
    )

    rows.append({
        "locus": locus,
        **s,
        "informative_ref_rows": informative,
        "ref_match_rate": rate,
        "status": status,
    })


with (
    OUT /
    "BBJ_GRCH37_REF_AUDIT.tsv"
).open("w", newline="") as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            rows[0].keys()
        ),
    )

    w.writeheader()
    w.writerows(rows)


with (
    OUT /
    "BBJ_GRCH37_REF_MISMATCHES.tsv"
).open("w", newline="") as f:

    fields = [
        "locus",
        "chr",
        "pos",
        "bbj_ref",
        "bbj_alt",
        "reference",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(
        mismatches
    )


for r in rows:
    print(r)


if any(
    x["status"] != "PASS"
    for x in rows
):
    print(
        "BBJ_GRCH37_REF_AUDIT=REVIEW"
    )
else:
    print(
        "BBJ_GRCH37_REF_AUDIT=PASS"
    )

PY

}


# =============================================================================
# 04. GIGASTROKE RAW ALLELE → GRCh37 REFERENCE AUDIT
# =============================================================================

audit_gigastroke_ref() {

python3 - \
  "$REFDIR/human_g1k_v37.fasta" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import gzip
import subprocess
import sys
import re


FASTA = Path(sys.argv[1])
OUT = Path(sys.argv[2])

BASE = Path(
    "/srv/is-analysis/data/is/reference/"
    "gigastroke/eas"
)


FILES = {
    "AS": BASE /
        "GCST90104544_buildGRCh37.tsv.gz",

    "AIS": BASE /
        "GCST90104545_buildGRCh37.tsv.gz",

    "CES": BASE /
        "GCST90104546_buildGRCh37.tsv.gz",

    "LAS": BASE /
        "GCST90104547_buildGRCh37.tsv.gz",

    "SVS": BASE /
        "GCST90104548_buildGRCh37.tsv.gz",
}


LOCI = {
    "BBJ_IS_L001": ("4", 80681072, 81684341),
    "BBJ_IS_L002": ("10", 104839048, 105960401),
    "BBJ_IS_L003": ("12", 110823939, 113617897),
    "BBJ_IS_L004": ("13", 110454237, 111459643),
}


for phenotype, p in FILES.items():

    if not p.exists():
        raise RuntimeError(
            f"missing {phenotype}: {p}"
        )


# ------------------------------------------------------------------
# Resolve contigs
# ------------------------------------------------------------------

contigs = {}

with Path(
    str(FASTA) + ".fai"
).open() as f:

    for line in f:

        x = line.rstrip("\n").split("\t")

        if len(x) >= 2:
            contigs[x[0]] = int(x[1])


def resolve_contig(chrom):

    if chrom in contigs:
        return chrom

    c = "chr" + chrom

    if c in contigs:
        return c

    raise RuntimeError(chrom)


seqs = {}

for locus, (
    chrom,
    start,
    end,
) in LOCI.items():

    contig = resolve_contig(
        chrom
    )

    region = (
        f"{contig}:{start}-{end}"
    )

    p = subprocess.run(
        [
            "samtools",
            "faidx",
            str(FASTA),
            region,
        ],
        capture_output=True,
        text=True,
    )

    if p.returncode != 0:
        raise RuntimeError(
            p.stderr
        )

    seq = "".join(
        x.strip()
        for x in p.stdout.splitlines()[1:]
    ).upper()

    seqs[locus] = seq


def nchrom(x):

    x = str(x).strip()

    if x.lower().startswith(
        "chr"
    ):
        x = x[3:]

    return x


summary = []
bad = []


for phenotype, path in FILES.items():

    print(
        "AUDIT",
        phenotype,
        path,
    )

    stats = {
        locus: {
            "rows": 0,
            "snv_rows": 0,
            "effect_is_ref": 0,
            "other_is_ref": 0,
            "neither_is_ref": 0,
            "both_is_ref": 0,
            "complex_or_non_snv": 0,
        }
        for locus in LOCI
    }

    with gzip.open(
        path,
        "rt",
        errors="replace",
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t",
        )

        required = [
            "chromosome",
            "base_pair_location",
            "effect_allele",
            "other_allele",
        ]

        missing = [
            x
            for x in required
            if x not in reader.fieldnames
        ]

        if missing:
            raise RuntimeError(
                f"{phenotype}: "
                f"missing {missing}"
            )

        for r in reader:

            chrom = nchrom(
                r["chromosome"]
            )

            try:
                pos = int(
                    float(
                        r[
                            "base_pair_location"
                        ]
                    )
                )
            except Exception:
                continue

            locus = None

            for loc, (
                c,
                lo,
                hi,
            ) in LOCI.items():

                if (
                    chrom == c
                    and lo <= pos <= hi
                ):
                    locus = loc
                    break

            if locus is None:
                continue

            s = stats[locus]
            s["rows"] += 1

            ea = (
                str(
                    r["effect_allele"]
                )
                .strip()
                .upper()
            )

            oa = (
                str(
                    r["other_allele"]
                )
                .strip()
                .upper()
            )

            if not (
                len(ea) == 1
                and len(oa) == 1
                and ea in "ACGT"
                and oa in "ACGT"
            ):

                s[
                    "complex_or_non_snv"
                ] += 1
                continue

            s["snv_rows"] += 1

            c, lo, hi = LOCI[locus]

            ref = seqs[locus][
                pos - lo
            ]

            e_ref = (
                ea == ref
            )

            o_ref = (
                oa == ref
            )

            if (
                e_ref
                and not o_ref
            ):

                s[
                    "effect_is_ref"
                ] += 1

            elif (
                o_ref
                and not e_ref
            ):

                s[
                    "other_is_ref"
                ] += 1

            elif (
                e_ref
                and o_ref
            ):

                s[
                    "both_is_ref"
                ] += 1

            else:

                s[
                    "neither_is_ref"
                ] += 1

                if len(bad) < 10000:

                    bad.append({
                        "phenotype": phenotype,
                        "locus": locus,
                        "chr": chrom,
                        "pos": pos,
                        "effect_allele": ea,
                        "other_allele": oa,
                        "reference_base": ref,
                    })


    for locus, s in stats.items():

        resolved = (
            s["effect_is_ref"]
            + s["other_is_ref"]
        )

        rate = (
            resolved /
            s["snv_rows"]
            if s["snv_rows"]
            else float("nan")
        )

        status = (
            "PASS"
            if (
                s["snv_rows"] > 0
                and s["neither_is_ref"] == 0
                and s["both_is_ref"] == 0
            )
            else "REVIEW"
        )

        summary.append({
            "phenotype": phenotype,
            "locus": locus,
            **s,
            "reference_resolved_snv": resolved,
            "reference_resolved_rate": rate,
            "status": status,
        })


with (
    OUT /
    "GIGASTROKE_GRCH37_REF_AUDIT.tsv"
).open("w", newline="") as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            summary[0].keys()
        ),
    )

    w.writeheader()
    w.writerows(summary)


with (
    OUT /
    "GIGASTROKE_GRCH37_REF_UNRESOLVED.tsv"
).open("w", newline="") as f:

    fields = [
        "phenotype",
        "locus",
        "chr",
        "pos",
        "effect_allele",
        "other_allele",
        "reference_base",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(bad)


for x in summary:

    print(
        x["phenotype"],
        x["locus"],
        "rows=",
        x["rows"],
        "snv=",
        x["snv_rows"],
        "EA_REF=",
        x["effect_is_ref"],
        "OA_REF=",
        x["other_is_ref"],
        "NEITHER=",
        x["neither_is_ref"],
        "COMPLEX=",
        x["complex_or_non_snv"],
        "rate=",
        x["reference_resolved_rate"],
        x["status"],
    )

PY

}


# =============================================================================
# 05. READINESS R2
# =============================================================================

build_readiness() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


def read_tsv(name):

    p = OUT / name

    if not p.exists():
        return []

    with p.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


ref_qc = read_tsv(
    "GRCH37_REFERENCE_CONTIG_QC.tsv"
)

bbj = read_tsv(
    "BBJ_GRCH37_REF_AUDIT.tsv"
)

giga = read_tsv(
    "GIGASTROKE_GRCH37_REF_AUDIT.tsv"
)


reference_ok = (
    len(ref_qc) > 0
    and all(
        x["status"] == "PASS"
        for x in ref_qc
    )
)


bbj_ok = (
    len(bbj) == 4
    and all(
        x["status"] == "PASS"
        for x in bbj
    )
)


giga_ok = (
    len(giga) == 20
    and all(
        x["status"] == "PASS"
        for x in giga
    )
)


rows = [

    (
        "PHASE8D_R2",
        "AUDIT_COMPLETE"
    ),

    (
        "GRCH37_REFERENCE",
        (
            "PASS_1000G_GRCH37_FINGERPRINT"
            if reference_ok
            else "REVIEW"
        ),
    ),

    (
        "BBJ_REF_FASTA",
        (
            "VERIFIED_PRIMARY_LOCI_GRCH37"
            if bbj_ok
            else "REVIEW_REQUIRED"
        ),
    ),

    (
        "GIGASTROKE_REF_ALT",
        (
            "REFERENCE_RESOLVED_PRIMARY_LOCI_BIALLELIC_SNV"
            if giga_ok
            else "PARTIAL_OR_REVIEW_REQUIRED"
        ),
    ),

    (
        "GIGASTROKE_ANCESTRY",
        "VERIFIED_LOCAL_METADATA_EAST_ASIAN"
    ),

    (
        "GIGASTROKE_SAMPLE_SIZE",
        "VERIFIED_LOCAL_METADATA_TOTAL_N"
    ),

    (
        "HARMONIZATION",
        "FINAL_AUDITED_FILES_MATCH_SWAP_ONLY_NO_COMPLEMENT"
    ),

    (
        "LD_ORIENTATION",
        "DOCUMENTED_SIGNED_REF_BASED_FROM_PHASE5_PROVENANCE"
    ),

    (
        "LD_STRUCTURAL_QC",
        "PASS_ALL_DISCOVERED_MATRICES"
    ),

    (
        "SUSIER_RSS_API",
        "AUDITED"
    ),

    (
        "EXACT_SUBSET_MATRIX_PROVENANCE",
        "NOT_REUSED_RECONSTRUCT_FRESH_FOR_NEW_ANALYSES"
    ),

    (
        "BBJ_PRIMARY_FINEMAP",
        "EXPLORATORY_EAS504"
    ),

    (
        "PHASE9_11_OUTPUTS",
        "EXPLORATORY_FROZEN"
    ),

    (
        "PHASE11B_FGF5",
        "EXPLORATORY_MOLECULAR_COLOC_SUPPORT"
    ),

    (
        "L003_CONDITIONING",
        (
            "READY_FOR_FRESH_EAS_RS671_ANALYSIS"
            if (
                reference_ok
                and bbj_ok
            )
            else
            "PENDING_REFERENCE_VALIDATION"
        ),
    ),

    (
        "TPMI",
        "WAITING_RATE_LIMIT"
    ),

    (
        "CKB",
        "WAITING_DECRYPTION_KEY"
    ),
]


p = (
    OUT /
    "SCIENTIFIC_READINESS_R2.tsv"
)

with p.open(
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
    )

    w.writerow([
        "component",
        "status",
    ])

    w.writerows(rows)


for a, b in rows:
    print(
        f"{a} = {b}"
    )


if (
    reference_ok
    and bbj_ok
):
    print(
        "L003_NEXT_GATE=OPEN"
    )
else:
    print(
        "L003_NEXT_GATE=CLOSED"
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
  "02_prepare_grch37_reference" \
  prepare_reference

run_step \
  "03_bbj_grch37_ref_audit" \
  audit_bbj_ref

run_step \
  "04_gigastroke_grch37_ref_audit" \
  audit_gigastroke_ref

run_step \
  "05_build_readiness_r2" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE8D-R2 COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="
column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"

echo
echo "===== REFERENCE CONTIG QC ====="
cat \
  "$OUT/GRCH37_REFERENCE_CONTIG_QC.tsv" \
  2>/dev/null \
  || true

echo
echo "===== BBJ REF AUDIT ====="
cat \
  "$OUT/BBJ_GRCH37_REF_AUDIT.tsv" \
  2>/dev/null \
  || true

echo
echo "===== GIGASTROKE REF AUDIT ====="
cat \
  "$OUT/GIGASTROKE_GRCH37_REF_AUDIT.tsv" \
  2>/dev/null \
  || true

echo
echo "===== SCIENTIFIC READINESS R2 ====="
column -t -s $'\t' \
  "$OUT/SCIENTIFIC_READINESS_R2.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/SCIENTIFIC_READINESS_R2.tsv"

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

