#!/usr/bin/env bash

# ============================================================
# IS PHASE 8D — CORE SCIENTIFIC VALIDATION AUDIT
#
# READ-ONLY scientific audit.
#
# DOES NOT:
#   - delete/move existing outputs
#   - download external data
#   - modify LD matrices
#   - symmetrise LD
#   - run SuSiE / coloc
#   - perform complement rescue
#   - advance Phase9/10/11
#
# ============================================================

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

OUT="$ROOT/results/is/stage4_cross_eas/phase8d_validation_audit"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_eas_phase8d/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

mkdir -p "$OUT" "$LOGDIR"

printf "step\tstatus\trc\n" > "$STATUS"

run_step() {
  local name="$1"
  shift
  local LOG="$LOGDIR/${name}.log"

  echo
  echo "================================================================================"
  echo "${name}"
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


# ============================================================
# 1. GIGASTROKE LOCAL METADATA AUDIT
# ============================================================

metadata_audit() {

python3 - <<'PY'
from pathlib import Path
import json
import re
import csv

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

roots = [
    ROOT / "data/is/reference/gigastroke",
    ROOT / "data/is",
]

files = []

for root in roots:
    if not root.exists():
        continue

    for p in root.rglob("*"):
        if not p.is_file():
            continue

        if p.suffix.lower() in {".json", ".yaml", ".yml"}:
            files.append(p)

files = sorted(set(files))

wanted = {
    "accession",
    "ancestry",
    "sample_ancestry",
    "population",
    "trait",
    "phenotype",
    "sample_size",
    "n",
    "cases",
    "controls",
    "n_cases",
    "n_controls",
    "genome_build",
    "build",
    "is_harmonised",
    "is_harmonized",
    "file_type",
    "source",
}

rows = []
raw_extract = []

def emit(source, key, value, parsed_via):
    sval = "" if value is None else str(value)

    rows.append({
        "source": str(source),
        "accession": "",
        "key": str(key),
        "value": sval,
        "parsed_via": parsed_via,
    })

def flatten(obj, prefix="", source=None, parsed_via=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)

            low = str(k).lower()

            if (
                low in wanted
                or any(x in low for x in [
                    "ancestry",
                    "sample_size",
                    "case",
                    "control",
                    "build",
                    "accession",
                    "trait",
                    "phenotype",
                    "harmonis",
                ])
            ):
                emit(source, key, v, parsed_via)

            flatten(v, key, source, parsed_via)

    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            flatten(v, f"{prefix}[{i}]", source, parsed_via)

try:
    import yaml
    has_yaml = True
except Exception:
    has_yaml = False

for p in files:

    text = ""

    try:
        text = p.read_text(errors="replace")
    except Exception:
        continue

    for ln, line in enumerate(text.splitlines(), start=1):

        low = line.lower()

        if any(k in low for k in [
            "accession",
            "ancestry",
            "sample_ancestry",
            "sample_size",
            "case",
            "control",
            "genome_build",
            "build:",
            "trait",
            "phenotype",
            "harmonis",
        ]):
            raw_extract.append(
                f"{p}\tLINE={ln}\t{line}"
            )

    parsed = False

    if p.suffix.lower() == ".json":
        try:
            obj = json.loads(text)
            flatten(obj, source=p, parsed_via="json")
            parsed = True
        except Exception:
            pass

    elif p.suffix.lower() in {".yaml", ".yml"} and has_yaml:
        try:
            docs = list(yaml.safe_load_all(text))
            for obj in docs:
                flatten(obj, source=p, parsed_via="yaml")
            parsed = True
        except Exception:
            pass

    if not parsed:
        for ln, line in enumerate(text.splitlines(), start=1):

            m = re.match(
                r"^\s*([A-Za-z0-9_.-]+)\s*:\s*(.*?)\s*$",
                line
            )

            if not m:
                continue

            k, v = m.groups()
            low = k.lower()

            if (
                low in wanted
                or any(x in low for x in [
                    "ancestry",
                    "sample_size",
                    "case",
                    "control",
                    "build",
                    "accession",
                    "trait",
                    "phenotype",
                    "harmonis",
                ])
            ):
                emit(
                    p,
                    f"{k}@line{ln}",
                    v,
                    "raw_regex"
                )

out = OUT / "GIGASTROKE_METADATA_LOCAL_AUDIT.tsv"

with out.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "source",
            "accession",
            "key",
            "value",
            "parsed_via",
        ],
    )
    w.writeheader()
    w.writerows(rows)

(OUT / "GIGASTROKE_METADATA_RAW_EXTRACT.txt").write_text(
    "\n".join(raw_extract) + "\n"
)

print("metadata_files =", len(files))
print("metadata_rows  =", len(rows))
print("raw_extract    =", len(raw_extract))

ancestry_rows = [
    x for x in rows
    if "ancestry" in x["key"].lower()
]

print()
print("===== ANCESTRY VALUES =====")

for x in ancestry_rows[:100]:
    print(
        x["source"],
        x["key"],
        repr(x["value"])
    )

PY

}


# ============================================================
# 2. FASTA CANDIDATE INVENTORY — NO AUTO SELECTION
# ============================================================

fasta_inventory() {

python3 - <<'PY'
from pathlib import Path
import csv
import os

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

roots = [
    ROOT / "data",
    ROOT / "reference",
    ROOT / "resources",
]

MAX_DEPTH = 6

rows = []

for root in roots:

    if not root.exists():
        continue

    root_depth = len(root.parts)

    for current, dirs, files in os.walk(root):

        cur = Path(current)

        depth = len(cur.parts) - root_depth

        if depth >= MAX_DEPTH:
            dirs[:] = []

        for fn in files:

            low = fn.lower()

            if not (
                low.endswith(".fa")
                or low.endswith(".fasta")
                or low.endswith(".fa.gz")
                or low.endswith(".fasta.gz")
            ):
                continue

            p = cur / fn

            fai_candidates = [
                Path(str(p) + ".fai"),
            ]

            if str(p).endswith(".gz"):
                fai_candidates.append(
                    Path(str(p)[:-3] + ".fai")
                )

            fai = next(
                (
                    x for x in fai_candidates
                    if x.exists()
                ),
                None
            )

            contigs = ""

            if fai is not None:
                try:
                    contigs = ",".join(
                        [
                            line.split("\t")[0]
                            for line in fai.read_text(
                                errors="replace"
                            ).splitlines()[:10]
                       ]
                    )
                except Exception:
                    pass

            rows.append({
                "path": str(p),
                "bytes": p.stat().st_size,
                "indexed": int(fai is not None),
                "fai_path": "" if fai is None else str(fai),
                "contig_preview": contigs,
                "selection_status": "CANDIDATE_ONLY_NOT_SELECTED",
            })

with (OUT / "FASTA_CANDIDATES.tsv").open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "path",
            "bytes",
            "indexed",
            "fai_path",
            "contig_preview",
            "selection_status",
        ],
    )

    w.writeheader()
    w.writerows(rows)

print("FASTA candidates =", len(rows))

for r in rows:
    print(
        r["path"],
        "indexed=",
        r["indexed"],
        "contigs=",
        r["contig_preview"]
    )

PY

}


# ============================================================
# 3. susieR API AUDIT — NO ANALYTICAL CALLS
# ============================================================

susier_api_audit() {

Rscript - <<'RS'

OUT <- "/srv/is-analysis/results/is/stage4_cross_eas/phase8d_validation_audit"
F <- file.path(OUT, "SUSIER_API_AUDIT.txt")

sink(F)
sink(F, type="message")

ok <- TRUE

cat("===== susieR API AUDIT =====\n\n")

if (!requireNamespace("susieR", quietly=TRUE)) {

    cat("susieR=NOT_INSTALLED\n")
    ok <- FALSE

} else {

    cat(
        "susieR_version=",
        as.character(packageVersion("susieR")),
        "\n\n",
        sep=""
    )

    funcs <- c(
        "susie_rss",
        "estimate_s_rss",
        "kriging_rss"
    )

    for (fn in funcs) {

        cat(
            "\n==================================================\n",
            fn,
            "\n==================================================\n",
            sep=""
        )

        obj <- tryCatch(
            getExportedValue("susieR", fn),
            error=function(e) NULL
        )

        if (is.null(obj)) {

            cat("FUNCTION_NOT_FOUND\n")
            ok <- FALSE
            next
        }

        cat("\nFORMALS\n")
        print(formals(obj))

        cat("\nBODY_PREVIEW\n")

        bb <- tryCatch(
            deparse(body(obj)),
            error=function(e) character()
        )

        if (length(bb) > 0) {
            cat(
                paste(
                    head(bb, 80),
                    collapse="\n"
                ),
                "\n"
            )
        }

        cat("\nGETANYWHERE\n")

        ga <- tryCatch(
            getAnywhere(fn),
            error=function(e) NULL
        )

        if (!is.null(ga)) {
            print(ga)
        }
    }
}

sink(type="message")
sink()

marker <- file.path(
    OUT,
    "SUSIER_API_STATUS.txt"
)

writeLines(
    if (ok) "AUDITED" else "PENDING_API_AUDIT",
    marker
)

cat(
    "SUSIER_API=",
    if (ok) "AUDITED" else "PENDING_API_AUDIT",
    "\n",
    sep=""
)

RS

}


# ============================================================
# 4. HARMONIZED INPUT AUDIT
# ============================================================

harmonized_audit() {

python3 - <<'PY'
from pathlib import Path
import gzip
import csv
import re

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

search_roots = [
    ROOT / "results/is/stage4_cross_eas",
    ROOT / "data/is/reference/gigastroke",
]

candidates = []

for root in search_roots:

    if not root.exists():
        continue

    for p in root.rglob("*"):

        if not p.is_file():
            continue

        low = str(p).lower()

        if not (
            low.endswith(".tsv")
            or low.endswith(".tsv.gz")
            or low.endswith(".txt")
            or low.endswith(".txt.gz")
        ):
            continue

        if any(
            x in low
            for x in [
                "harmon",
                "aligned",
                "phase6",
                "gigastroke",
           ]
        ):
            candidates.append(p)

candidates = sorted(set(candidates))

alignment_names = {
    "alignment_status",
    "align_status",
    "harmonization_status",
    "harmonisation_status",
    "alignment",
    "harmonization",
    "harmonisation",
    "match_type",
}

se_names = {
    "aligned_se",
    "se",
    "standard_error",
    "stderr",
}

rows = []
headers_out = []

for p in candidates:

    opener = gzip.open if str(p).endswith(".gz") else open

    try:
        fh = opener(
            p,
            "rt",
            errors="replace"
        )
    except Exception:
        continue

    with fh:

        first = fh.readline().rstrip("\n")

        if not first:
            continue

        delim = "\t" if "\t" in first else None

        if delim is None:
            continue

        header = first.split("\t")

        headers_out.append(
            f"\n===== {p} =====\n"
            + "\t".join(header)
            + "\n"
        )

        lower = {
            c.lower(): c
            for c in header
        }

        align_col = next(
            (
                lower[x]
                for x in alignment_names
                if x in lower
            ),
            None
        )

        se_col = next(
            (
                lower[x]
                for x in se_names
                if x in lower
            ),
            None
        )

        align_counts = {}
        n = 0

        if align_col is not None:

            idx = header.index(align_col)

            for line in fh:

                if not line.strip():
                    continue

                parts = line.rstrip("\n").split("\t")

                if idx >= len(parts):
                    continue

                v = parts[idx].strip()
                align_counts[v] = align_counts.get(v, 0) + 1
                n += 1

        rows.append({
            "path": str(p),
            "alignment_column": "" if align_col is None else align_col,
            "se_column": "" if se_col is None else se_col,
            "rows_scanned": n,
            "alignment_counts": ";".join(
                f"{k}:{v}"
                for k, v in sorted(
                    align_counts.items()
                )
            ),
            "contains_complement_status": int(
                any(
                    "COMPLEMENT" in k.upper()
                    for k in align_counts
                )
            ),
            "audit_status": (
                "ALIGNMENT_COLUMN_FOUND"
                if align_col is not None
                else "NO_ALIGNMENT_COLUMN"
            ),
        })

with (OUT / "HARMONIZED_INPUT_AUDIT.tsv").open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "path",
            "alignment_column",
            "se_column",
            "rows_scanned",
            "alignment_counts",
            "contains_complement_status",
            "audit_status",
        ],
    )

    w.writeheader()
    w.writerows(rows)

(OUT / "HARMONIZED_HEADERS.txt").write_text(
    "".join(headers_out)
)

print("candidate_files =", len(candidates))
print("audited_files   =", len(rows))

for r in rows:
    print(
        r["path"],
        r["audit_status"],
        r["alignment_counts"],
        "complement=",
        r["contains_complement_status"]
    )

PY

}


# ============================================================
# 5. RAW LD STRUCTURAL QC
# ============================================================

ld_structural_qc() {

python3 - <<'PY'
from pathlib import Path
import numpy as np
import csv
import math

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

search_roots = [
    ROOT / "results/is/stage3_finemap",
    ROOT / "results/is/stage4_cross_eas",
    ROOT / "data/is",
]

bins = []

for root in search_roots:

    if not root.exists():
        continue

    for p in root.rglob("*.bin"):

        low = p.name.lower()

        if not any(
            x in low
            for x in [
                "l001",
                "l002",
                "l003",
                "l004",
            ]
        ):
            continue

        bins.append(p)

bins = sorted(set(bins))

rows = []

for p in bins:

    var_candidates = [
        Path(str(p) + ".vars"),
        p.with_suffix(".vars"),
    ]

    vars_path = next(
        (
            x for x in var_candidates
            if x.exists()
        ),
        None
    )

    if vars_path is None:

        rows.append({
            "matrix_path": str(p),
            "vars_path": "",
            "n": "",
            "actual_bytes": p.stat().st_size,
            "expected_float64_bytes": "",
            "format_status": "VARS_NOT_FOUND",
            "max_asym": "",
            "diag_min": "",
            "diag_max": "",
            "max_diag_absdev1": "",
            "matrix_min": "",
            "matrix_max": "",
            "nonfinite": "",
        })

        continue

    with vars_path.open(
        errors="replace"
    ) as f:
        n = sum(
            1 for x in f
            if x.strip()
        )

    actual = p.stat().st_size

    exp64 = n * n * 8
    exp32 = n * n * 4

    if actual == exp64:
        dtype = np.dtype("<f8")
        fmt = "FLOAT64_EXACT"

    elif actual == exp32:
        dtype = None
        fmt = "FLOAT32_SIZE_MATCH_NO_MEMMAP"

    else:
        dtype = None
        fmt = "FORMAT_UNRESOLVED"

    stats = {
        "max_asym": "",
        "diag_min": "",
        "diag_max": "",
        "max_diag_absdev1": "",
        "matrix_min": "",
        "matrix_max": "",
        "nonfinite": "",
    }

    if dtype is not None:

        M = np.memmap(
            p,
            dtype=dtype,
            mode="r",
            shape=(n, n)
        )

        d = np.asarray(
            M[
                np.arange(n),
                np.arange(n)
            ]
        )

        stats["diag_min"] = float(
            np.nanmin(d)
        )

        stats["diag_max"] = float(
            np.nanmax(d)
        )

        stats["max_diag_absdev1"] = float(
            np.nanmax(
                np.abs(d - 1.0)
            )
        )

        matrix_min = math.inf
        matrix_max = -math.inf
        nonfinite = 0
        max_asym = 0.0

        block = 512

        for i0 in range(0, n, block):

            i1 = min(n, i0 + block)

            A = np.asarray(
                M[i0:i1, :]
            )

            finite = np.isfinite(A)

            nonfinite += int(
                A.size - finite.sum()
            )

            if finite.any():

                matrix_min = min(
                    matrix_min,
                    float(
                        np.nanmin(A)
                    )
                )

                matrix_max = max(
                    matrix_max,
                    float(
                        np.nanmax(A)
                    )
                )

            B = np.asarray(
                M[:, i0:i1]
            ).T

            diff = np.abs(
                A - B
            )

            if np.isfinite(diff).any():

                max_asym = max(
                    max_asym,
                    float(
                        np.nanmax(diff)
                    )
                )

        stats["matrix_min"] = matrix_min
        stats["matrix_max"] = matrix_max
        stats["nonfinite"] = nonfinite
        stats["max_asym"] = max_asym

    rows.append({
        "matrix_path": str(p),
        "vars_path": str(vars_path),
        "n": n,
        "actual_bytes": actual,
        "expected_float64_bytes": exp64,
        "format_status": fmt,
        **stats,
    })

with (OUT / "LD_STRUCTURAL_QC.tsv").open(
    "w",
    newline=""
) as f:

    fields = [
        "matrix_path",
        "vars_path",
        "n",
        "actual_bytes",
        "expected_float64_bytes",
        "format_status",
        "max_asym",
        "diag_min",
        "diag_max",
        "max_diag_absdev1",
        "matrix_min",
        "matrix_max",
        "nonfinite",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t,
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(rows)

print("LD matrices audited =", len(rows))

for r in rows:
    print(
        r["matrix_path"],
        r["format_status"],
        "n=",
        r["n"],
        "asym=",
        r["max_asym"],
        "diagdev=",
        r["max_diag_absdev1"]
    )

PY

}


# ============================================================
# 6. SCRIPT PROVENANCE / LD ORIENTATION / EXACT SUBSET / BAD ROWS
# ============================================================

provenance_and_orientation() {

python3 - <<'PY'
from pathlib import Path
import hashlib
import csv
import re
import os

ROOT = Path("/srv/is-analysis")
REPO = ROOT / "IS_Analysis_V3"
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

# ------------------------------------------------------------
# Script provenance
# ------------------------------------------------------------

scripts = []

for base in [
    REPO,
    ROOT / "server",
]:

    if not base.exists():
        continue

    for p in base.rglob("*"):

        if not p.is_file():
            continue

        low = p.name.lower()

        if any(
            x in low
            for x in [
                "phase5",
                "phase8b",
                "phase8c",
                "finemap",
                "ld",
                "harmon",
            ]
        ):
            if p.suffix.lower() in {
                ".sh",
                ".py",
                ".r",
                ".R".lower(),
            }:
                scripts.append(p)

scripts = sorted(set(scripts))

prov_lines = []

orientation_rows = []

patterns = [
    "plink",
    "vcor",
    "r-unphased",
    "ref",
    "alt",
    "a1",
    "a2",
    "effect",
    "other",
    "allele",
    "swap",
    "complement",
    "subset",
    "l003",
    "l004",
]

for p in scripts:

    try:
        b = p.read_bytes()
        text = b.decode(
            errors="replace"
        )
    except Exception:
        continue

    sha = hashlib.sha256(
        b
    ).hexdigest()

    prov_lines.append(
        f"\n===== {p} =====\n"
        f"sha256={sha}\n"
        f"bytes={len(b)}\n"
        f"mtime={p.stat().st_mtime}\n"
    )

    for ln, line in enumerate(
        text.splitlines(),
        start=1
    ):

        low = line.lower()

        if any(
            pat in low
            for pat in patterns
        ):

            prov_lines.append(
                f"{ln}: {line}\n"
            )

            if any(
                x in low
                for x in [
                    "plink",
                    "vcor",
                    "r-unphased",
                    "ref",
                    "alt",
                    "a1",
                    "a2",
                    "effect",
                    "other",
                    "swap",
                    "complement",
                ]
            ):

                orientation_rows.append({
                    "evidence_type": "SCRIPT_LINE",
                    "source": str(p),
                    "evidence": f"line {ln}: {line}",
                    "status": "OBSERVED_NOT_INFERRED",
                })

(OUT / "SCRIPT_PROVENANCE.txt").write_text(
    "".join(prov_lines)
)

# ------------------------------------------------------------
# PVAR orientation evidence
# ------------------------------------------------------------

for base in [
    ROOT / "results/is/stage3_finemap",
    ROOT / "data/is",
]:

    if not base.exists():
        continue

    for p in base.rglob("*.pvar"):

        low = str(p).lower()

        if not any(
            x in low
            for x in [
                "l001",
                "l002",
                "l003",
                "l004",
                "1000g",
                "bbj",
            ]
        ):
            continue

        lines = []

        try:
            with p.open(
                errors="replace"
            ) as f:
                for _ in range(8):
                    line = f.readline()
                    if not line:
                        break
                    lines.append(
                        line.rstrip("\n")
                    )
        except Exception:
            continue

        orientation_rows.append({
            "evidence_type": "PVAR_PREVIEW",
            "source": str(p),
            "evidence": " || ".join(lines),
            "status": "OBSERVED_NOT_INFERRED",
        })

with (OUT / "LD_ORIENTATION_INVENTORY.tsv").open(
    "w",
    newline=""
) as f:

    fields = [
        "evidence_type",
        "source",
        "evidence",
        "status",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(
        orientation_rows
    )

# ------------------------------------------------------------
# Exact diagnostic subset matrices
# ------------------------------------------------------------

extensions = {
    ".bin",
    ".npy",
    ".rds",
    ".rdata",
    ".mtx",
}

targets = {
    "L003_BBJ": [],
    "L004_AIS": [],
}

for base in [
    ROOT / "results/is/stage4_cross_eas",
    ROOT / "results/is/stage3_finemap",
]:

    if not base.exists():
        continue

    for p in base.rglob("*"):

        if not p.is_file():
            continue

        if p.suffix.lower() not in extensions:
            continue

        low = str(p).lower()

        if "l003" in low and "bbj" in low:
            targets["L003_BBJ"].append(p)

        if "l004" in low and "ais" in low:
            targets["L004_AIS"].append(p)

exact_rows = []

for target, paths in targets.items():

    if paths:

        for p in sorted(
            set(paths)
        ):

            exact_rows.append({
                "target": target,
                "path": str(p),
                "bytes": p.stat().st_size,
                "status": "CANDIDATE_MATRIX_FOUND_NEEDS_SCRIPT_CONFIRMATION",
            })

    else:

        exact_rows.append({
            "target": target,
            "path": "",
            "bytes": "",
            "status": "NOT_PERSISTED_RECONSTRUCTION_PENDING_SCRIPT_AUDIT",
        })

with (OUT / "EXACT_SUBSET_MATRIX_AUDIT.tsv").open(
    "w",
    newline=""
) as f:

    fields = [
        "target",
        "path",
        "bytes",
        "status",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(exact_rows)

# ------------------------------------------------------------
# BBJ bad-row artifacts
# ------------------------------------------------------------

bad_rows = []

for base in [
    ROOT / "data/is/processed/japan/bbj",
    ROOT / "results/is/stage3_finemap/japan/bbj",
    ROOT / "data/is/east_asia/japan/bbj",
]:

    if not base.exists():
        continue

    for p in base.rglob("*"):

        if not p.is_file():
            continue

        low = p.name.lower()

        if any(
            x in low
            for x in [
                "bad",
                "reject",
                "invalid",
                "fail",
            ]
        ):

            bad_rows.append({
                "path": str(p),
                "bytes": p.stat().st_size,
                "mtime": p.stat().st_mtime,
                "status": "FOUND",
            })

if not bad_rows:

    bad_rows.append({
        "path": "",
        "bytes": "",
        "mtime": "",
        "status": "NO_NAMED_BAD_ROW_ARTIFACT_FOUND",
    })

with (OUT / "BBJ_BAD_ROW_ARTIFACT_INVENTORY.tsv").open(
    "w",
    newline=""
) as f:

    fields = [
        "path",
        "bytes",
        "mtime",
        "status",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(
        bad_rows
    )

print("scripts =", len(scripts))
print(
    "orientation evidence =",
    len(orientation_rows)
)

for r in exact_rows:
    print(
        "exact subset:",
        r["target"],
        r["status"],
        r["path"]
    )

PY

}


# ============================================================
# 7. SAMPLE SIZE AUDIT
# ============================================================

sample_size_audit() {

python3 - <<'PY'
from pathlib import Path
import csv
import re

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

metadata = OUT / "GIGASTROKE_METADATA_LOCAL_AUDIT.tsv"

rows = []

if metadata.exists():

    with metadata.open() as f:

        r = csv.DictReader(
            f,
            delimiter="\t"
        )

        for x in r:

            key = x["key"].lower()

            if any(
                z in key
                for z in [
                    "sample_size",
                    "case",
                    "control",
                    "ancestry",
                    "accession",
                    "trait",
                    "phenotype",
                ]
            ):

                rows.append({
                    "source": x["source"],
                    "key": x["key"],
                    "value": x["value"],
                    "interpretation": "LOCAL_METADATA_OBSERVATION_ONLY",
                })

with (OUT / "SAMPLE_SIZE_AUDIT.tsv").open(
    "w",
    newline=""
) as f:

    fields = [
        "source",
        "key",
        "value",
        "interpretation",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(
        rows
    )

print("sample-size metadata rows =", len(rows))

PY

}


# ============================================================
# 8. SCIENTIFIC READINESS
# ============================================================

write_readiness() {

python3 - <<'PY'
from pathlib import Path
import csv

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/is/stage4_cross_eas/phase8d_validation_audit"

api_status = "PENDING_API_AUDIT"

p = OUT / "SUSIER_API_STATUS.txt"

if p.exists():
    x = p.read_text().strip()
    if x == "AUDITED":
        api_status = "AUDITED"

phase10_output = False

for base in [
    ROOT / "results/is",
]:

    if not base.exists():
        continue

    for p in base.rglob("*"):

        if not p.is_file():
            continue

        low = str(p).lower()

        if "phase10" in low and p.stat().st_size > 0:
            phase10_output = True
            break

    if phase10_output:
        break

rows = [
    (
        "PHASE7B_GIGASTROKE_SUSIE",
        "EXPLORATORY_SAMPLE_SIZE_MISSING"
    ),
    (
        "PHASE8_GIGASTROKE_NAWARE",
        "INVALID_N_ASSIGNMENT"
    ),
    (
        "PHASE8B_GIGASTROKE_NAWARE",
        "EXPLORATORY_CORRECT_N_BUT_VALIDATION_PENDING"
    ),
    (
        "PHASE8C_MASTER_REBUILD",
        "BOOKKEEPING_COMPLETE"
    ),
    (
        "BBJ_GIGASTROKE_CS_COMPARE",
        "EXPLORATORY_ANCESTRY_REF_PENDING"
    ),
    (
        "BBJ_PRIMARY_FINEMAP",
        "EXPLORATORY_EAS504"
    ),
    (
        "GIGASTROKE_ANCESTRY",
        "PENDING_OFFICIAL_VERIFICATION"
    ),
    (
        "GIGASTROKE_REF_ALT",
        "PENDING_REFERENCE_AUDIT"
    ),
    (
        "BBJ_REF_FASTA",
        "PENDING_REFERENCE_AUDIT"
    ),
    (
        "SUSIER_RSS_API",
        api_status
    ),
    (
        "L003_CONDITIONING",
        "PENDING"
    ),
    (
        "TPMI",
        "WAITING_RATE_LIMIT"
    ),
    (
        "CKB",
        "WAITING_DECRYPTION_KEY"
    ),
    (
        "PHASE9A_OUTPUTS",
        "EXPLORATORY_FROZEN_CORE_VALIDATION_PENDING"
    ),
    (
        "PHASE9B_OUTPUTS",
        "EXPLORATORY_FROZEN_CORE_VALIDATION_PENDING"
    ),
]

if phase10_output:

    rows.append(
        (
            "PHASE10_OUTPUTS",
            "EXPLORATORY_FROZEN_DO_NOT_INTERPRET"
        )
    )

else:

    rows.append(
        (
            "PHASE10",
            "DO_NOT_RUN_BEFORE_CORE_VALIDATION"
        )
    )

# Phase11B is also frozen because it was executed before core validation.
rows.append(
    (
        "PHASE11B_OUTPUTS",
        "EXPLORATORY_FROZEN_CORE_VALIDATION_PENDING"
    )
)

with (OUT / "SCIENTIFIC_READINESS.tsv").open(
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

for k, v in rows:
    print(
        f"{k} = {v}"
    )

PY

}


# ============================================================
# RUN AUDIT STEPS
# ============================================================

run_step "01_metadata_audit" \
  metadata_audit

run_step "02_fasta_inventory" \
  fasta_inventory

run_step "03_susier_api_audit" \
  susier_api_audit

run_step "04_harmonized_input_audit" \
  harmonized_audit

run_step "05_ld_structural_qc" \
  ld_structural_qc

run_step "06_provenance_orientation" \
  provenance_and_orientation

run_step "07_sample_size_audit" \
  sample_size_audit

run_step "08_scientific_readiness" \
  write_readiness


# ============================================================
# FINAL SUMMARY
# ============================================================

echo
echo "================================================================================"
echo "PHASE 8D AUDIT COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="
column -t -s $'\t' "$STATUS" 2>/dev/null || cat "$STATUS"

echo
echo "===== SCIENTIFIC READINESS ====="
column -t -s $'\t' \
  "$OUT/SCIENTIFIC_READINESS.tsv" \
  2>/dev/null \
  || cat "$OUT/SCIENTIFIC_READINESS.tsv"

echo
echo "===== OUTPUT FILES ====="

for f in \
  GIGASTROKE_METADATA_LOCAL_AUDIT.tsv \
  GIGASTROKE_METADATA_RAW_EXTRACT.txt \
  FASTA_CANDIDATES.tsv \
  HARMONIZED_INPUT_AUDIT.tsv \
  HARMONIZED_HEADERS.txt \
  SUSIER_API_AUDIT.txt \
  LD_STRUCTURAL_QC.tsv \
  LD_ORIENTATION_INVENTORY.tsv \
  BBJ_BAD_ROW_ARTIFACT_INVENTORY.tsv \
  SCRIPT_PROVENANCE.txt \
  SAMPLE_SIZE_AUDIT.tsv \
  SCIENTIFIC_READINESS.tsv \
  EXACT_SUBSET_MATRIX_AUDIT.tsv
do

  P="$OUT/$f"

  if [[ -f "$P" ]]; then
    printf "FOUND\t%s\t%s bytes\n" \
      "$P" \
      "$(stat -c%s "$P" 2>/dev/null || echo NA)"
  else
    printf "MISSING\t%s\n" "$P"
  fi

done

echo
echo "NOTE:"
echo "STEP_STATUS=PASS means the audit code executed successfully."
echo "It does NOT mean scientific validation passed."
echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
