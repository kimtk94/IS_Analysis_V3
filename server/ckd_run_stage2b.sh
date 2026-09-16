#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STAGE2="${CKD_STAGE2_ROOT:-/srv/is-analysis/results/ckd/stage2}"
OUT="${CKD_STAGE2B_ROOT:-/srv/is-analysis/results/ckd/stage2b_coloc}"
PQTL="${CKD_STAGE2_PQTL_ROOT:-/srv/is-analysis/data/ckd/stage2_pqtl}"
REFROOT="${CKD_REFERENCE_ROOT:-/srv/is-analysis/data/reference}"
COORDS="${CKD_GENE_COORDS:-$REFROOT/gene_coordinates_hg38.tsv}"
CHAIN="${CKD_HG38_TO_HG19_CHAIN:-$REFROOT/hg38ToHg19.over.chain}"
EGFR="${CKD_EGFR_RAW:-/srv/is-analysis/data/ckd/rawdata/outcome/EUR/ckdgen_stanzick2021_egfr_eur/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz}"
PYTHON="${CKD_PYTHON:-$ROOT/.venv-ckd/bin/python}"
R_LIB="${CKD_R_LIB:-/srv/is-analysis/.Rlib}"

test -x "$PYTHON"
test -f "$STAGE2/stage2_candidates.tsv"
test -f "$EGFR"
mkdir -p "$OUT" "$REFROOT" "$R_LIB"

if [[ ! -f "$COORDS" ]]; then
  command -v rclone >/dev/null
  rclone copyto "gdrive:IS_Analysis_V3/data/reference/gene_coordinates_hg38.tsv" "$COORDS"
fi

if [[ ! -f "$CHAIN" ]]; then
  command -v curl >/dev/null
  tmp="$CHAIN.gz"
  curl -L --fail --retry 5 -o "$tmp" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
  gzip -dc "$tmp" > "$CHAIN"
  rm -f "$tmp"
fi

if ! "$PYTHON" -c 'import pyliftover' >/dev/null 2>&1; then
  "$PYTHON" -m pip install --disable-pip-version-check pyliftover
fi

"$PYTHON" "$ROOT/scripts/prepare_ckd_stage2b_coloc.py" \
  --candidates "$STAGE2/stage2_candidates.tsv" \
  --pqtl-root "$PQTL" \
  --gene-coordinates "$COORDS" \
  --chain "$CHAIN" \
  --outcome-egfr "$EGFR" \
  --output-root "$OUT"

command -v Rscript >/dev/null || { echo "Rscript is required for coloc" >&2; exit 7; }
export R_LIBS_USER="$R_LIB"
Rscript -e 'if (!requireNamespace("coloc", quietly=TRUE)) install.packages("coloc", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'
Rscript "$ROOT/scripts/run_ckd_stage2b_coloc.R" "$OUT/coloc_input" "$OUT/coloc_results"

(
  cd "$OUT"
  find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)

if [[ "${SYNC_DRIVE:-0}" == "1" ]]; then
  RCLONE_REMOTE="${RCLONE_REMOTE:-gdrive}"
  rclone copy "$OUT" "${RCLONE_REMOTE}:IS_Analysis_V3/results/ckd/stage2b_coloc" \
    --checksum --transfers 1 --checkers 2 --progress
  echo "CKD_STAGE2B_DRIVE_SYNC_PASS"
fi

echo "CKD_STAGE2B_PASS"
echo "RESULT=$OUT"
