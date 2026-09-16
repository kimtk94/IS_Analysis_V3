#!/usr/bin/env bash
set -euo pipefail

# Reproducible public-data acquisition stage for the CKD thesis.
# Restricted KoGES data are intentionally excluded and must remain in the
# approved analysis environment.

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RAW="${CKD_DATA_ROOT:-$ROOT/data/rawdata/ckd}"
READY="${CKD_READY_ROOT:-$ROOT/data/analysis_ready/ckd}"
BASE_PYTHON="${PYTHON:-python3}"
VENV="${CKD_VENV:-$ROOT/.venv-ckd}"

command -v curl >/dev/null
command -v "$BASE_PYTHON" >/dev/null

if [[ ! -x "$VENV/bin/python" ]]; then
  echo "[env] creating CKD virtual environment: $VENV"
  if ! "$BASE_PYTHON" -m venv "$VENV"; then
    cat >&2 <<'EOF'
ERROR: Python venv support is not installed.
On Ubuntu/Debian install it once, then rerun this script:

  sudo apt update
  sudo apt install -y python3-venv

Do not use --break-system-packages.
EOF
    exit 2
  fi
fi

PYTHON="$VENV/bin/python"

if ! "$PYTHON" -c 'import openpyxl' 2>/dev/null; then
  echo "[env] installing openpyxl into $VENV"
  "$PYTHON" -m pip install --disable-pip-version-check 'openpyxl==3.1.5'
fi

echo "[env] python=$PYTHON"
"$PYTHON" -c 'import openpyxl; print("[env] openpyxl=" + openpyxl.__version__)'

if [[ "${CKD_ENV_ONLY:-0}" == "1" ]]; then
  echo "CKD_ENV_PASS"
  exit 0
fi

mkdir -p "$RAW" "$READY" "$RAW/eas" "$RAW/instruments"

echo "[1/6] EUR public kidney outcomes"
"$PYTHON" "$ROOT/scripts/download_ckd_public_data.py" --dest "$RAW" \
  --data-id ckdgen_stanzick2021_egfr_eur \
  --data-id ckdgen_wuttke2019_ckd_eur \
  --data-id ckdgen_wuttke2019_bun_eur \
  --data-id ckdgen_gorski2017_egfrcys_eur \
  --data-id ckdgen_teumer2019_uacr_eur

echo "[2/6] EAS eGFR/BUN"
curl -L --fail --retry 5 --retry-all-errors https://ndownloader.figshare.com/files/42774049 -o "$RAW/eas/eGFR.gz"
curl -L --fail --retry 5 --retry-all-errors https://ndownloader.figshare.com/files/43218600 -o "$RAW/eas/BUN.gz"
gzip -t "$RAW/eas/eGFR.gz"
gzip -t "$RAW/eas/BUN.gz"

echo "[3/6] UKB-PPP Sun 2023 cis instruments"
curl -L --fail --retry 5 --retry-all-errors \
  'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx' \
  -o "$RAW/instruments/Sun2023.xlsx"
"$PYTHON" "$ROOT/scripts/extract_ukbppp_st16_cis.py" \
  --xlsx "$RAW/instruments/Sun2023.xlsx" \
  --output "$READY/UKBPPP_ST16_cis_independent.tsv.gz" \
  --summary "$READY/UKBPPP_ST16_cis_independent.summary.json"

I="$READY/UKBPPP_ST16_cis_independent.tsv.gz"
match () {
  local schema="$1" input="$2" phenotype="$3" ancestry="$4"
  mkdir -p "$READY/$ancestry"
  "$PYTHON" "$ROOT/scripts/extract_ckd_matched_outcomes.py" \
    --schema "$schema" --input "$input" --instruments "$I" \
    --output "$READY/$ancestry/$phenotype.tsv.gz" \
    --summary "$READY/$ancestry/$phenotype.summary.json" \
    --phenotype "$phenotype" --ancestry "$ancestry"
}

echo "[4/6] EUR instrument-matched outcomes"
match stanzick_egfr "$RAW/outcome/EUR/ckdgen_stanzick2021_egfr_eur/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz" eGFRcrea EUR
match wuttke_ckd "$RAW/outcome/EUR/ckdgen_wuttke2019_ckd_eur/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz" CKD EUR
match wuttke_bun "$RAW/support/EUR/ckdgen_wuttke2019_bun_eur/BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt.gz" BUN EUR
match gorski_egfrcys "$RAW/support/EUR/ckdgen_gorski2017_egfrcys_eur/CKDGen_1000Genomes_DiscoveryMeta_eGFRcys_overall.csv.gz" eGFRcys EUR
match teumer_uacr "$RAW/support/EUR/ckdgen_teumer2019_uacr_eur/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz" UACR EUR

echo "[5/6] EAS instrument-matched outcomes"
match chen_egfr "$RAW/eas/eGFR.gz" eGFRcrea EAS
match chen_bun "$RAW/eas/BUN.gz" BUN EAS

echo "[6/6] Integrity inventory"
find "$READY" -type f ! -name 'SHA256SUMS.txt' -print0 | sort -z | xargs -0 sha256sum > "$READY/SHA256SUMS.txt"
"$PYTHON" - "$READY" <<'PY'
import json, sys
from pathlib import Path
root=Path(sys.argv[1])
summaries=[]
for p in sorted(root.rglob("*.summary.json")):
    if p.name.startswith("UKBPPP_"):
        continue
    x=json.loads(p.read_text())
    if "matched_rows" in x:
        summaries.append(x)
assert len(summaries)==7, len(summaries)
assert all(x["matched_rows"]>0 for x in summaries)
(root/"MATCH_SUMMARY.json").write_text(json.dumps(summaries,indent=2)+"\n")
print(json.dumps(summaries,indent=2))
PY

echo "CKD_DATA_STAGE_PASS"
echo "RAW=$RAW"
echo "READY=$READY"
