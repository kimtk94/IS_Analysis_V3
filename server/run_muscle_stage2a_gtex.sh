#!/usr/bin/env bash

ROOT="/srv/is-analysis/IS_Analysis_V3"
WORKTREE_ROOT="/srv/is-analysis/worktrees"
STAMP=$(date '+%Y%m%d_%H%M%S')

MUSCLE_REPO="$WORKTREE_ROOT/MUSCLE_STAGE2A_$STAMP"
STAGE1="/srv/is-analysis/results/muscle/stage1_gwas_audit/MUSCLE_STAGE1_GWAS_AUDIT.tsv"
OUT="/srv/is-analysis/results/muscle/stage2a_gtex"
LOG="/srv/is-analysis/results/muscle/logs"
RUNLOG="$LOG/MUSCLE_STAGE2A_$STAMP.log"

mkdir -p "$WORKTREE_ROOT" "$OUT" "$LOG"

echo "============================================================"
echo " MUSCLE STAGE 2A — GTEx MUSCLE eQTL/sQTL"
echo "============================================================"

echo
echo "===== 1. PRESERVE CURRENT CKD/IS WORKTREE ====="
git -C "$ROOT" status --short | head -100

echo
echo "===== 2. FETCH REMOTE MAIN ====="
git -C "$ROOT" fetch origin main
FETCH_RC=$?
echo "FETCH_RC=$FETCH_RC"
git -C "$ROOT" rev-parse origin/main

echo
echo "===== 3. CREATE CLEAN STAGE2A WORKTREE ====="
git -C "$ROOT" worktree add --detach "$MUSCLE_REPO" origin/main
WT_RC=$?
echo "WORKTREE_RC=$WT_RC"
echo "MUSCLE_REPO=$MUSCLE_REPO"

echo
echo "===== 4. INPUT / CODE CHECK ====="
for f in \
  "$STAGE1" \
  "$MUSCLE_REPO/scripts/muscle_stage2a_gtex.py" \
  "$MUSCLE_REPO/tests/test_muscle_stage2a_gtex.py"
do
  if [ -s "$f" ]; then
    echo "FOUND   $f"
  else
    echo "MISSING $f"
  fi
done

echo
echo "===== 5. API CONNECTIVITY ====="
python3 - <<'PY'
import json
import urllib.request

urls = {
    "ENSEMBL": "https://rest.ensembl.org/variation/human/rs74038095?content-type=application/json",
    "GTEX": "https://gtexportal.org/api/v2/dataset/variant?snpId=rs74038095&datasetId=gtex_v10&pageSize=5",
}

for name, url in urls.items():
    try:
        req = urllib.request.Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "IS_Analysis_V3-MUSCLE-Stage2A/1.0",
            },
        )
        with urllib.request.urlopen(req, timeout=30) as r:
            raw = r.read()
        json.loads(raw.decode("utf-8"))
        print(name, "OK", len(raw), "bytes")
    except Exception as e:
        print(name, "ERROR", repr(e))
PY

echo
echo "===== 6. COMPILE ====="
python3 -m py_compile "$MUSCLE_REPO/scripts/muscle_stage2a_gtex.py"
COMPILE_RC=$?
echo "COMPILE_RC=$COMPILE_RC"

echo
echo "===== 7. UNIT TEST ====="
python3 "$MUSCLE_REPO/tests/test_muscle_stage2a_gtex.py"
TEST_RC=$?
echo "TEST_RC=$TEST_RC"

echo
echo "===== 8. RUN STAGE 2A ====="
python3 \
  "$MUSCLE_REPO/scripts/muscle_stage2a_gtex.py" \
  --input "$STAGE1" \
  --outdir "$OUT" \
  --dataset gtex_v10 \
  --tissue Muscle_Skeletal \
  --pause 0.20 \
  2>&1 | tee "$RUNLOG"

RUN_RC=${PIPESTATUS[0]}
echo "RUN_RC=$RUN_RC"

echo
echo "===== 9. OUTPUT CHECK ====="
for f in \
  MUSCLE_STAGE2A_VARIANT_NORMALIZATION.tsv \
  MUSCLE_STAGE2A_EQTL.tsv \
  MUSCLE_STAGE2A_SQTL.tsv \
  MUSCLE_STAGE2A_CANDIDATE_GENES.tsv \
  MUSCLE_STAGE2A_SUMMARY.json
do
  if [ -s "$OUT/$f" ]; then
    echo "OK      $f ($(du -h "$OUT/$f" | cut -f1))"
  else
    echo "MISSING/EMPTY $f"
  fi
done

echo
echo "===== 10. SUMMARY ====="
python3 - <<'PY'
import json
from pathlib import Path

p = Path("/srv/is-analysis/results/muscle/stage2a_gtex/MUSCLE_STAGE2A_SUMMARY.json")
if not p.exists():
    print("SUMMARY_MISSING")
    raise SystemExit(0)

d = json.loads(p.read_text())
for k in [
    "version",
    "dataset",
    "tissue",
    "n_input_variants",
    "n_ensembl_grch38_resolved",
    "n_source_position_matches_grch38",
    "n_gtex_variant_resolved",
    "n_variants_with_significant_eqtl",
    "n_variants_with_significant_sqtl",
    "n_eqtl_associations",
    "n_sqtl_associations",
    "n_candidate_genes",
    "n_api_error_variants",
]:
    print(k, "=", d.get(k))

print()
print("TOP CANDIDATES")
for x in d.get("top_candidate_genes", [])[:20]:
    print(
        x.get("gene"),
        "score=" + str(x.get("score")),
        "rsids=" + str(x.get("rsids")),
        "eQTL=" + str(x.get("eqtl")),
        "sQTL=" + str(x.get("sqtl")),
        sep="\t",
    )

print()
if d.get("n_input_variants") == 20 and d.get("n_api_error_variants") == 0:
    print("STAGE2A_BASIC_QC=PASS")
else:
    print("STAGE2A_BASIC_QC=CHECK_REQUIRED")
PY

echo
echo "===== 11. TOP CANDIDATE GENE TABLE ====="
if [ -s "$OUT/MUSCLE_STAGE2A_CANDIDATE_GENES.tsv" ]; then
  head -21 "$OUT/MUSCLE_STAGE2A_CANDIDATE_GENES.tsv" | column -t -s $'\t'
fi

echo
echo "============================================================"
echo " MUSCLE STAGE 2A FINISHED"
echo "============================================================"
echo "WORKTREE=$MUSCLE_REPO"
echo "OUT=$OUT"
echo "LOG=$RUNLOG"
echo "FETCH_RC=$FETCH_RC"
echo "WORKTREE_RC=$WT_RC"
echo "COMPILE_RC=$COMPILE_RC"
echo "TEST_RC=$TEST_RC"
echo "RUN_RC=$RUN_RC"
echo "============================================================"
