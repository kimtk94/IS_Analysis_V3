#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STAGE1="${CKD_STAGE1_ROOT:-/srv/is-analysis/results/ckd/stage1}"
STAGE2="${CKD_STAGE2_ROOT:-/srv/is-analysis/results/ckd/stage2}"
STAGE2B="${CKD_STAGE2B_ROOT:-/srv/is-analysis/results/ckd/stage2b_coloc}"
STAGE2C="${CKD_STAGE2C_ROOT:-/srv/is-analysis/results/ckd/stage2c_susie}"
STAGE3A="${CKD_STAGE3A_ROOT:-/srv/is-analysis/results/ckd/stage3a_kidney}"
STAGE3B="${CKD_STAGE3B_ROOT:-/srv/is-analysis/results/ckd/stage3b_celltype}"
OUT="${CKD_MASTER_ROOT:-/srv/is-analysis/results/ckd/master}"
PYTHON="${CKD_PYTHON:-python3}"
RUN_STAGE4_SENSITIVITY="${RUN_STAGE4_SENSITIVITY:-0}"
RUN_PROVISIONAL_STAGE19="${RUN_PROVISIONAL_STAGE19:-0}"

mkdir -p "$OUT"

echo "============================================================"
echo " CKD MASTER — MATERIALIZE EXISTING EVIDENCE"
echo "============================================================"

test -f "$STAGE1/stage1_protein_summary.tsv"
test -f "$STAGE1/mr/EUR/eGFRcrea.tsv"
test -f "$STAGE1/mr/EAS/eGFRcrea.tsv"
test -f "$STAGE2/stage2_candidates.tsv"
test -f "$STAGE2B/coloc_results/STAGE2B_COLOC_DEFAULT.tsv"
test -f "$STAGE2C/susie_results/STAGE2C_SUSIE_DEFAULT.tsv"

echo "[Stage7] Cross-ancestry bridge"
"$PYTHON" "$ROOT/scripts/run_master_cross_ancestry.py" \
  --discovery-mr "$STAGE1/mr/EUR/eGFRcrea.tsv" \
  --replication-mr "$STAGE1/mr/EAS/eGFRcrea.tsv" \
  --discovery-coloc "$STAGE2B/coloc_results/STAGE2B_COLOC_DEFAULT.tsv" \
  --output "$OUT/stage07_cross_ancestry.tsv" \
  --discovery-ancestry EUR \
  --replication-ancestry EAS

echo "[Stage10] Phenotype expansion bridge"
PHENO_MANIFEST="$OUT/stage10_phenotype_manifest.tsv"
cat > "$PHENO_MANIFEST" <<EOF
mr_file	label	group	coloc_file
$STAGE1/mr/EUR/eGFRcrea.tsv	eGFRcrea	kidney_function	$STAGE2B/coloc_results/STAGE2B_COLOC_DEFAULT.tsv
$STAGE1/mr/EUR/CKD.tsv	CKD	kidney_disease	
$STAGE1/mr/EUR/BUN.tsv	BUN	kidney_function	
$STAGE1/mr/EUR/eGFRcys.tsv	eGFRcys	kidney_function	
$STAGE1/mr/EUR/UACR.tsv	UACR	kidney_damage	
EOF
"$PYTHON" "$ROOT/scripts/run_master_phenotype_matrix.py" \
  --manifest "$PHENO_MANIFEST" \
  --long-output "$OUT/stage10_phenotype_long.tsv" \
  --wide-output "$OUT/stage10_phenotype_matrix.tsv"

echo "[Stage14/15] Existing kidney localization bridge"
if [[ -f "$STAGE3A/STAGE3A_KIDNEY_EVIDENCE.tsv" && -f "$STAGE3B/STAGE3B_INTEGRATED_EVIDENCE.tsv" ]]; then
  "$PYTHON" "$ROOT/scripts/bridge_ckd_master_localization.py" \
    --stage3a "$STAGE3A/STAGE3A_KIDNEY_EVIDENCE.tsv" \
    --stage3b "$STAGE3B/STAGE3B_INTEGRATED_EVIDENCE.tsv" \
    --output "$OUT/stage14_15_localization.tsv"
else
  echo "WARN: Stage3A/3B localization source missing; Stage14/15 bridge not written." >&2
fi

if [[ "$RUN_STAGE4_SENSITIVITY" == "1" ]]; then
  echo "[Stage4] MR robustness sensitivity"
  test -f "$STAGE1/harmonized/EUR/eGFRcrea.tsv.gz"
  "$PYTHON" "$ROOT/scripts/run_master_mr_robustness.py" \
    --input "$STAGE1/harmonized/EUR/eGFRcrea.tsv.gz" \
    --output "$OUT/stage04_mr_robustness.tsv"
  cat > "$OUT/stage04_CAVEAT.txt" <<EOF
Stage4 robustness uses multiple UKB-PPP ST16 conditional cis signals.
Residual LD may remain among conditional signals. Treat weighted median, MR-Egger,
Cochran Q and leave-one-out as sensitivity analyses until ancestry-matched LD/covariance
is explicitly verified. Do not upgrade causal evidence solely from this file.
EOF
else
  echo "[Stage4] skipped by default. Set RUN_STAGE4_SENSITIVITY=1 to materialize with residual-LD caveat."
fi

echo "[Stage19] Evidence manifest template"
EVM="$OUT/stage19_evidence_manifest.tsv"
cat > "$EVM" <<EOF
stage	file	key_column
mr	$STAGE1/mr/EUR/eGFRcrea.tsv	gene_symbol
coloc	$STAGE2B/coloc_results/STAGE2B_COLOC_DEFAULT.tsv	gene_symbol
susie	$STAGE2C/susie_results/STAGE2C_SUSIE_DEFAULT.tsv	gene_symbol
ancestry	$OUT/stage07_cross_ancestry.tsv	gene_symbol
phenotype	$OUT/stage10_phenotype_matrix.tsv	gene_symbol
EOF
if [[ -f "$OUT/stage04_mr_robustness.tsv" ]]; then
  printf "robustness\t%s\tgene_symbol\n" "$OUT/stage04_mr_robustness.tsv" >> "$EVM"
fi
if [[ -f "$OUT/stage14_15_localization.tsv" ]]; then
  printf "localization\t%s\tgene_symbol\n" "$OUT/stage14_15_localization.tsv" >> "$EVM"
fi

CAND="$OUT/stage19_candidates.tsv"
awk -F "\t" 'NR==1 {for(i=1;i<=NF;i++) if($i=="gene_symbol") g=i; print "gene_symbol"; next} g && $g!="" && !seen[$g]++ {print $g}' \
  "$STAGE2/stage2_candidates.tsv" > "$CAND"

if [[ "$RUN_PROVISIONAL_STAGE19" == "1" ]]; then
  echo "[Stage19] PROVISIONAL integration from currently materialized evidence only"
  "$PYTHON" "$ROOT/scripts/run_master_evidence_integration.py" \
    --candidates "$CAND" \
    --manifest "$EVM" \
    --output "$OUT/stage19_PROVISIONAL_evidence_matrix.tsv"
  cat > "$OUT/stage19_PROVISIONAL_CAVEAT.txt" <<EOF
This evidence matrix is provisional. Missing Stage8/9/11/12/13/16/17/18 evidence
must not be scored as negative evidence. Do not use the provisional tier as the
manuscript final candidate classification.
EOF
else
  echo "[Stage19] provisional integration skipped. Set RUN_PROVISIONAL_STAGE19=1 only for pipeline QA."
fi

echo
echo "Materialized:"
find "$OUT" -maxdepth 1 -type f -name "stage*" -printf "%f\n" 2>/dev/null | sort || true
echo "CKD_MASTER_EXISTING_BRIDGE_PASS"
