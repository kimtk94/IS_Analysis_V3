#!/usr/bin/env bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
cd "$ROOT" || exit 1

set +e
set +u
set +o pipefail 2>/dev/null || true

R_LIB="$ROOT/.R-metabolic-resilience"
REG="$ROOT/results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv"
CIS="$ROOT/results/metabolic_resilience/stage3_full_pgwas/cis_marginal"
VCF="$ROOT/data/metabolic_resilience/stage2_gwas/ld_reference_1kg_eur/stage3_regional_vcf"

mkdir -p "$R_LIB"

echo "===================================================="
echo "STAGE 3-E0 : COLOCALIZATION PREFLIGHT"
echo "===================================================="

if [ -s "$REG" ]; then
  echo "[PASS] $REG"
else
  echo "[FAIL] $REG"
fi

echo
echo "cis files:"
find "$CIS" -maxdepth 1 -name '*cis1Mb_marginal.tsv.gz' -type f | wc -l

echo "regional VCFs:"
find "$VCF" -maxdepth 1 -name '*.vcf.gz' -type f | wc -l

echo
echo "R package audit:"

Rscript - "$R_LIB" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1], .libPaths()))
for (p in c("coloc","susieR")) {
  cat(p, requireNamespace(p, quietly=TRUE), "\n")
}
RS

if [ "${INSTALL_COLOC_PACKAGES:-0}" = "1" ]; then
  echo
  echo "Installing coloc + susieR into $R_LIB"
  Rscript - "$R_LIB" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
lib <- args[1]
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
need <- c("coloc","susieR")
missing <- need[!sapply(need, requireNamespace, quietly=TRUE)]
if (length(missing)) {
  install.packages(missing, repos="https://cloud.r-project.org", lib=lib)
}
RS
fi

echo
echo "NOTE:"
echo "- coloc.abf can proceed after regional outcome extraction."
echo "- T2D case-control coloc remains HOLD until no-UKB case proportion is verified."
echo "- Do not assume binary sample fraction."
