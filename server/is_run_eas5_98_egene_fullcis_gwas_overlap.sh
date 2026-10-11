#!/usr/bin/env bash
# Read-only offline analysis. Requires previously SHA256/EBI-checked source files.
# No trading/service changes. Non-idempotent derivative output is never overwritten.
set -euo pipefail
ROOT=/srv/is-analysis
SCRIPT="$ROOT/worktrees/is-broad-discovery-20261009/scripts/is/audit_is_eas5_98_egene_fullcis_gwas_exact.py"
OUT="$ROOT/results/is/audits/is_eas5_98_egene_fullcis_gwas_overlap_20261011_v2"
LOG=/tmp/is_eas5_98_fullcis_gwas_overlap_20261011_v2.log
if [ -e "$OUT" ];then
 echo "REFUSE_OVERWRITE_EXISTING_SCIENCE_OUTPUT $OUT" >&2;exit 1
fi
python3 "$SCRIPT" \
 --fdr "$ROOT/results/is/audits/is_eas5_studywide_BH_egenes_20261011_v1/IS_EAS5_720_OFFICIAL_STUDYWIDE_BH_EGENE_SOURCE_TESTS.tsv" \
 --master "$ROOT/results/is/audits/is_eas5_full_cis_master_readiness_20261011_v1/IS_EAS5_720_QTL_FULLCIS_BEFORE_AFTER_READINESS.tsv" \
 --older "$ROOT/results/is/audits/is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1/IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv" \
 --bbj "$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz" \
 --ais "$ROOT/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz" \
 --chain "$ROOT/data/reference/hg38ToHg19.over.chain" \
 --fasta "$ROOT/reference/grch37_1000g/human_g1k_v37.fasta" \
 --vendor "$ROOT/data/is/reference/vendor_pyliftover" \
 --qtl-g0015 "$ROOT/data/is/qtl/eqtl_catalogue/eas5_G0015_extended_cis_20261011_v1" \
 --qtl-t1 "$ROOT/data/is/qtl/eqtl_catalogue/eas5_T1_extended_cis_20261011_v1" \
 --qtl-t4 "$ROOT/data/is/qtl/eqtl_catalogue/eas5_T4_extended_cis_20261011_v1" \
 --qtl-g0018 "$ROOT/data/is/qtl/eqtl_catalogue/eas5_G0018_edge_cis_20261011_v1" \
 --out-dir "$OUT" > "$LOG" 2>&1 || { tail -n 30 "$LOG"; exit 1; }
echo "PASS $OUT"
tail -n 3 "$LOG"
