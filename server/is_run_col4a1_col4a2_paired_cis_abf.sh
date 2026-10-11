#!/usr/bin/env bash
# Reproducible, read-only COL4A1 vs COL4A2 paired SNP coloc.abf sensitivity.
# Usage: bash server/is_run_col4a1_col4a2_paired_cis_abf.sh 20261012_v2
# No canonical GWAS/QTL inputs or original 646 coloc results are modified.
set -euo pipefail
if [ "$#" -eq 0 ];then ID=20261011_v2;else ID="$1";fi
if ! [[ "$ID" =~ ^[A-Za-z0-9_-]{3,50}$ ]];then
 echo "INVALID_RUN_ID" >&2;exit 2
fi
ROOT=/srv/is-analysis
REPO="$ROOT/worktrees/is-broad-discovery-20261009"
A="$ROOT/results/is/audits"
INPUT="$A/is_col4a1_col4a2_paired_coloc_$ID"
ABF="$A/is_col4a1_col4a2_paired_abf_sensitivity_$ID"
SUM="$A/is_col4a1_col4a2_paired_abf_science_gate_$ID"
for p in "$INPUT" "$ABF" "$SUM";do
 if [ -e "$p" ];then
  echo "REFUSE_ANY_EXISTING_OUTPUT $p" >&2;exit 3
 fi
done
python3 "$REPO/scripts/is/prepare_is_col4a1_col4a2_paired_cis_coloc.py" \
 --qtl-root "$ROOT/data/is/qtl/eqtl_catalogue/eas5_T1_extended_cis_20261011_v1" \
 --bbj "$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz" \
 --ais "$ROOT/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz" \
 --fasta "$ROOT/reference/grch37_1000g/human_g1k_v37.fasta" \
 --chain "$ROOT/data/reference/hg38ToHg19.over.chain" \
 --vendor "$ROOT/data/is/reference/vendor_pyliftover" \
 --prior-98-file "$A/is_eas5_98_egene_fullcis_gwas_overlap_20261011_v3_fastaanchored/IS_EAS5_98_BH05_FDR_EGENE_FULLCIS_GWAS_ALLELE_READINESS.tsv" \
 --out-dir "$INPUT"
Rscript --vanilla "$REPO/scripts/is/run_is_col4a1_col4a2_paired_abf_sensitivity.R" \
 --input "$INPUT/IS_COL4A1_COL4A2_TISSUE_PAIRED_COMMON_SNP_INPUT.tsv" \
 --source-qc "$INPUT/IS_COL4A1_COL4A2_THREE_TISSUE_SNP_QC.tsv" \
 --out-dir "$ABF" \
 > "/tmp/is_col4a1_col4a2_paired_abf_$ID.log" 2>&1
python3 "$REPO/scripts/is/summarize_is_col4a1_col4a2_paired_abf.py" \
 --abf "$ABF/IS_COL4A1_COL4A2_THREE_QTL_P12_AF_N_ABF_SENSITIVITY.tsv" \
 --paired-input "$INPUT/IS_COL4A1_COL4A2_TISSUE_PAIRED_COMMON_SNP_INPUT.tsv" \
 --source-fdr "$A/is_eas5_studywide_BH_egenes_20261011_v1/IS_EAS5_720_OFFICIAL_STUDYWIDE_BH_EGENE_SOURCE_TESTS.tsv" \
 --legacy "$A/legacy_coloc_p12_sensitivity_20261010_v1/IS_LEGACY_646_ABF_P12_GRID.tsv" \
 --out-dir "$SUM"
Rscript --vanilla "$REPO/scripts/is/render_is_col4a1_col4a2_abf_prior_grid.R" \
 --input "$ABF/IS_COL4A1_COL4A2_THREE_QTL_P12_AF_N_ABF_SENSITIVITY.tsv" \
 --out "$ABF/IS_COL4A1_COL4A2_THREE_TISSUE_H3_H4_PRIOR_SENSITIVITY.svg"
Rscript --vanilla "$REPO/scripts/is/render_is_col4a1_col4a2_locus_associations.R" \
 --input "$INPUT/IS_COL4A1_COL4A2_TISSUE_PAIRED_COMMON_SNP_INPUT.tsv" \
 --out "$ABF/IS_COL4A1_COL4A2_ROSMAP_EAS_AIS_GWAS_CIS_LOCUS_PLOT.svg"
echo "IS_COL4A1_COL4A2_REPRODUCIBLE_PAIR_SOURCE_COLOC_P12_QC_COMPLETE"
echo "INPUT=$INPUT ABF=$ABF SUMMARY=$SUM"
echo "STATUS=EXPLORATORY_NOT_UNIQUE_GENE_CAUSALITY"
