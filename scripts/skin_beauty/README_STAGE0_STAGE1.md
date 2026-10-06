# SKIN_BEAUTY_MR — Stage 0/1 runbook

## Operational outcomes

- `ebi-a-GCST90094903`: Facial wrinkles (under eye), Korean/East Asian, OpenGWAS N=11,079, 4,939,316 variants.
- `ebi-a-GCST90094904`: Facial wrinkles (crow's feet), Korean/East Asian, OpenGWAS N=11,079, 4,939,316 variants.

The 2022 paper reports 17,019 Korean women across discovery + replication. OpenGWAS corresponds to the discovery GWAS (N=11,079).

## 1. Sync scripts from Google Drive

```bash
cd /srv/is-analysis/IS_Analysis_V3

mkdir -p scripts/skin_beauty
rclone copy \
  "gdrive:MASTER_DEGREE/SKIN_BEAUTY_MR" \
  scripts/skin_beauty \
  --include '*.py' \
  --include '*.R' \
  --include '*.sh' \
  --include '*.md' \
  --progress

chmod +x scripts/skin_beauty/*.sh scripts/skin_beauty/*.py scripts/skin_beauty/*.R 2>/dev/null || true
```

## 2. OpenGWAS authentication

Current OpenGWAS protected endpoints require a JWT. Generate one from the OpenGWAS profile page and keep it only in the shell environment.

```bash
export OPENGWAS_JWT='PASTE_TOKEN_HERE'
```

Do **not** save the JWT in Git, Drive, scripts, `.env` files that are committed, or logs.

## 3. Stage 0

```bash
cd /srv/is-analysis/IS_Analysis_V3/scripts/skin_beauty

ROOT=/srv/is-analysis/IS_Analysis_V3 \
  bash run_skin_stage0.sh
```

Expected outputs:

```text
results/skin_beauty/stage0_audit/STAGE0_SKIN_GWAS_AUDIT.json
results/skin_beauty/stage0_audit/STAGE0_EXPOSURE_FILE_CANDIDATES.tsv
results/skin_beauty/stage0_audit/STAGE0_OPENGWAS_METADATA.tsv
results/skin_beauty/stage0_audit/STAGE0_OPENGWAS_FILES.tsv
results/skin_beauty/stage0_audit/STAGE0_WRINKLE_TOPHITS.tsv
```

## 4. Stage 1 cis-pQTL MR

If the Stage 0 candidate list identifies the correct existing UKB-PPP cis-pQTL instrument file, set it explicitly:

```bash
export SKIN_PQTL_EXPOSURE='/srv/is-analysis/.../UKBPPP_...cis....tsv.gz'

cd /srv/is-analysis/IS_Analysis_V3/scripts/skin_beauty
ROOT=/srv/is-analysis/IS_Analysis_V3 \
  bash run_skin_stage1.sh
```

If `SKIN_PQTL_EXPOSURE` is omitted, the R script performs a best-effort search under the project root. Explicit selection is preferred once Stage 0 identifies the intended file.

Stage 1 queries only instrument rsIDs against each wrinkle outcome with `proxies=0`. This is deliberate: OpenGWAS proxy search is European-reference based, which should not be silently used for a Korean/East-Asian outcome.

Expected outputs:

```text
results/skin_beauty/stage1_mr/STAGE1_EXPOSURE_INSTRUMENTS.tsv
results/skin_beauty/stage1_mr/STAGE1_OUTCOME_UNDER_EYE.tsv
results/skin_beauty/stage1_mr/STAGE1_OUTCOME_CROWS_FEET.tsv
results/skin_beauty/stage1_mr/STAGE1_HARMONISED.tsv
results/skin_beauty/stage1_mr/STAGE1_MR_RESULTS.tsv
results/skin_beauty/stage1_mr/STAGE1_MR_RESULTS_WITH_MULTIPLE_TESTING.tsv
results/skin_beauty/stage1_mr/STAGE1_HETEROGENEITY.tsv
results/skin_beauty/stage1_mr/STAGE1_EGGER_INTERCEPT.tsv
```

## 5. Decision gate before Stage 2

Proceed to coloc/SuSiE only after confirming:

- exact exposure file and protein identifiers,
- F > 10 instruments,
- successful allele harmonisation,
- number of tested proteins and multiple-testing threshold,
- MR-positive / prioritized proteins,
- EAS LD reference availability for each candidate locus.

For Stage 2, query only regional windows around prioritized protein loci and use an ancestry-matched EAS LD matrix locally.
