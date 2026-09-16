# CKD Stage 1: proteome-wide MR screen

## Scope

Stage 1 screens UKB-PPP Sun et al. 2023 cis-pQTL instruments against the compact CKD outcome set.

Primary discovery outcome:
- EUR eGFRcrea

Supporting EUR phenotypes:
- CKD
- BUN
- eGFRcys
- UACR

Cross-ancestry support:
- EAS eGFRcrea
- EAS BUN

## Exposure orientation

Sun et al. Supplementary Table 16 encodes variants as `CHR:POS:REF:ALT` and reports
`ALT freq`. The conditional exposure estimate `Beta(cond)` is treated as the ALT
allele effect. This is consistent with the published rs4985556 example: the A allele
is associated with lower IL34, while ST16 encodes rs4985556 as C>A with a negative
conditional beta.

Protein measurements in the UKB-PPP analysis were inverse-rank normalized before
association testing. MR estimates should therefore be interpreted as effects per
genetically predicted increase in normalized plasma protein level, not raw NPX units.

## Instrument QC

`F = beta_cond^2 / se_cond^2`

Default threshold:
- F >= 10

ST16 already contains SuSiE/conditional independent signals, but multiple signals
within one cis region may retain residual LD. Therefore:

- Primary screen: strongest-F harmonized cis instrument per protein, Wald ratio.
- Sensitivity: fixed-effect IVW across all available ST16 signals.
- Do not treat IVW sensitivity as definitive until ancestry-matched LD/covariance
  is explicitly checked.

## Harmonization

Exposure effect allele:
- ALT

Outcome beta is reoriented to ALT.

Rules:
- direct allele match: keep beta
- swapped allele match: flip beta
- non-palindromic strand complement: keep/flip as appropriate
- palindromic A/T or C/G:
  - require exposure ALT frequency and outcome EAF
  - drop if exposure MAF > 0.42
  - use frequency concordance within 0.10 to choose orientation
  - otherwise drop
- EAS compact files currently do not contain EAF, so palindromic variants are
  conservatively excluded there.

## Multiple testing

EUR eGFRcrea is the only Stage 1 discovery endpoint.

Per-protein strongest-instrument Wald P values are adjusted with Benjamini-Hochberg
FDR across tested proteins. Bonferroni status is also emitted.

CKD/BUN/eGFRcys/UACR are supporting phenotypes and are not counted as additional
discovery endpoints.

## Run on server

```bash
cd /srv/is-analysis/IS_Analysis_V3
git pull --ff-only origin main

CKD_READY_ROOT=/srv/is-analysis/data/ckd/analysis_ready \
CKD_STAGE1_ROOT=/srv/is-analysis/results/ckd/stage1 \
SYNC_DRIVE=1 \
RCLONE_REMOTE=gdrive \
DRIVE_RESULTS_BASE='IS_Analysis_V3/results/ckd/stage1' \
  server/ckd_run_stage1_mr.sh
```

Expected terminal markers:

```text
CKD_STAGE1_PASS
CKD_STAGE1_LOCAL_PASS
CKD_STAGE1_DRIVE_SYNC_PASS
```

## Main outputs

```text
stage1/
├── instrument_qc.tsv.gz
├── harmonized/
│   ├── EUR/
│   └── EAS/
├── mr/
│   ├── EUR/
│   └── EAS/
├── stage1_protein_summary.tsv
├── stage1_eur_egfr_fdr05_hits.tsv
├── STAGE1_SUMMARY.json
└── SHA256SUMS.txt
```

## Interpretation

Stage 1 signals are screening results, not causal claims.

Required downstream filters:
1. phenotype consistency
2. EAS replication / variant recovery
3. reverse MR for kidney function -> plasma protein
4. candidate-specific full cis-region retrieval
5. colocalization
6. kidney pQTL/eQTL and kidney cell-type localization
7. drug-target triangulation

Kidney dysfunction can alter circulating protein concentrations through filtration,
clearance and systemic effects. Reverse-direction analysis is therefore essential.
