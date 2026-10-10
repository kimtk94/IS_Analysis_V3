# IS existing-data-only progress: 2,225-gene matrix and 646/646 SNP-level coloc recovery

**Date:** 2026-10-11 (KST). **No new controlled/restricted/public GWAS database was acquired.** Work uses previously saved server evidence and the user's pre-existing Google Drive `MASTER_DEGREE/IS_COLAB/coloc_ready_v2/abf` backup.

## Completed work: 80-region, 2,225-gene evidence matrix

- **80 provisional regions**, comprising 30 Japanese/EAS discovery components and 50 EUR positional clusters, **not** 80 statistically independent causal loci.
- **2,225 unique stable Ensembl gene IDs**, **2,425** distinct positional gene–region records: **869 EAS/Japanese**, **1,556 EUR**.
- One `NEURL1` legacy-anchor-only row in the expanded input is held in a separate ledger. It maps by existing GTEx `gene_base=ENSG00000107954` to the GENCODE v19 positional symbol **`NEURL`**, not to an invented 2,226th candidate. All eight previously prioritized gene–tissue assays attach through exact stable IDs and BBJ locus.
- The 646 archived coloc assays were associated with 43 EAS positional gene IDs; none were transferred to an EUR-only region as independent QTL evidence.
- Human healthy-reference donor data (GSE256493, **80,515 cells, six donors**) are descriptive expression localization, **not ischemic-stroke disease differential expression**.
- Updated readiness tier from the **646 new direct ABF replays**:
  - **43** genes with historic EAS GTEx SNP-level ABF repro contexts, **not** proven causal genes;
  - **1,049** genes prioritized for positional GWAS region follow-up;
  - **1,133** genes retained as other positional/suggestive follow-up.
  - **2,225 total; none excluded or declared causally validated.**

Primary derivative:
`/srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_646replayed_20261011_v2/`

Files:
- `IS_2225_GENE_EVIDENCE_MATRIX.tsv`;
- `IS_2425_GENE_REGION_EVIDENCE.tsv`;
- `IS_LEGACY_ONLY_ANCHORS.tsv`;
- `IS_GENE_EVIDENCE_MATRIX_MANIFEST.json` with immutable source SHA-256 references.

## Original 646 GWAS/QTL per-SNP input files: backup recovered and verified

Previous on-server SNP replay had **30 PASS**, **616 MISSING_INPUT**. A separate Google Drive folder already held exactly **646** correctly named assay input files. An authorized `rclone` copy populated a **new** derivative directory:

`/srv/is-analysis/results/is/audits/legacy_coloc_drive_recovery_20261011_v1/input/`

All **646 files (669,756,087 bytes combined)** were present after copying. Remote-original versus local backup **rclone MD5 checksum check: 646 matching files, zero differences**. No original 30 archived inputs, historic master, or other canonical paths were overwritten. A separate file-by-file source schema/row/allele audit returned:
- **30 `LOCAL_SOURCE_QC_PASS`** (originally available on server);
- **616 `DRIVE_RECOVERED_QC_PASS`** (recovered from existing Drive);
- **0 local input QC conflicts**.

Audit:
`/srv/is-analysis/results/is/audits/is_646_drive_source_final_20261011_v1/`

Note: source file recovery by itself is **not** independent statistical or biological validation.

## Direct per-SNP `coloc.abf` recomputation: 646/646 PASS

Using recovered SNP-level `gwas_beta, gwas_se, gwas_maf, eqtl_beta, eqtl_se, eqtl_maf, match_key` and **coloc v5.2.3**:

- Archived priors: `p1=1e-4, p2=1e-4, p12=1e-5`.
- BBJ IS binary-outcome model: `N=174,686`, `cases=22,664` (`s=22,664/174,686`).
- GTEx v8 quantitative QTL model: historically used the **first SNP's scalar `qtl_n_scalar`** for each assay.
- **646/646 all five H0–H4 posteriors reproduced to max absolute error 6.71685e-15**, far below the `1e-8` pass threshold.
- **Zero / 646 archived baseline assays had `PP.H4 >= 0.8`** at these pre-specified priors. Prior choice materially affects H4; do not select favorable priors post hoc.
- No biological mechanism or causal gene is validated merely by numerical reproduction.

Main derivative:
`/srv/is-analysis/results/is/audits/is_646_coloc_replay_eqtl_per_snp_n_20261011_v1/IS_646_GWAS_QTL_SNP_ABF_REPLAY.tsv`

Scientific QC:
`/srv/is-analysis/results/is/audits/is_646_coloc_replay_science_qc_per_snp_n_20261011_v1/IS_646_COLOC_REPLAY_SCIENCE_READINESS.json`

## Major scientific QC flags requiring follow-up

### GWAS vs GTEx molecular-QTL MAF discordance

Absolute `|GWAS MAF - QTL MAF| > 0.1` in:
- **More than 25% and at most 50% of SNPs in 588/646 assays**;
- **More than 50% of SNPs in 58/646 assays**.

The **median per-assay fraction was 0.414403**, across all 646 assays. This reflects substantial cross-cohort/ancestry frequency differences (with possible remaining harmonization or QC issues), not proof of an allele flip or proof of artifact. Both source phenotype/ancestry assumptions and palindromic/allele-frequency handling deserve a pre-specified independent review. **Do not retroactively remove variants only to raise H4.**

### Per-SNP GTEx sample N varies; historic model uses first SNP N

**All 646 assays** have variable per-SNP `eqtl_n`; median **8 distinct N values per assay**, maximum **16**. The former readout that suggested `0` was an instrumentation error caused by examining the fixed `qtl_n_scalar` rather than `eqtl_n`, and has been corrected.

### Full 646-assay QTL scalar-N sensitivity

Using the same SNP-level summary data under four explicit **scalar** GTEx N values per assay: historical first SNP, minimum, rounded median, and maximum (`646 * 4 = 2,584` coloc model runs):

- **Max absolute shift in H4** vs historical first-N: **0.0429922895**.
- **Median assay maximum absolute shift:** approximately **0.0014**.
- **0/646 assays shifted H4 by more than 0.05**.
- **1/646 crossed the descriptive H4=0.5 threshold**; **0/646 crossed H4=0.8**.
- The crossing assay was **INA (`ENSG00000148798`)**, `BBJ_IS_L002`, `GTEx_V8__Brain_Anterior_cingulate_cortex_BA24`: `H4_first = 0.53636854`, `H4_min_QTL_N = 0.49337625`.
- **CALHM2** `BBJ_IS_L002`, cerebellar hemisphere: `H4_first = 0.79768355`, `H4_min_QTL_N = 0.76664289`; neither satisfies `H4>=0.8` under this test.

These are *alternative scalar-N assumptions*, **not actual SNP-specific N models**, and they do **not** address LD heterogeneity, multi-signal assumptions, or frequency discordance.

Derivative:
`/srv/is-analysis/results/is/audits/is_646_qtl_scalar_n_sensitivity_20261011_v1/`:
- `IS_646_COLOC_QTL_N_SCALAR_SENSITIVITY.tsv`;
- `IS_646_COLOC_QTL_N_SENSITIVITY_MANIFEST.json`.

## Software, tests and inference gate

- Source code on branch `research/is-broad-discovery-20261009`:
  - `scripts/is/build_existing_is_gene_evidence_matrix.py`
  - `scripts/is/audit_local_drive_coloc_abf_recovery.py`
  - `scripts/is/replay_existing_646_coloc_abf.R`
  - `scripts/is/summarize_646_coloc_replay_science_readiness.py`
  - `scripts/is/run_existing_646_coloc_qtl_n_sensitivity.R`
- **Ten additional offline Python regression tests passed on the server** for gene ID joins, input integrity and posterior QC. Direct R coloc reproduction and scalar-N sensitivity were executed with complete 646-assay input. The CI workflow has been updated, but **remote GitHub Actions run status was not separately verified**.
- QTL matched cohort LD is not available; ABF's single-causal-signal assumption has not been invalidated or validated through multi-signal matched-LD fine-mapping. `PP.H4` is a posterior under those assumptions; none of 646 deserves `VERIFIED_CAUSAL_GENE`.
- GIGASTROKE EAS contains BBJ; BBJ and GIGASTROKE EAS do not form independent replication. This pipeline is a GTEx expression-QTL versus BBJ ischemic-stroke coloc branch, **not an alcohol-consumption MR proof**.
- **ALDH2 alcohol fine-mapping remains BLOCKED, ADH1B alcohol fine-mapping EXPLORATORY**.
- Original 2,225 positional genes are retained.

## Next local-only actions

1. Design a **pre-specified allele-frequency harmonization audit** using matched effect-allele AF, palindromic exclusions, and ancestry-aware interpretation, preserving full unfiltered 646-assay results in parallel. Significant MAF differences warrant scientific review, not automated evidence deletion.
2. Analyze **646 coloc H3 vs H4** across tissues, gene IDs and locus contexts; add prior grids for all 646 using direct SNP-level analysis and avoid favorable-prior selection. Explicitly distinguish H4 model support from signal colocalization under multi-signal LD.
3. Produce publication-ready Figure-ready outputs: full 2,225-gene coverage heatmap, per-assay MAF difference vs H4, prior sensitivity, 6-donor healthy reference expression and competing mechanisms. Preserve uncertainty flags in all figures.
4. Reassess gene mechanisms and rank **follow-up readiness**, not causality; keep all 2,225 positional gene IDs.
