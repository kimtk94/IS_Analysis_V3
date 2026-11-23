# IS P3A — ROSMAP 15-gene original QTL coloc sdY estimator and N sensitivity

**2026-11-23 KST.** Ischemic stroke only. 36 original gene×QTL tissue assays (15 ROSMAP gene expression QTD000434 and 21 GTEx), 20 distinct genes, one chr10 stroke locus. Formal coloc.abf/susie remains **0/36 HOLD**.

## What is being verified

Previous merged [PR #323](https://github.com/kimtk94/IS_Analysis_V3/pull/323) recovered the exact original publisher QTD000434 SNP associations for **all 15 ROSMAP candidate genes**: 367,753 original chr10 rows in the audit union, 100,822 for 15 genes, 89,844 gene-wise deduplicated records. Every observed candidate QTL SNP has **AN=1120**, and the official eQTL Catalogue says the autosomal sample-size proxy is AN/2 = 560. This number is NOT individually validated regression N.

The original [coloc R/claudia.R](https://github.com/chr1swallace/coloc/blob/main/R/claudia.R) (Git blob SHA1 043783ac5ddec3477ff3b91e3e1faf271c51ce09) implements sdY.est(vbeta, maf, n): for each SNP set varbeta = SE squared, oneover = 1/varbeta, nvx = 2*N*maf*(1-maf), regress nvx on oneover WITHOUT INTERCEPT and return the square root of its slope. The original coloc package emits a warning because this is an ESTIMATED trait SD, not a directly observed normalized phenotype SD. [Publisher reference](https://chr1swallace.github.io/coloc/reference/sdY.est.html).

Original eQTL Catalogue [batch coloc implementation](https://github.com/eQTL-Catalogue/colocalisation/blob/master/bin/run_coloc_batch.R) uses eQTL beta, se squared, maf and first variant AN/2 as quantitative QTL N. This published example is not validation of exact gene-wise regression N for Release7 QTD000434 or binary stroke GWAS case fraction.

## Original publisher source integrity

Fetch official Ensembl105 phenotype metadata ZIP/gzip (Zenodo pinned MD5 b999b21dcffeecbfe0bf886a5a91997b), the **single public QTL index** (SHA256 6cbf12826e19a61594662972a1205ba38cdd3071984e51c24261bff8ace9806a), and one original QTD000434 BGZF range (SHA256 5a54c1faadb0530ecdca37497f82f2e28f1003c000d7a24bb9fa8422de3fa2f1). Never download the full ~4GB source, never repeat 15 separate remote queries.

Require each original gene and GRCh38 REF/ALT, numeric beta/SE/MAF/AN, original association count and all 15 source per-gene unique variant counts exactly matching PR #323. Source expected gene IDs and variant counts: SLK 5854, GSTO2 6093, NEURL 5961, CCDC147 6213, CALHM2 6018, INA 6031, GSTO1 6103, ITPRIP 6202, CNNM2 5809, C10orf32 5728, USMG5 6041, WDR96 6099, AS3MT 5723, RP11-127L20.3 6163, RP11-724N1.1 5806.

## Diagnostics and sensitivity

Calculate per gene:
1. Full real source cis variant estimate of sdY, identical to the R zero-intercept OLS estimator, **under N=AN/2=560 assumption**.
2. Descriptive regression R squared, median/P10/P90 single-SNP implied SD.
3. Independent estimations for **MAF below 0.10** and **MAF at least 0.10** where at least 30 variants per stratum; ratio and reported counts, not a permission threshold.
4. Odd/even splits after sorting variant IDs. These are **not statistically independent**, because of LD and genomic position; purely a numerical stability indicator.
5. Hypothetical source-N sensitivity at **540 / 560 / 580**, with exact analytical scaling sdY(N)/sdY(560) = sqrt(N/560). The two alternative N values are intentionally NOT claimed to be real sample counts and do not permit effective-N substitution.

**Interpretation:** Original expression SD remains unmeasured without the normalized QTL input phenotype. Model residual, MAF subgroup variation and splitting are diagnostics. The 2023 phenotype position and ±1Mb default are not independently verified Release7 runtime windows. Model estimates are not automatically approved for formal coloc.

## Files and review gates

- Python: scripts/is/audit_p3a_rosmap_real_sdY_est_15gene.py
- Tests: tests/test_is_p3a_rosmap_real_sdY_est_15gene.py
- Actions: .github/workflows/is-p3a-rosmap-15gene-sdy-est-sensitivity.yml
- Artifact: only 15 gene aggregate numerical diagnostics, source hashes, 36 HOLD statuses. No per-SNP effect records or donor IDs.

Original GIGASTROKE EAS+EUR chr10 source beta = ln(OR) is directly verified on 725,707 records in [PR #325](https://github.com/kimtk94/IS_Analysis_V3/pull/325). **Next scientific gate:** exact source ancestry/study case fraction s and appropriate per-file sample N, QTL N and actual sdY assay normalization validation, and harmonization of GRCh38 QTL with GRCh37 GWAS. Do not invent s/N_eff, set sdY=1 or infer H4. Formal Coloc remains **0/36 HOLD**.
