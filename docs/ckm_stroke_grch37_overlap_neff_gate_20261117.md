# CKM × Stroke — Independent GRCh37 allele reference, sample overlap, and LDSC N gate

**2026-11-17 KST.** This document follows the original-data SUCCESS evidence from [GitHub Actions #43889295894](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43889295894), the verified aggregate review [#43895387493](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43895387493), and the BMI 1,007,899-variant three-way panel.

## 1. Orthogonal reference test — Ensembl GRCh37 FASTA, not a GWAS source

Official independent genomic sequence source: [Ensembl GRCh37 REST API: sequence/region](https://grch37.rest.ensembl.org/documentation/info/sequence_region_post). The API accepts batches of at most 50 GRCh37 genomic intervals and returns the forward-strand 1-bp reference nucleotide.

The local Orca test **actually completed** by streaming the complete pinned corrected original GIANT BMI source again and validating the original 1,007,899-variant candidate gzip/plaintext SHA-256. We selected **the two smallest SHA256(`ckm-stroke-grch37-ref-forward-v1|rsID`) values per autosome**, yielding **44 deterministic, non-significance-selected rsIDs**, 2 on each of chromosomes 1–22. Only these 44 public chr:bp:allele loci were submitted to Ensembl GRCh37 in a single POST request. The complete original BMI GWAS was not retained, and no GIGASTROKE original outcome was downloaded.

**Orca observed result: 44/44 GRCh37 forward-reference nucleotides belong to the source BMI allele pair**, 0 reverse-complement-only, 0 incompatible, 0 unresolved. Status `PASS_STRATIFIED_44_SENTINELS_ONLY`. 44/44 includes both cases in which the reference nucleotide equals GIANT Tested_Allele and equals GIANT Other_Allele; this is not an effect-direction association test.

- Code: `scripts/ckm_stroke/audit_bmi_grch37_reference_sentinels.py`.
- Offline regression tests: `tests/test_ckm_stroke_bmi_grch37_reference_sentinels.py` — 8/8 PASS on Orca.
- Local gitignored evidence: `artifacts/ckm_stroke/bmi_grch37_independent_sentinels_20261117.json`.
- Reproducible `GitHub Actions`: `.github/workflows/ckm-stroke-bmi-independent-grch37-alleles.yml`. Lightweight offline tests on code changes; only explicit `[bmi-ensembl44-full]` commit marker or manual dispatch retrieves reference sources and queries the public Ensembl API. Full original GWAS and temporary candidate gzip are never uploaded; compact JSON evidence only.

**Important limits:** This is *sampled* independent forward-reference evidence and cannot certify all 1,007,899 candidates, all 4 outcome allele directions, GWAS allele sign, strand harmonization on unsampled loci, or biological effects. The previous complete BMI↔GIGASTROKE join found all observed allele pairs `same` (AIS 1,006,932; LAA 956,929; CE 996,662; SVD 967,356). This unusually uniform orientation remains **REVIEW**, plausibly due to shared input processing/reference orientation.

## 2. Sample overlap: clear risk, unknown actual shared count

**Verified publisher evidence:**
- [Yengo et al. (2018), BMI GWAS](https://doi.org/10.1093/hmg/ddy271) combined previous GIANT cohorts with approximately **450,000 UK Biobank European participants**, reaching approximately 700,000 total participants. The exact source original is the GIANT consortium's corrected 2018 BMI summary-statistic release.
- [Mishra et al. (2022), GIGASTROKE](https://www.nature.com/articles/s41586-022-05165-3) explicitly discusses **limited overlap between European GIGASTROKE samples and exposure GWAS cohorts** and implemented sensitivity analyses excluding overlapping cases from UK Biobank exposure instruments. This statement is not proof of zero UKB overlap for the four exact publisher GRCh37 summaries.

**Risk status = OVERLAP_POSSIBLE_UNQUANTIFIED**. Do not assert independence or assign an overlap count from the number of cases, sample size or case/control totals. GIGASTROKE cohort lists and per-variant/phenotype cohort overlap should be compared to the exact GIANT/UKB exposure sample; if impossible, LDSC rg intercept and constrained/unconstrained intercept sensitivities must be presented as assumption-dependent and methodologically caveated, not interpreted as empirical proof of zero overlap.

## 3. Case-control N and effective N — numerical scenario, not source N per SNP

The outcome GWAS Catalog `buildGRCh37.tsv.gz` has **no N field**. Reported total case/control sample counts below are *aggregated phenotype-level counts* and can include cohort-specific missingness and varying SNP coverage; do not silently treat these as per-variant N.

| EUR GIGASTROKE | Case count | Control count | Total N | Illustrative `4*cases*controls/(cases+controls)` |
|---|---:|---:|---:|---:|
| AIS | 62,100 | 1,234,808 | 1,296,908 | ~236,506 |
| LAA | 6,399 | 1,234,808 | 1,241,207 | ~25,464 |
| CE | 10,804 | 1,234,808 | 1,245,612 | ~42,841 |
| SVD | 6,811 | 1,234,808 | 1,241,619 | ~27,095 |

The formula gives **illustrative effective sample sizes for a single homogeneous case-control study**, not a validated GIGASTROKE meta-analysis effective N. Analysis needs ascertainment/phenotype-scale specification, per-cohort mix, effective N convention for LDSC/munge, and documented methods from the outcome meta-analysis. Do not use raw total N uncritically in LDSC for highly imbalanced case/control data.

### LDSC methodological gate (explicit, fail-closed)

- **Technical input QC**: GIANT full original SHA256 and duplicated/P0 handling PASS; 1,007,899 BMI candidate vs EUR LDscore chr:BP PASS; 4 original GIGASTROKE EUR raw MD5 PASS; full BMI vs four outcome coordinate/allele category PASS; Ensembl 44 sentinel forward-reference PASS.
- **No automatic permission** to generate or run final LDSC h2, rg, constrained intercept, or trait causal MR based solely on technical QC.
- **Unresolved scientific gates**: actual cohort overlap estimates, correct N_eff convention and per-variant heterogeneity, outcome effect size/SE and binary trait scale, independent larger-sample genome/reference orientation QC, LD Score ancestry/build compatibility, h2 regression and intercept diagnostic criteria.
- Once those are addressed, run **single-trait h² and intercept diagnostics first**, then only consider pairwise rg if h²/SE/intercept diagnostics are satisfactory; do not infer causality from rg.
- **Scientific decision: HOLD**, `ready_for_ldsc=false`, no h2/rg/MR/coloc performed. Existing PR #284 remains Draft and unmerged.

## 4. Reproduction in Orca

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_bmi_grch37_reference_sentinels.py -v
py -3 -m scripts.ckm_stroke.audit_bmi_grch37_reference_sentinels --candidate artifacts/ckm_stroke/bmi_3way_conservative.sumstats.gz --provenance artifacts/ckm_stroke/bmi_ldsc_candidate_qc_20261116.json --out artifacts/ckm_stroke/bmi_grch37_independent_sentinels_20261117.json
```

The local report is not itself a direct validation of the GIGASTROKE allele sign: it checks the **BMI allele pair versus an independent genomic reference base at 44 predetermined loci only**.
