# CKM × Stroke — GIGASTROKE four EUR outcome full-source evidence review

**Reviewed:** 2026-11-17 KST · **Source of empirical results:** GitHub Actions [#43889295894](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43889295894), run completed SUCCESS at 2026-11-16 15:08:19 UTC (2026-11-17 00:08:19 KST). Exact ZIP artifact [#14237113287](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43889295894/artifacts/14237113287), containing `bmi_candidate_four_original_ci.json` and `gigastroke_four_full_md5_bmi_allele_ci.json`.

## Primary-source results: actual complete original gzip files

The audit streamed all 4 **official original GRCh37 EUR** GWAS Catalog files through a bounded gzip parser and checked each entire compressed stream's computed MD5 against its corresponding official GWAS Catalog `md5sum.txt` and YAML sidecar. Total original compressed bytes read **1,013,815,950**. Overall technical status `PASS_ALL_REQUESTED_ORIGINAL_MD5_AND_ALLELE_COUNTS`; scientific `HOLD`.

| Phenotype | Accession | Full gzip bytes | MD5 (computed = publisher) | Outcome rows | BMI matched unique variants | BMI SNP absent in outcome | Match rate |
|---|---|---:|---|---:|---:|---:|---:|
| AIS | GCST90104540 | 283,482,604 | `1bfe4ae8a5042fb24cb0a562b05b0d2f` | 7,482,032 | **1,006,932** | 967 | **99.9041%** |
| LAA | GCST90104542 | 230,156,887 | `71ea9fd7c06791a831861053b6a1cfbe` | 5,784,788 | **956,929** | 50,970 | **94.9429%** |
| CE | GCST90104541 | 263,057,516 | `c252262586a6ea84b985f30dd7ab32be` | 6,659,793 | **996,662** | 11,237 | **98.8851%** |
| SVD | GCST90104543 | 237,118,943 | `f41bc5f55b3cab238e305e3e86363f69` | 5,974,028 | **967,356** | 40,543 | **95.9775%** |

The 1,007,899-SNP BMI reference panel was independently rebuilt in Linux from GIANT corrected original and strict MRCIEU/Zenodo EUR reference intersection. Candidate canonical plaintext SHA256 **`ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf`** matched the prior Orca and Linux results.

### Explicit allele-orientation accounting

The four original phenotype files use `chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error` among their fields (without rsID or per-SNP N). At the uniquely matching BMI coordinate:
- AIS: `same` **1,006,932**; swapped 0; complement_same 0; complement_swapped 0; mismatch 0; collision 0.
- LAA: `same` **956,929**; swapped 0; complement_same 0; complement_swapped 0; mismatch 0; collision 0.
- CE: `same` **996,662**; swapped 0; complement_same 0; complement_swapped 0; mismatch 0; collision 0.
- SVD: `same` **967,356**; swapped 0; complement_same 0; complement_swapped 0; mismatch 0; collision 0.

**Review flag:** 100% `same` allele orientation in all 4 outcomes may reflect a shared standardized/forward-strand reference system. It is **not independent proof of original allele orientation**. Confirm source processing and a third-party allele reference, do not declare allele-direction validation complete for all purposes. The current workflow reports source-orientation classification and a few sanity examples; it **does not export full per-SNP harmonized effects**.

All four had **0 overlapping BMI candidate coordinate collisions**, and the outcome-original rows partition as `candidate coordinate matches + non-BMI rows = original outcome rows`. No false `0 matches` due to a missing archive was counted: full raw MD5 passed.

## Study case/control numbers, with strict limits

The GIGASTROKE EUR paper and subsequent method tables report **1,234,808 controls** for each listed phenotype and cases **AIS 62,100, LAA 6,399, CE 10,804, SVD 6,811**. These match the GWAS Catalog YAML's `sample_size` fields after addition. Useful sources:
- Mishra et al. (2022), *Nature*, `https://www.nature.com/articles/s41586-022-05165-3`
- [Peer-reviewed stroke GWAS descriptions, 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13279166/)
- [Published phenotype case/control table](https://cdn-links.lww.com/permalink/cm9/c/cm9_2024_11_18_tian_cmj-2024-1505r2_sdc1.pdf)

A common arithmetic case-control harmonic proxy is `N_eff = 4*N_cases*N_controls/(N_cases+N_controls)`:

| Phenotype | Total N | Approx. harmonic N_eff (illustration only) |
|---|---:|---:|
| AIS | 1,296,908 | 236,505.83 |
| LAA | 1,241,207 | 25,464.04 |
| CE | 1,245,612 | 42,841.16 |
| SVD | 1,241,619 | 27,094.55 |

**Do not substitute these values blindly as per-variant LDSC N.** The original files have no per-SNP N, meta-analysis studies can differ in case/control sample composition by SNP, and binary-trait liability-scale assumptions, cohort overlap and ascertainment require a justified analysis plan. The absolute total N is not equivalent to case/control effective N.

## Independent evidence guard

New `scripts/ckm_stroke/review_gigastroke_full_evidence.py` checks all four original pins/MD5 status, expected allele and coordinate row accounting, exact candidate hash and outcome counts, direction-sign partitions, and scientific HOLD. It writes a **small aggregate-only JSON** and fails closed on stale or inconsistent sources. Eight offline regression tests in `tests/test_ckm_stroke_review_gigastroke_full_evidence.py` cover missing outcome, bad checksum, unexpected rows, candidate hash changes, invalid scientific GO, and allele-direction accounting.

This check is deterministic and does not redownload the 1.014 GB originals. The existing full-original GitHub workflow can run it immediately after the source scan whenever the heavy job is deliberately triggered.

## Scientific interpretation: technical PASS ≠ h2 or rg approval

**Technical**: four official compressed originals MD5 checked, parsed all rows, matched and classified BMI candidate loci. **Scientific**: `HOLD`, `ready_for_ldsc=false`, `h2/rg/MR/coloc = NOT_RUN`.

Remaining gates in priority order:
1. Independent allele orientation check against a separate EUR reference and GIGASTROKE normalization documentation, given uniform `same` labels.
2. Decide how to set case-control effective N in LDSC and heterogeneity/error checks. Per-SNP N is absent; no N_eff was entered into a munged outcome file.
3. Quantify GIANT/UK Biobank sample overlap with GIGASTROKE and assess intercept bias/sensitivity.
4. Produce separately approved, allele-aware signed outcome Z input, with full attrition accounting and exact source provenance.
5. Only if per-phenotype h2/intercept and regression quality meet preregistered thresholds, consider LDSC rg; never interpret an h2/rg estimate created solely because the GitHub job was SUCCESS.

**PR #284 remains Draft, unmerged.**
