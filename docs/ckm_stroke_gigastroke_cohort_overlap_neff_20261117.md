# CKM × Stroke — GIGASTROKE corrected primary cohort overlap and LDSC N gate

**Evidence:** 2026-11-17 KST · original publisher-corrected Mishra et al. (2022) Supplementary Tables 1–56. **Source SHA-256:** `7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4`, **756,910 bytes**. Original archive not committed; processed locally in Orca and released only as aggregate JSON.

**Primary source:** [Publisher-corrected GIGASTROKE 2022 Supplementary Tables](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/41586_2022_5492_MOESM2_ESM.xlsx) · [Nature original](https://www.nature.com/articles/s41586-022-05165-3) · [2022 publisher correction](https://doi.org/10.1038/s41586-022-05492-5) (correction specifically addressed Table 55, but this is the corrected 56-table bundle).

## 1. Original Supplementary Table 1 — verified UK Biobank overlap hazard

The `Suppl.Table1` original worksheet's **row 15 (UKBiobank)**, in the **Incident & prevalent combined** panel (columns J–P), contains:
- **AIS case count: 2,211**, **control count: 369,419** (also AS cases 2,932).
- CES, SVS and LAS columns **NA**, rather than zero observed events or a demonstration of zero participant overlap.
- The UK Biobank row is classified **Population-based, incident**.
- This is one listed discovery EUR cohort among the **44 named EUR rows** in Table 1 (rows 5–48). The publisher's **TOTAL EUR** row 49 reports AIS 62,100, CES 10,804, SVS 6,811, LAS 6,399, and controls 1,234,808.
- UK Biobank's AIS cases are **3.5604% of the EUR AIS case total**, and its listed controls are **29.9171% of the EUR control total**. These are cohort contribution percentages, **not percentages of genetically identical UKB subjects overlapping with the BMI exposure**.

The selected GIANT/Yengo (2018) BMI GWAS includes **456,426 UK Biobank EUR participants**, as reported in its [published Methods](https://pubmed.ncbi.nlm.nih.gov/30124842/) (total BMI meta-analysis average N 681,275). The two GWAS thus explicitly share a recruitment source. Individual identifiers and exact common-person membership are unavailable in these public summary data.

**Overlap status: `POSSIBLE_CONFIRMED_SHARED_COHORT_BUT_INDIVIDUAL_COUNT_UNKNOWN`.** Never set a numeric overlap fraction or claim independent samples on the basis of cohort names alone. Direct UKB exposure overlap is an evident risk for AIS; for LAA/CE/SVD, the UK Biobank row's NA subtype counts do not exclude overlap through other cohorts. Other shared cohorts (e.g. HRS, MEGASTROKE/GIANT participants) must be checked separately if evidence permits.

## 2. Official accession mapping — primary publisher precedence

The original workbook uses the column order **AIS, CES, SVS, LAS**, not the shorthand order AIS, LAA, CE, SVD. The official EBI GWAS Catalog YAML was independently checked on **2026-11-17**, including accession, `trait_description`, GRCh37 and total sample size. No LAA/CE swap is permitted.

| Phenotype | Official EBI accession | EUR cases | EUR controls | TOTAL EUR N |
|---|---|---:|---:|---:|
| AIS | `GCST90104540` | 62,100 | 1,234,808 | 1,296,908 |
| LAA/LAS | **`GCST90104542`** | **6,399** | 1,234,808 | 1,241,207 |
| CE/CES | **`GCST90104541`** | **10,804** | 1,234,808 | 1,245,612 |
| SVD/SVS | `GCST90104543` | 6,811 | 1,234,808 | 1,241,619 |

Some 2026 secondary supplementary tables appear to reverse 41/42; they do not overrule the **actual publisher metadata** and corrected 2022 Table 1. The existing `config/ckm_stroke/gigastroke_eur_official.tsv` was correct and remains unchanged.

## 3. Two alternative case–control N_eff scenarios — *neither approved*

Original GWAS files do **not** contain per-variant N. To expose the uncertainty, two arithmetical **illustrations** (not LDSC input recommendations) were computed from the **corrected 2022 primary cohort table**:

- **Total-first scenario:** `4 * TOTAL_EUR_CASES * TOTAL_EUR_CONTROLS / (CASES + CONTROLS)`.
- **Cohort-by-cohort sum scenario:** for each EUR study row having both numeric case and control counts, `sum_i(4*cases_i*controls_i/(cases_i+controls_i))`. **Excludes**, rather than imputes, study rows missing controls, and ignores per-SNP coverage, cohort-specific imputation and meta-analysis weighting.

| Trait | Total-first N_eff scenario | Cohortwise sum scenario | Cases excluded due to missing cohort controls |
|---|---:|---:|---:|
| AIS | 236,505.83 | **195,570.72** | 103 (HVH2) |
| LAA | 25,464.04 | **23,149.57** | 0 |
| CE | 42,841.16 | **39,605.59** | 28 (HVH2) |
| SVD | 27,094.55 | **23,758.72** | 39 (HVH2) |

The differing values **demonstrate meaningful heterogeneity in ascertainment and study proportions**; neither is proven correct as an LDSC per-variant effective sample size. Missing HVH2 controls are **missing**, not zero. These numbers are specifically labelled `illustrative` in the evidence JSON and must **never** populate a final `N` column automatically. Next scientific task: review consortium meta-analysis method, ascertainment and per-variant available studies; seek cohort-level effective N convention or stratified data, and preregister N sensitivity scenarios.

## 4. Original Supplementary Table 3B and 3C — reference benchmarks, not newly run LDSC

| Trait | Published EUR LDSC intercept (Table 3B) | Published observed-scale Cov-LDSC SNP-h² (SE) (Table 3C) | Published h²/SE |
|---|---:|---:|---:|
| AIS | 1.10 | 0.0049 (0.0005) | 9.80 |
| LAA | 1.09 | 0.0020 (0.0012) | **1.67** |
| CE | 1.06 | 0.0074 (0.0012) | 6.17 |
| SVD | 1.08 | 0.0033 (0.0011) | 3.00 |

LAA's weak published h² precision makes an rg denominator potentially unstable. SVD is also lower precision than AIS/CE. The study's **Cov-LDSC**, regression variants and source preprocessing may differ from any future standard LDSC implementation. These are **published historical benchmarks only**: do not present as new results or require exact numerical equality. The univariate intercept cannot by itself determine the number of shared participants between BMI and stroke GWAS.

[Official LDSC FAQ](https://github.com/bulik/ldsc/wiki/FAQ) warns that low heritability and limited SNP coverage destabilize rg and that constraining the genetic covariance intercept depends on whether overlap can be quantified or independence is established. **Do not automatically constrain the intercept.**

## 5. Reproducible evidence contract

- Code: `scripts/ckm_stroke/audit_gigastroke_supplement_cohort_neff.py` (stdlib read-only OOXML, signed original XLSX SHA256, direct Table1+Table3 source cells, official manifest accession crosscheck, two N_eff scenarios, fail-closed HOLD).
- Tests: `tests/test_ckm_stroke_gigastroke_supplement_cohort_neff.py` (8 tests, including CE/LAA accidental accession swap and source corruption).
- Workflow: `.github/workflows/ckm-stroke-gigastroke-cohort-neff-audit.yml`. Runs offline fixtures on relevant commits; complete 0.757 MB publisher source only on explicit `[gigastroke-cohort-full]` commit marker or manual dispatch. Uploads **JSON only**, not original XLSX or GWAS.
- Local Orca evidence: `artifacts/ckm_stroke/gigastroke_corrected_cohort_ldsc_gate_20261117.json`. Source XLSX is stored **only** under ignored local `artifacts/ckm_stroke/references/`.

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_gigastroke_supplement_cohort_neff.py -v
py -3 -m scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff --out artifacts/ckm_stroke/gigastroke_corrected_cohort_ldsc_gate_20261117.json
```

**Technical evidence gate = PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY. Scientific gate = HOLD**, `ready_for_ldsc=false`, no LDSC h²/rg/MR/coloc run. PR #284 stays Draft, unmerged.

## 6. Next study decision

1. Prefer a **UK Biobank-excluded BMI GWAS sensitivity exposure** (e.g. a distinct pre-UKB GIANT BMI summary statistic) where comparable EUR meta-analysis files are independently verified; reduced power and phenotype differences must be explicitly discussed. Do not replace the main BMI exposure without a preregistered design.
2. Resolve case-control meta-analysis N / effective N using original methodology, per-cohort coverage and, where possible, available study-level N; carry alternative Neff choices as **assumption sensitivity**, not an arbitrary fixed constant.
3. Develop an individual-variant QC report for 4 GIGASTROKE outcome sign effects while keeping raw data access and metadata provenance reproducible.
4. Run **h² / mean χ² / intercept QC only after** legal/reference and effective-N assumptions are documented; prioritize AIS/CE for initial power checks, flag LAA and potentially SVD for limited genetic heritability precision. Prohibit rg inference when diagnostics indicate poor h² estimability.
