# CKM × Stroke — single-trait LDSC h²/intercept *pre-regression* gate (2026-11-17)

**Scientific state: `HOLD`, `ready_for_ldsc=false`.** This is a source-confirmed **preflight**, not LDSC h², intercept, confounding ratio, partitioned h², rg, MR or colocalization. Input descriptive `Z²` aggregation **does not fit LD Score regression** or estimate heritability. All raw and derived association files remain in ignored local/CI temporary directories.

## 1. Scope and exact source provenance

**GIANT BMI primary 2018 (Yengo et al.)** independently built original EUR GRCh37 candidate `SNP N Z A1 A2`: 1,007,899 SNPs; canonical full plaintext SHA256 `ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf`. N field is a per-SNP integer written from the verified source; Z is source BETA/SE. Not yet LDSC-munged/approved. The corrected original BMI full source, HM3 pair and LDscore provenance were independently verified by earlier CI.

**GIANT BMI alternate 2015 (Locke et al.)** independently defined source-selected GRCh37 EUR candidate `SNP N_REPORTED Z_FROM_BETA_SE A1 A2`: 1,070,753 SNPs, canonical full plaintext SHA256 `a1cc5ec06cc81d2bb0a583eb76c0f7e81784cd5dc3f622270de61023f95de4bf`. Signed Z is calculated from original **rounded** 2015 BETA and SE. Published source P was genomic-control-adjusted and may not match Z-derived P. The 2015 N is published per-SNP, sometimes **fractional**, and is preserved as given. **Never silently integer-round 2015 N** or rename this intermediate into an LDSC-approved `N` column. See `docs/ckm_stroke_locke_2015_independent_panel_20261117.md`.

**Stroke outcomes**: publisher GWAS Catalog original EUR GRCh37 official `GCST90104540` AIS, `GCST90104542` LAA, `GCST90104541` CE, `GCST90104543` SVD. Each original gzip was fully downloaded and original MD5 verified in prior GitHub CI [#43889295894](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43889295894). Original 11-column header includes beta, odds_ratio, standard_error, effect_allele, other_allele, p_value and allele frequency, **but lacks rsID and per-variant N**. Coordinate/allele harmonization with BMI 2018 had been validated across source rows previously.

**New scale QC** (2026-11-17): official original first 131,072 compressed bytes per stroke subtype, HTTP Range 206 and original GWAS Catalog archive host. **13,201 complete association rows** were checked: AIS 3,413; LAA 3,256; CE 3,273; SVD 3,259. `beta≈log(odds_ratio)` holds across all sampled valid rows with maximum discrepancy **4.9730×10⁻¹⁵**. This verifies **sampled log-odds scale**, not complete effect orientation or per-SNP sample N. Four different subtype original source accessions are checked against a strict primary manifest; never swap CE/LAA.

**Cohort provenance:** Original Locke 2015 Nature supplement and corrected Mishra 2022 Nature supplement were SHA256-verified and independently re-read on Linux [#44000371973](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44000371973). They share **10 explicitly named cohorts and 6 source-family links**, without person-level overlap records. Neither pre-UKB Locke BMI exposure nor Yengo BMI is proven independent of GIGASTROKE. Exact overlap count and sample-overlap covariance are **unknown**.

## 2. New descriptive source Z² results (Orca 2026-11-17; *no LDSC run*)

| Source | SNPs | Mean Z² | Median Z² | Lambda GC (median Z² / 0.4549364) | SNP with Z²>80 | Per-SNP N median | Fractional N rows |
|---|---:|---:|---:|---:|---:|---:|---:|
| Locke 2015 BMI EUR | **1,070,753** | **1.25082047** | 0.48501276 | 1.06611107 | 159 | 233,648.0 | **11,506** |
| Yengo 2018 BMI EUR | **1,007,899** | **3.98119481** | 1.24408284 | 2.73463011 | 2,041 | 688,787.0 | **0** |

Both means exceed the LDSC FAQ's weak-signal floor of about **1.02**, but this is only a **diagnostic observation**, not LDSC scientific GO. Differences arise from sample size, source inclusion, phenotype processing, rounding/GC correction and SNP panels. **A larger mean Z² or lambda does NOT prove larger heritability or population stratification**, which requires regressing chi2 against LD Scores under a properly specified model. The candidate sets are independently selected, not identical SNP panels.

## 3. Per-outcome h² input assessment and published reference only

The primary corrected 2022 Supplementary Table 1 provides EUR counts AIS 62,100 / CE 10,804 / SVD 6,811 / LAA 6,399, with 1,234,808 listed EUR controls. **These are not per-variant effective sample sizes**. The global and sum-of-cohort harmonic N_eff scenarios differ materially, and at least the HVH2 reported cases lack standalone control counts. Two arithmetic scenarios in `docs/ckm_stroke_gigastroke_cohort_overlap_neff_20261117.md` are explicitly **NOT APPROVED** for automatically writing an LDSC outcome N column.

| Stroke subtype | Primary GWAS | Published 2022 *Cov-LDSC* observed-scale h² (SE) | Published h²/SE | Present decision |
|---|---|---|---:|---|
| AIS | GCST90104540 | 0.0049 (0.0005) | 9.80 | `HOLD_N_EFF` |
| LAA | GCST90104542 | 0.0020 (0.0012) | **1.67** | `HOLD_N_EFF` + `LOW_H2_PRECISION_REVIEW` |
| CE | GCST90104541 | 0.0074 (0.0012) | 6.17 | `HOLD_N_EFF` |
| SVD | GCST90104543 | 0.0033 (0.0011) | **3.00** | `HOLD_N_EFF` + `LOW_H2_PRECISION_REVIEW` |

The published values are **historical benchmarks from Mishra et al. Cov-LDSC**; no new h² estimates are produced here. The h²/SE threshold **4.0** is a conservative local prespecified *review flag*, **not a universal LDSC inclusion rule**. Other future diagnostics and sensitivity analyses may change decisions.

## 4. Proposed explicit diagnostic quality gates before LDSC h², then rg

1. **Source identity and association statistics**: original MD5 or independent SHA256, correct EUR ancestry/GRCh37, consistent A1 effect direction, rsID↔chr:BP matching and unique variants; source-derived quantitative Z validated against reported P/GC adjustment. Outcome log(OR) sampled scale is supported but full effect and standard-error QC remains.
2. **Effective sample N**: use original per-SNP N for continuous BMI where defensible, document reported fractional values and any rounding convention. For AIS/LAA/CE/SVD document **meta-analysis cohort weights, sample missingness and case-control effective N assumptions**, with alternative sensitivity inputs recorded separately. Do **not** auto-fill total N or infer per-SNP N from the corrected Table1 totals.
3. **Reference**: independently pinned MRCIEU/Zenodo EUR GRCh37 HapMap3 SNPs, 22 autosomal LD score sets and `.M_5_50` scalars; verify final regression SNP count after matching to GWAS. Per [official LDSC FAQ](https://github.com/bulik/ldsc/wiki/FAQ) **fewer than 200k regression SNPs is a warning**, 600k+ is preferable for precision (especially rg); the current source candidate count is **not yet the final regression SNP count**.
4. **One-trait h² diagnostics when authorized**: log number of regression SNPs, mean χ², lambda GC, observed-scale h² with SE, intercept with SE, h²/SE and intercept-related warnings; review mean χ² ≲1.02 for low polygenic signal, negative h², or unreliable h² denominator. **Do not fix the univariate LDSC intercept to 1** without specific justification, particularly for GC-corrected 2015 BMI. Explicitly differentiate observed vs disease liability scale; prevalence inputs need documented sources.
5. **Cross-trait rg only after single-trait evidence**: review low-h² subtype results (especially LAA/SVD), record study overlap and outcome effective-N assumptions; do **not constrain cross-trait intercept to 0** when overlap is unquantified. A zero sample-overlap intercept is not justified by pre-UKB publication or shared cohort names. No rg/MR is run in this preflight.

**Interpretation:** Official LDSC guidance notes that GC correction can bias h² downward and intercept below 1, that incorrect N columns bias h², and that mean χ² ≤1.02 can indicate insufficient polygenic signal. These are methodological cautions, **not automatic statistical thresholds that override source validity**. See [official LDSC heritability tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) and [FAQ](https://github.com/bulik/ldsc/wiki/FAQ).

## 5. Executable contract and safe CI scope

- `scripts/ckm_stroke/audit_gigastroke_logor_prefix.py`: bounded Range original GWAS prefix parser, four separate EUR GWAS accession pins, sampled log(OR) check, 5 offline fail-closed unit tests.
- `scripts/ckm_stroke/audit_ldsc_preregression_readiness.py`: verifies source selection provenance and entire independent BMI candidate plaintext checksums, duplicate IDs, 2015 fractional N vs 2018 integral N, aggregates mean/median Z² and lambda only, cross-checks original 4 outcome log(OR) and corrected h² reference, rejects unsubstantiated overlap and writes `HOLD`. 10 offline fail-closed unit tests.
- CI `.github/workflows/ckm-stroke-ldsc-singletrait-preregression.yml`: lightweight 15 tests on file changes, source-heavy independent replication on opt-in `[ldsc-preregression-full]` commit message or `workflow_dispatch`. Recreates 2015/2018 panels from original sources into ignored artifacts, independently scans source outcome logOR prefix, verifies original corrected two cohort supplements; exports only small aggregate JSON, then deletes temporary GWAS candidate files. **No LDSC executable**, no auto scientific GO, no remote server.
- Local Orca ignored reports: `artifacts/ckm_stroke/gigastroke_logor_prefix_20261117.json`, `artifacts/ckm_stroke/ldsc_h2_preregression_readiness_20261117.json`.

**Completion criteria for this task:** strict original provenance + source χ² descriptive metrics + effect scale + scientific HOLD gate + 15/15 tests + GitHub Linux independent reproduction. **Not completed / not attempted:** per-SNP outcome effective N, cohort participant overlap count, verified final LDSC munging/scaling, single-trait LDSC regression and interpretation, any rg or MR. Draft PR #284 stays unmerged.
