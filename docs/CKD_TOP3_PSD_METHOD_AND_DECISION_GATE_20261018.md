# CKD Top3 LD PSD eigenspectrum audit — numerical method and acceptance

The initial [ARPACK run](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39050802209) continued running without a verified eigenvalue result at audit review time; this is **not** PSD PASS. Near-singular large reference-LD matrices can slow iterative smallest-eigenvalue solvers. The read-only inspector was changed to use `scipy.linalg.eigh(... subset_by_index=(0,0), driver='evr')` (LAPACK direct dense symmetric eigenpair) in [new run](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39051010599).

Source files: preexisting `/srv/is-analysis/data/ckd/stage2c_ld/ld/{SDCCAG8,INHBC,UMOD}.ld.bin` and `.ld.meta.tsv`. All copied to a transient GitHub Actions runner over authorized read-only SSH; research server files and live settings unchanged.

**What the check tests:** exact minimum eigenpair, residual `||Rv - lambda v||`, numerical PSD tolerance `max(1e-5, 10*residual)` and reference LD source SHA256. A diagonal matrix `D` of ±1 allele-orientation signs transforms signed LD by `DRD` and **preserves every eigenvalue exactly** (since `D=D^-1`). The test need not diagonalize the signed matrix again.

**What it cannot establish:** cause of UMOD SuSiE nonconvergence, true QTL/GWAS allele harmonization, ancestry-matched LD, donor-matched LD, overlap or protein-level causality. PSD of a reference correlation matrix is necessary numerical consistency but not sufficient scientific validity. Near-zero negative eigenvalues may be float32 rounding artifacts and must be interpreted against residual/tolerance, not as a hard biological rejection.

**Decision gate:** label `PSD_NUMERIC_PASS` or `PSD_NUMERIC_FAIL` only after GitHub Actions finishes and provides actual minimum eigenvalues. Until then keep `PSD_PENDING`. Follow-up comparison with original signed effect alleles remains separate.

## Verified execution outcome (2026-10-17 16:49 UTC)

[Actual completed GitHub Actions run 39051010599](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39051010599) **SUCCESS**.

| Locus | n variants | Minimum eigenvalue (reference float32 LD) | Eigenpair residual norm | Tolerance | Status |
|---|---:|---:|---:|---:|---|
| SDCCAG8 | 3,854 | -1.0811689244e-6 | 2.0619700416e-14 | 1e-5 | PSD_NUMERIC_PASS |
| INHBC | 2,679 | -2.1349083682e-6 | 2.8797066793e-14 | 1e-5 | PSD_NUMERIC_PASS |
| UMOD | 3,090 | -1.7592639593e-6 | 2.8940721734e-14 | 1e-5 | PSD_NUMERIC_PASS |

The slightly negative minimum eigenvalues are **within the predeclared numerical tolerance**, not strictly nonnegative in full-precision mathematical terms. These values may reflect float32 rounding of an otherwise positive-semidefinite estimated correlation matrix. Thus the interpretation is **PSD within numerical tolerance**, not proof of perfect exact PSD.

**UMOD causal diagnostic conclusion:** No large negative eigenvalue was detected that would by itself explain the SuSiE convergence issue. The convergence problem remains **unresolved**; reference/population LD compatibility, model specification, sample overlap and trait signal complexity remain relevant explanations and must be examined separately. These tests do not establish any such alternative as the cause.

**Remaining gate:** independently verify effect allele coding/orientation and summary-statistic-vs-reference LD compatibility. Maintain reference LD n=503 limitation and do not promote manuscript causal status.
