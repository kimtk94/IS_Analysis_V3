# CKD Top3 LD PSD eigenspectrum audit — numerical method and acceptance

The initial [ARPACK run](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39050802209) continued running without a verified eigenvalue result at audit review time; this is **not** PSD PASS. Near-singular large reference-LD matrices can slow iterative smallest-eigenvalue solvers. The read-only inspector was changed to use `scipy.linalg.eigh(... subset_by_index=(0,0), driver='evr')` (LAPACK direct dense symmetric eigenpair) in [new run](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39051010599).

Source files: preexisting `/srv/is-analysis/data/ckd/stage2c_ld/ld/{SDCCAG8,INHBC,UMOD}.ld.bin` and `.ld.meta.tsv`. All copied to a transient GitHub Actions runner over authorized read-only SSH; research server files and live settings unchanged.

**What the check tests:** exact minimum eigenpair, residual `||Rv - lambda v||`, numerical PSD tolerance `max(1e-5, 10*residual)` and reference LD source SHA256. A diagonal matrix `D` of ±1 allele-orientation signs transforms signed LD by `DRD` and **preserves every eigenvalue exactly** (since `D=D^-1`). The test need not diagonalize the signed matrix again.

**What it cannot establish:** cause of UMOD SuSiE nonconvergence, true QTL/GWAS allele harmonization, ancestry-matched LD, donor-matched LD, overlap or protein-level causality. PSD of a reference correlation matrix is necessary numerical consistency but not sufficient scientific validity. Near-zero negative eigenvalues may be float32 rounding artifacts and must be interpreted against residual/tolerance, not as a hard biological rejection.

**Decision gate:** label `PSD_NUMERIC_PASS` or `PSD_NUMERIC_FAIL` only after GitHub Actions finishes and provides actual minimum eigenvalues. Until then keep `PSD_PENDING`. Follow-up comparison with original signed effect alleles remains separate.
