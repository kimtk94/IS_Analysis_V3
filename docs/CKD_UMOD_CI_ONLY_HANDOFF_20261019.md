# CKD UMOD CI-only handoff and data gate — 2026-10-19

## Verified serverless execution

[GitHub Actions run 39220440127](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39220440127) executed entirely on GitHub-hosted Ubuntu without contacting the failed research server. The job correctly FAILED in `Verify actual UMOD dataset and declare blocker` with `BLOCKED_MISSING_PRIMARY_DATA`. This is deliberate fail-closed behavior and is not scientific nonconvergence. Previous workflow [39216291866](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39216291866) misleadingly returned success with SuSiE skipped.

## Corrected workflow

`.github/workflows/ckd-umod-serverless-susie.yml` now retains validation and execution in **one job on the same runner** to avoid file loss between jobs. Absent files fail the check before R or SuSiE is installed. If actual files exist, the workflow installs R/susieR, runs exploratory LD shrinkage, and publishes a log as a GitHub Actions artifact. Passing CI is not equivalent to SuSiE convergence.

## Required validated research bundle

Place in runner workspace `input/ckd/umod/`:

- `UMOD.tsv.gz` — Stage2C aligned harmonized eGFR and pQTL summary data
- `UMOD.ld.meta.tsv` — original SNP metadata / signs
- `UMOD.ld.vars` — ordered LD variant identifiers
- `UMOD.ld.bin` — signed-interpretation float32 LD matrix in original ordering

These files must be sourced from a trusted backup or reproducibly reconstructed reference-LD workflow. Do not fabricate, substitute demonstration data, or commit large/private research datasets into Git.

## Data blocker

The files are not committed to this Git branch, and Drive filename searches have not located a directly usable archive. Until these exist in the runner, scientific reproduction is **BLOCKED_MISSING_PRIMARY_DATA**. No server-side files or past CKD results were changed.

## Remaining engineering task

Secure a retrievable original data bundle from an authorized backup or re-derive from public 1KG reference and original inputs, verify its checksums and metadata, then stage the bundle in the *same* CI runner job before the input gate. For scientific comparability verify assembly, effect-allele orientation, SNP ID/order, LD provenance, selection filters, and code/parameter versions.
