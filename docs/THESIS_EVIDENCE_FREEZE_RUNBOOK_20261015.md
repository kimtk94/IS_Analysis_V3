# Thesis evidence freeze — existing-results read-only audit

**Date:** 2026-10-15  
**Status:** IMPLEMENTED_FOR_REVIEW; science results NOT RERUN.  
**Policy:** Follow [Hong-Hee Won / high-impact thesis direction](../docs/RESEARCH_DIRECTION_WON_HIGH_IMPACT_20261015.md) and [scope coverage audit](../docs/RESEARCH_SCOPE_COVERAGE_AUDIT_20261015.md). **No new statistical models.**

## Goal

Transform *existing* IS and CKD aggregate result matrices into an immutable-by-default, SHA256-traceable **source-readiness inventory**. This script does not recalculate effect sizes, GWAS, Mendelian randomization, fine mapping, molecular QTL, or colocalization. Numeric PP.H4 is kept as observed source output but not upgraded to a causal finding or PASS.

- IS: preserve **2,225 unique stable Ensembl IDs**; optionally cross-validate **2,425** gene–region rows against the frozen 2026-10-11 manifest. Missing QTL/LD does not mean biological negativity.
- CKD: require exactly **ACP1, CPVL, F12, GSTA1, GSTA3, HLA-E, INHBC, SDCCAG8, UMOD**; annotate missing EUR/EAS colocalization, unavailable kidney cell context and technical exceptions. A KoGES rsID is **not** proof of participant-level longitudinal validation.
- Each input is hashed and every output is reproducible from the original file. Preserve original source files and controlled data separately. Outputs contain candidate-level aggregates only.

## Original source inputs

Existing-source paths on the research server (verify presence; they are **not** checked in this implementation review):

```text
/srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/
  IS_2225_GENE_EVIDENCE_MATRIX.tsv
  IS_2425_GENE_REGION_EVIDENCE.tsv
  IS_GENE_EVIDENCE_MATRIX_MANIFEST.json

/srv/is-analysis/results/ckd/integrated_evidence/
  CKD_9GENE_EVIDENCE_MATRIX.tsv
```

The IS source is generated in `research/is-broad-discovery-20261009`. The CKD 9-gene matrix producer is `scripts/build_ckd_integrated_evidence.py` on `research/ckd-thesis-20260927`. Both research branches have diverged from main. **Do not merge either blindly or re-run analysis just to feed this audit.**

## Run on an authorized computer after read-only input inventory

The following script is intentionally **not** a server deployment script. It only reads input files and creates new derived files in a **new, empty** output folder. It never runs `git pull`, `set -e`, destructive cleanup, or controlled individual-level genotype access.

```bash
python3 scripts/audit_thesis_evidence_freeze.py \
  --is-input /srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2225_GENE_EVIDENCE_MATRIX.tsv \
  --is-regions /srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2425_GENE_REGION_EVIDENCE.tsv \
  --is-manifest /srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_GENE_EVIDENCE_MATRIX_MANIFEST.json \
  --ckd-input /srv/is-analysis/results/ckd/integrated_evidence/CKD_9GENE_EVIDENCE_MATRIX.tsv \
  --out-dir /srv/is-analysis/results/thesis_evidence_freeze/20261015_v1
```

If one domain's input is unavailable, run *only the available domain* with its own new output folder; do not manufacture placeholder rows. CLI default IS count is **2225**. `--expected-is` exists for synthetic tests, not permission to quietly change the frozen discovery denominator.

Outputs:

```text
IS_EVIDENCE_READINESS_FROZEN.tsv
CKD_EVIDENCE_READINESS_FROZEN.tsv
EVIDENCE_FREEZE_PROVENANCE.json
```

The provenance includes original source paths, SHA256, bytes, row counts and dataset-level coverage. It clearly sets `scientific_validation_completed: false`. Existing nonempty output directories are rejected.

## Validity / interpretation gates

| Check | Audit behavior | What remains independently necessary |
|---|---|---|
| IS gene denominator | Fail on count mismatch, invalid or repeated Ensembl ID | Provenance of GWAS outcome, ancestry and lead alleles |
| IS region joins | Fail if gene missing/unexpected or 2,425 frozen pair mismatch | Genome build, lead vs LD independence, candidate attrition |
| QTL donor LD | Preserve missingness in every gene row | Study-appropriate signed LD and matched alleles |
| Posterior | Fail on non-numeric or impossible PP.H4 outside [0,1] | Priors, SNP coverage, allele alignment, MR consistency |
| CKD 9-gene integrity | Fail on missing, duplicate or extra gene | Full MR summary, tissue resources, prior discovery overlap |
| KoGES | Always `koges_longitudinal_validated=FALSE` here | Actual controlled access, eligible genotypes and repeated eGFR modeling |
| Manuscript claim | Always label source evidence `SOURCE_DATA_PRESENT_UNVERIFIED` | Independent source-based scientific audit before promoting claim |

### Do not mistake technical passes for biological results

A synthetic GitHub Actions PASS confirms the audit script's fail-closed behavior, but **does not** establish whether existing MR/coloc estimates satisfy their biological assumptions. Review source outputs/figures before marking manuscript gates PASS.

## CI

`python -m unittest discover -s tests -p test_thesis_evidence_freeze.py -v`

[Workflow](../.github/workflows/thesis-evidence-freeze-ci.yml) runs synthetic denominator, missingness, nine-gene, posterior, manifest hash and non-overwrite tests. It has no access to server data by design.

## Next acceptance actions

1. Verify that both source files actually exist and record their SHA256 and git HEAD from the authorized server.
2. Run the read-only freeze **only on existing data**; review missingness and QC flags without inventing any causal labels.
3. Produce IS whole-universe supplementary source table plus representative locus figures; produce CKD nine-gene integrative evidence figure and prior-paper comparison.
4. Stop new discovery/model development after high-impact relevant workflow coverage is demonstrated. Record missing BBJ donor/cohort LD and KoGES controlled data as limits, not invisible PASS.
