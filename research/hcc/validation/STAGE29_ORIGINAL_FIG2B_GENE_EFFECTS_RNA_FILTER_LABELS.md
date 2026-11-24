# HCC Stage29 — original Figure2b gene-wise effect benchmark, RNA 588-feature filter, original patient-label direction

Date: 2026-10-31. GitHub code PR #251 and governance PR #252 remain DRAFT. Original Figures and the public RDS sources are SHA256 pinned. CI is **technical-source evidence**, not scientific/clinical significance.

## Why Stage29 follows Stages 28/28B

Stage28 (#41186682016) independent source feature train/holdout found the **same bijective 50-sample column permutation** in T, P and PP, all 50/50 positions, with the original same-name coordinates matching only 2/50. All heldout Spearman medians restored to 1. Stage28B (#41187009479) showed all heldout T/P/PP numeric z-score residuals <1e-6 after the same sample reorder, with `HCCfullz` reproducing feature-wise `scale(HCCfull)` on shared features. **No original clinical/biopsy label orientation proven.**

Corrected Ng 2022 original Figure2b XLSX (#41185931131) establishes official 5,481 Ensembl IDs and Spearman rho median 0.287059, 5,253 Ensembl string IDs common to current one-to-one 5,343. The original published source file SHA256 `16376438aa16a09a13db417a4ae1f4cc16b87ec2c47630aeae651004501be420`; `Fig2b` physical `sheet8.xml` resolved through `workbook.xml.rels` (numeric sheet mapping).

## Three independent GitHub Actions tasks

1. **P0: Original Fig2b 5,481 gene-specific published rho vector.** Source XLSX real `Fig2b` row A Ensembl / B numeric Spearman, strict SHA256 + 5,481 IDs + median 0.287059 + prior original gene-set digest; store only in runner-private temp. Strict current UniProt reviewed 1:1 5,343 Ensembl–Protein pairs + public pinned `HCCfull` 50-biopsy matrix / official 48-patient index; intersect 5,253 exact ENSG IDs. Compute per-gene RNA–protein Spearman on 50 biopsies, 48-patient averaged and 46 singleton patients, then aggregate compare published *same gene* rho vector with each current vector: effect-vector Spearman/Pearson concordance, directional sign agreement, median absolute effect delta, quantiles. No per-gene vectors in artifacts or logs. **Even perfect agreement is technical corroboration, not reproduction of 2022 original RSEM/normal-reference assay and covariates.**
2. **P0: RNA 588-feature exclusion.** Compare 20,356 `HCCfull` T row ENSGs to 19,768 `HCCfullz`, quantify 588 missing, variance and abundance distributions, actual perfect SD/mean threshold separation, how many omitted are among 588 lowest SD/mean rows. Do **not** infer original author filtering method without independent source-code evidence or direct reproducible selection threshold.
3. **P0: Original biopsy label source orientation.** Pin official Ng Supplementary Data 1 SHA `edf1d2afd3baaedd56edee8536c16f10e4a62cb6cdac929863a062671ee73d25`, verify 50 public RDS sample labels join 48 original patients with both assay flags; group names only as ephemeral 1–48 numeric IDs. Separately infer common T/P/PP cross-RDS column permutation and report number of same-label positions / same vs different original 48 patient groups *conditional on published named biopsy axes*. No automatic label correction; official independent RDS creation and clinical/assay batch lineage needed to decide which orientation is actually correct.

Source `HCCfull.rds` SHA `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`; `HCCfullz.rds` SHA `9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007`; original repo source commit `438f52271d0a86a6701baa7a96e80c516ba575c2`.

## Prespecified interpretation and GO/HOLD

- **GO_SOURCE:** gene-level published effect vector source provenance validated; 5,253 common ENSGs, source cohort 50/48/46, and cross-cohort effect-vector compare returns finite values and verified exact source hashes.
- **GO_TECHNICAL_FILTER:** excluded 588 source feature distribution reproduced; `PASS_FILTER_METHOD` remains HOLD unless actual author filter rule independently found.
- **HOLD_LABEL_ORIENTATION:** both RDS share apparent string biopsy labels but source values differ in original indexing; seeing a cross-sample permutation with correct patient metadata only indicates their labeling cannot *both* be correct under the inferred source transformation. The correct canonical orientation requires an independently known numeric/phenotype source. Never overwrite label axes before proof.
- **HOLD_V07:** comparing original and current per-gene rho vectors is *effect-level source agreement* not independent replication (same Ng/Suter cohort). Published Ng biopsy counts/covariate adjustment and tumor-normal normalization may differ, not covered by matching aggregate median.
- **NO_GO_BIOLOGY:** novel pathway/etiology/prognostic effects, FDR, imaging association and clinical model are not established. Independent Gao 2019 public RNA+Proteome data and verified 1:1 patient/biopsy join remain a hard prerequisite.

Workflow: [Stage29 parallel GitHub Actions](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41187890718). New result numbers and decisions must be entered in `DECISION_LOG.md` **only after** completed job logs and artifact checks. Current V03 and V07 **PARTIAL**.

Expected output: exactly 3 aggregate artifacts, no source sample IDs, no unpublished per-gene rho, no clinical identifiers; runner-private files removed even on failure.
