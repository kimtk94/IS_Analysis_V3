# HCC Stage28 — source sample-axis crossmodal permutation and disjoint-feature holdout validation

Date: 2026-10-31. Protocol finalized before interpreting real Stage28 CI output. All source hashes pinned, two PRs remain DRAFT, no patient/feature IDs or inferred index permutation exported.

## Motivating evidence

Stage27b source #41185636861 PASS (four jobs). Source `cbg-ethz/HCC` commit `438f52271d0a86a6701baa7a96e80c516ba575c2`; `HCCfull.rds` SHA256 `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`; `HCCfullz.rds` SHA256 `9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007`.

Same gene/protein feature identity, named same 50 biopsy columns but across-source feature-rank Spearman medians: T RNA -0.002431585, P protein 0.024057623, PP phosphoprotein 0.037118848 (control M and CN=1). Cross-sample best-match same-column winners: T/P/PP each 2/50 and each fitted row-wise maximum is bijective (a possible 50-sample permutation); M 50/50 and CN 46/50. Matrix similarity was `crossprod(within-feature z-scored matrix A,B)/feature_count`: an unbounded aggregation score, not Pearson or a probability.

## Planned independent tests and interpretation

1. **Source immutable + exact two-way feature/sample identity checks.** Align source features and 50 sample IDs by exact string, fail closed on duplicate, missingness or unexpected shape. Keep source full/fullz unmodified.
2. **Independent feature split:** sort source feature IDs and split deterministically by odd vs even feature positions. Infer full→fullz sample permutation only from training half, using the same row-standardized aggregate sample–sample similarity; reverse infer from held-out half to test stability of mapping. Never select features or tune mapping based on observed correlation restoration.
3. **Strict uniqueness:** report whether each train and reverse fit is a 50-element permutation, unchanged sample positions and agreement of the independently inferred maps. A perfect one-to-one match in one run is necessary but NOT proof of correct patient identity.
4. **Evaluate *only held-out features***: aggregate per-feature Spearman rho across 50 biopsies before and after candidate reorder, including median, 10/90% quantiles, median within-feature delta and proportion rho>0.95 or negative. Reverse held-out evaluation checks symmetry. M and CN source layers serve as identity/negative controls.
5. **Cross-omic concordance:** compare the exact 50-element permutation objects in runner-private memory and export **only their Hamming match counts** for T–P, T–PP, P–PP, and M–CN. If continuous modalities agree 50/50 and held-out rho≈1, source column order/label alignment is a strong explanatory mechanism; otherwise test independent normalization, transformation, or source label changes. Always verify original Ng sample clinical/biopsy IDs before correcting or exporting sample names.
6. **Synthetic safety regression:** create a 600×10 simulated matrix and reorder sample columns by known nontrivial bijection. Train on 300 features, evaluate on independent 300 features and require exact inferred mapping and held-out Spearman 1. No source clinical IDs in output or Github logs.
7. **Fail-closed postflight:** sources and dimensions exact, bounds on ranks/proportions, no per-gene/per-patient values, even a positive restoration cannot automatically PASS V03 or V07 because raw Ng 2022 tumor/normal preprocessing and external cohort replication remain unverified.

## Operational GO/HOLD/STOP

- **GO_TECHNICAL:** independent train/test cross-RDS matches agree, all 50 continuous-omic sample maps agree and holdout rank concordance meaningfully restored. **Not** GO_BIOLOGICAL.
- **REVIEW:** per-modality bijections differ, train/test unstable, holdout rank concordance remains low, or original source patient label metadata contradictory. Do not overwrite source RDS or declare one true permutation.
- **HOLD:** original Ng 2022 patient sample manifest and original RDS creation/normal reference mapping not independently reconciled.
- **NEXT:** Preserve Stage24 primary HCCfull RNA–protein comparison, then original Fig2b 5,481 gene-wise rho-vector benchmark with 5,253 exact shared genes and 48 independent patient source controls; independently obtain Gao 2019 actual RNA+protein same-patient matrices before mechanism/clinical claims.

Source validated corrected original `Fig2b` 5,481 gene list and 5,253 intersection in Actions #41185931131 (three jobs PASS; individual gene rho still untested).

Source code / CI Stage28: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186682016 — pending at initial protocol timestamp. Any numerical findings must be appended to D-041 AFTER job logs inspected.

## Actual CI findings / Stage28B (2026-10-31)

**Stage28 PASS_TECHNICAL:** CI https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186682016, artifact #13037146465. On disjoint heldout source feature sets, same 50-sample map learned on 1st 50% and reversed 2nd 50% feature splits (50/50 positions agree). RNA–Protein, RNA–PP and Protein–PP individually inferred maps all agree **50/50**. Original same-column sample matches only **2/50**. Heldout Spearman medians T=−0.001872749→1, P=0.024345738→1, PP=0.036110444→1; fraction of independent test features with rho>0.95 after mapping **100%** for each.

**Stage28B PASS_NUMERIC_TECHNICAL:** CI https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41187009479, artifact #13037476122. On disjoint heldout rows after independent sample map training, actual `HCCfullz[,inferred_column_order]` reconstructed by source `(HCCfull - rowMean)/rowSD` for all common continuous features. Absolute residual median before vs after: T 0.780477824→0 at 9-decimal precision, P 0.824260984→0, PP 0.795386541→0; all tested cells have error <1e-6, and all heldout features have row-max error <1e-5. The exact floating-point error is not asserted zero beyond output precision.

**Interpretation update:** The previously nearly zero cross-RDS rank concordance is not evidence of distinct biological expression patterns; one common source sample-column permutation and subsequent feature-wise standardization explain both rank and numerical discrepancies. However source author code for the permanent fullz construction, true Ng sample-label orientation, why 588 RNA features were omitted, and original Ng 2022 reference-normal transformations remain unverified. **V03 PARTIAL**. No original source files or sample identifiers were modified, uploaded, or retrospectively repaired.

**Immediate following work:** Stage29 original Ng Figure2b correlation effect comparison and 50 biopsy / 48 patient group sensitivity; check source original sample-label map and 588 RNA gene filtering. Original published Fig2b source gene identity (5481) vs current strict (5343) share 5253 Ensembl IDs in corrected v3 CI https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185931131 (3/3 successful).
