# HCC Stage29 RNA missing-feature audit — independent technical source check

Audit date: 2026-10-31 (KST). Source Suter `cbg-ethz/HCC` immutable commit `438f52271d0a86a6701baa7a96e80c516ba575c2`.

## Source identity

- `HCCfull.rds` SHA256: `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`.
- `HCCfullz.rds` SHA256: `9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007`.
- RNA `HCCfull$T` has 20,356 features × 50 biopsy columns; `HCCfullz$T` has 19,768 × 50; all 19,768 fullz feature names are contained in full and exactly 588 full-only names remain.

## Direct independent verification

GitHub Actions [#41187760168](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41187760168) source SHA-pinned **SUCCESS**, artifact [#13037775954](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41187760168/artifacts/13037775954).

- Excluded 588/588 RNA features have exactly zero sample-wise standard deviation across the 50 source biopsy measurements.
- Retained 19,768/19,768 RNA features have positive sample-wise SD; retained minimum sample SD **0.008226619**; retained SD median **0.828504693**.
- Hence **the exact membership partition of HCCfullz RNA features equals filtering HCCfull RNA to positive-variance rows**, without exporting the gene list. Because no author-written original fullz construction code has been independently found, the **operational data rule** is proven, not the author's historical motivation or complete generation workflow.

## Link to Stage28

Stage28 CI #41186682016: independent training/held-out feature splits prove the same full→fullz 50-column permutation in RNA, protein and PP and restore all three heldout gene-wise Spearman rho medians to 1.

Stage28B CI #41187009479: after candidate column order, held-out feature-wise `row-zscore(HCCfull)` matches `HCCfullz` at **all tested cells absolute error <10⁻⁶**. Together with this Stage29 zero-variance feature filter verification, continuous layer *quantitative transformation* is fully **reconstructible for all retained features**; this does not prove true biopsy label orientation or original Ng 2022 raw RNA RSEM/normal reference preprocessing.

## GO / HOLD / next

- **GO_TECHNICAL**: source-verified 588 RNA feature omission = all zero-variance features.
- **HOLD_META**: identify whether all `HCCfullz` column names truly match their biological biopsies. Patient-group and clinical/sample labels must be reconciled against the original official Ng Supplementary Data 1, with no silent relabeling.
- **HOLD_ORIGINAL_NG_EFFECTS**: current `HCCfull` patient-aware rho is not the per-gene original Ng 2022 Figure2b rho numerical reproduction.
- **NEXT**: wait for independent multi-job Stage29 [#41187890718](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41187890718) (source patient group orientation, precise original 5,481 per-gene rho vs common 5,253 current estimates and separate feature filter).
- Both PR #251 / #252 remain Draft. V03 stays PARTIAL, V07 stays PARTIAL. Biological or clinical mechanism PASS none.
