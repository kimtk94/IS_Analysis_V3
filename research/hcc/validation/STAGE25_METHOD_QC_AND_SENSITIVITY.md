# HCC Stage25 — Preprocessing provenance and patient-level sensitivity gates

**Status 2026-10-31:** Workflows created and queued; do not assign statistical outcomes until job artifacts and logs are inspected. Source review is separate from numerical resampling evidence.

## Source documents (ranked by epistemic priority)

1. **Original Ng 2022 Nature Communications (published final):** https://www.nature.com/articles/s41467-022-29960-8 . Primary quantitative methods: GRCh37 STAR / RSEM; gene-level upper-quartile-normalized expected counts to 1000; log2 tumor-versus-median-normal liver comparisons; SWATH-MS proteome intensity normalized against pooled healthy liver sample and KNN missing-value treatment as described; raw expected counts with TMM/edgeR quasi-likelihood F-test for RNA differential expression. mRNA–protein across-tumor Spearman benchmark n=5481, median 0.287.
2. **Original Suter 2022 PLOS Computational Biology (published final):** https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009767 . Original normalized data from Ng; extra imputation and batch correction **where applicable**; hypermutated exclusion; 50 biopsies from 48 patients. Final preprocessing appendix: https://doi.org/10.1371/journal.pcbi.1009767.s013 (S5 Appendix).
3. **Suter 2022 ETH doctoral thesis (related earlier manuscript version, not final S5):** https://www.research-collection.ethz.ch/server/api/core/bitstreams/00e2e04a-78e0-414b-821e-36cf12f3e3c1/content . Earlier S7 Appendix on 2022 manuscript pp. 31–32 of chapter indicates: RNA upper-quartile-normalized log2 counts; protein log2 and median-subtraction, filter below half-sample detection; **imputation only for clustering** while differential expression used unimputed values; phosphoproteome additionally batch-adjusted. **Important**: this is informative but NOT a substitute for confirming which methods apply to which layers and HCCfull.rds in the final 2022 S5 appendix; treat as `PREPRINT_METHOD_CONTEXT_ONLY`.

## Explicit scientific guardrails

- Avoid applying the generic phrase "all omics were KNN-imputed and batch corrected" to any individual source layer. The earlier S7 states different rules for T, P and PP. The 2022 published source simply says imputation/batch correction *where applicable*.
- The Stage24 script names its source `Suter_preprocessed_normalized_imputed_batch_adjusted_not_raw_counts`; this is a **dataset-level provenance warning**, not proof that imputation and correction were applied to every RNA/protein matrix element. Specific T/P transformation validation remains V03 PARTIAL.
- For Ng original differential expression, raw RSEM expected counts and the correct normal/reference groups are required. Processed HCCfull matrices (50 biopsies) cannot support a count-based raw edgeR reproduction.
- 5,343 current strict reviewed UniProt 1:1 pairs are NOT verified as the original 5,481 correlated genes or 5,490 comparable features. A median match to 0.287 is only a technical sanity check. No new HCC biological claim from matching medians.
- Leave-one-independent-patient-out (48 repeats) and cluster-level patient bootstrap (80 repeats; fixed seed 20261031) are **descriptive sensitivity diagnostics**, not gene-level FDR, a causal mechanism or an inferential guarantee. The bootstrap interval should be described as an empirical patient-resampling percentile interval; gene-set uncertainty, source imputation and batch effects remain.

## Implementation and gates

| Gate | Evidence | GO condition | Status at creation |
|---|---|---|---|
| S5 provenance | SHA256-verified bounded S5 PDF and extracted keyword counts; do not export appendix text | Actual final appendix retrievable; individual T/P/PP step and exact matrix applicability manually reviewed | HOLD until corrected S5 CI succeeds and methods reviewed |
| Sample influence | Patient-aware 48x LOPO; max absolute median-rho change; 80x patient bootstrap percentile interval | Source hash, strict pairs, 48 patients; no individual IDs exported; confirm finite and guardrail PASS | HOLD until analysis CI completes |
| Ng original 5,481-gene set | Original publication historical exact IDs and annotation lineage | Gene membership identity, not merely 5,343/5,481 count ratio | HOLD |
| Clinical validity | Cohort, batch, tumor-purity, etiology and outcome-codebook verification | Prespecified inferential design with independent replication | HOLD |

Expected workflows:
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179644025 — final S5 source verification.
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179643931 — patient LOPO/bootstrap.
- Historical Stage24 actual data: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173779578 (median rho 0.284086, 5,343 strict pairs, 48 patients).
- Implementation PR #251 (Draft); validation log PR #252 (Draft).

**No automatic gate escalation:** V03 PARTIAL and V07 PARTIAL remain until outcomes are reviewed. Source audit CI success without PDF validation still means HOLD, and source-paper scientific methods cannot be marked PASS by keyword counts alone.
