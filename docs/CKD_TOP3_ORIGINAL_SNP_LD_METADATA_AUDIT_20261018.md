# CKD Top3 Stage2B-to-Stage2C SNP / allele / LD metadata audit — 2026-10-18

**Status:** METADATA_PASS_FILTERED_REFERENCE_LD_ONLY for **SDCCAG8, INHBC, UMOD**. Scientific causal inference and full numerical LD validation **NOT PASSED**.

Authenticated read-only server execution: [GitHub Actions 39043430679](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39043430679). Source: live server's existing Stage2B harmonized coloc inputs, Stage2C SuSiE prefilter inputs, and 1KG EUR signed reference LD metadata. Code: [audit_ckd_top3_snp_ld_readonly.py](../scripts/audit_ckd_top3_snp_ld_readonly.py).

## Results

| QC | SDCCAG8 | INHBC | UMOD |
|---|---:|---:|---:|
| Stage2B regional ABF rows | 6,269 | 4,777 | 6,048 |
| Stage2C prefilter SuSiE-input rows | 4,409 | 3,171 | 3,376 |
| Final LD metadata/LD matrix dimension | 3,854 | 2,679 | 3,090 |
| Stage2C prefilter rows without LD variant | 555 | 492 | 286 |
| LD variants missing from Stage2C prefilter input | 0 | 0 | 0 |
| Final LD variants also in Stage2B ABF | 3,854 | 2,679 | 3,090 |
| Allele or genomic position mismatches in LD-matched variants | 0 | 0 | 0 |
| LD binary wrong dimensions | 0 | 0 | 0 |
| LD diagonal outside tolerance / nonfinite | 0 | 0 | 0 |

## Interpretation of first-run failure and successful corrected audit

The initial auditor erroneously expected exact equality between the **unfiltered Stage2C SuSiE input** and **post-reference-match LD matrix**. This caused [run 39043336013](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39043336013) to fail. The read-only follow-up established the discrepancy counts above, and the original Stage2C R runner explicitly uses `match(meta$snp, dat$snp)` to order and subset summary statistics to the LD metadata SNP list. Thus, retaining extra Stage2C input variants before LD matching is expected and is not itself a scientific error.

Audit code was adjusted to require **LD SNP set ⊆ Stage2C input SNP set**, unique SNP IDs, exact LD.vars ordering with metadata, exact 4*n*n binary size, signed-orientation flags in ±1, matching effect/other alleles and positions versus Stage2B, and finite near-one LD matrix diagonal. [Success 39043430679](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39043430679).

## Explicit limitations

1. This audit **did not** validate the whole LD matrix (off-diagonal symmetry, eigenvalues, PSD correction), only dimensions and diagonal. Nor did it validate full signed alignment between genotype dosages and all pQTL/eGFR beta signs.
2. It **did not** recompute ABF/SuSiE; **ABF and SuSiE use different region SNP sets** (6,269 vs 3,854 for SDCCAG8, etc.).
3. LD comes from **1KG EUR n=503 reference samples**, not molecular cohort matched LD. Reference-LD sensitivity remains the appropriate claim.
4. UMOD eGFR SuSiE previously failed to converge and stays **NA**. INHBC ABF and SuSiE discordance remains; no technical metadata pass overrides biological inferential mismatch.
5. Sample overlap, pQTL measurement effects, ABF prior sensitivity and controlled Korean longitudinal validation remain unresolved.

**Paper claim gate:** Original SNP/allele *metadata subset* integrity PASS. Full signed-LD numerical integrity, model sensitivity and causal validation **PENDING**. Figures must use the latter distinction, not label this a final research replication.
