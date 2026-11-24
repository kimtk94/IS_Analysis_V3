# CKD TOP3 pQTL–1KG EUR MAF compatibility screen — 2026-10-18

**Actual read-only execution:** [GitHub Actions 39111484302](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39111484302), SUCCESS.

The audit joined the original Stage2C prefilter pQTL summary-statistic `maf_pqtl` to the LD metadata `reference_alt_freq`, computing the difference of **minor allele frequencies**, `abs(maf_pqtl - min(reference_alt_freq, 1-reference_alt_freq))`. This is an allele-frequency SCREEN, not a comparison of effect-allele dosage direction or an empirical study LD correlation matrix.

| Metric | SDCCAG8 | INHBC | UMOD |
|---|---:|---:|---:|
| LD-matched SNPs with MAF | 3,854 | 2,679 | 3,090 |
| Unmatched / missing frequency | 0 / 0 | 0 / 0 | 0 / 0 |
| Median abs MAF difference | 0.009857 | 0.008202 | 0.007475 |
| P95 abs MAF difference | 0.034505 | 0.042884 | 0.030213 |
| Maximum abs MAF difference | 0.113072 | 0.105364 | 0.133779 |
| abs MAF difference > 0.1 | 2 SNPs | 1 SNP | 1 SNP |
| Source/reference frequency availability | PASS | PASS | PASS |

**Scientific meaning:** Median differences below ~0.01 and comparatively small upper tails indicate broadly similar *minor allele frequency* spectra for this sampled matched SNP set. This does not establish LD equivalence or independent validation of original source allele harmonization; selected SNPs are already filtered to match the LD reference and may mask omitted allele/problematic variants.

**Remaining required verification:**
1. Inspect exact SNP identities of the 2+1+1 frequency outliers and assess strand/variant mapping/effect allele harmonization without disclosing individual-level data.
2. Validate raw source pQTL and kidney GWAS effect alleles against original consortium files, not solely pre-harmonized intermediate tables.
3. Evaluate LD-to-summary-statistics compatibility (signed LD against association Z-scores, sample-size differences, strand alignment, outcome-specific locus coverage and study overlap). Do not conclude that UMOD nonconvergence or INHBC discordance has been explained.
4. Leave KoGES longitudinal individual-level validation as NOT VERIFIED.

**Conclusions:** MAF availability and aggregate comparison computed successfully; `LD_COHORT_COMPATIBILITY_VERIFIED = FALSE`, `RAW_EFFECT_ALLELES_VERIFIED = FALSE`. Figures 1–4 are reproducible descriptive aggregates, not evidence of fully established causal mediation.
