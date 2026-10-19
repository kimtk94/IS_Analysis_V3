# CKD UMOD rs11860566 paired exploratory SuSiE sensitivity — 2026-10-18

## Verified results

- [30-iteration completed run 39136744972](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39136744972) SUCCESS
- [100-iteration completed run 39136804380](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39136804380) SUCCESS
- Earlier full 3090 SNP attempt [39126433876](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39126433876) **FAILED** due to SSH broken pipe before either arm reported a result; it does not establish nonconvergence.

Both successful runs used identical fixed selection of top 180 |aligned eGFR Z| SNPs plus rs11860566, selected once before either fit. Included arm 181 SNPs; excluded arm 180. Corresponding signed LD rows and columns and Z-score rows were subset identically. SuSiE `susie_rss`, `L=3`, `estimate_residual_variance=FALSE`, `refine=FALSE`, and the same N convention were used.

| Maximum iterations | rs11860566 included | rs11860566 excluded | Sets returned |
|---|---|---|---|
| 30 | not converged at 30 | not converged at 30 | 3 / 3 |
| 100 | not converged at 100 | not converged at 100 | 3 / 3 |

The 100-iteration exploratory study shows **dropping rs11860566 alone does not restore convergence in this fixed signal-enriched subset**. It does not prove that the SNP is unrelated to full-universe 3090-variant nonconvergence, nor does it identify the cause.

## Scientific guardrails

The subset is a signal-selected diagnostic and differs from the full Stage2C SuSiE analysis (3090 SNPs; original 1000 iterations; possibly different method configuration). The returned three sets must not be presented as validated credible sets from a converged fit. This analysis used the existing external 1KG EUR LD matrix, not cohort-matched LD.

No source datasets, frozen Stage2B/Stage2C outputs, PR base or LIVE settings were changed. The code was committed only to the research branch.

Next: isolate numeric conditioning / signed LD and evaluate convergence with stronger controls (e.g., reference LD regularization within a separate diagnostic fit, alternate smaller unsupervised contiguous locus windows) before any scientific inference.
