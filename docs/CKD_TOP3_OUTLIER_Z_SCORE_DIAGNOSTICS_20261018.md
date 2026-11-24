# CKD Top3 high-MAF outliers — harmonized Z-score screen (2026-10-18)

Verified read-only execution: [39125699255](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39125699255) SUCCESS. `scripts/audit_ckd_top3_outlier_z_diagnostics.py`; isolated workflow `.github/workflows/ckd-outlier-z-diagnostics.yml`.

| Gene | rsID | pQTL beta/SE Z | harmonized eGFR beta/SE Z |
|---|---|---:|---:|
| SDCCAG8 | rs35646972 | -0.49801 | 0.00000 |
| SDCCAG8 | rs6663222 | -2.06701 | +1.79921 |
| INHBC | rs11834807 | +12.41242 | -3.29835 |
| UMOD | rs11860566 | +1.96993 | -0.40000 |

Stage2C valid pairs / total: SDCCAG8 4409/4409, INHBC 3171/3171, UMOD 3376/3376; 0 missing SE or unusable beta/SE records.

**Interpretation boundary:** Z-score availability, sign and local relative rank do not independently establish that 1KG EUR LD agrees with GWAS Z covariance. Nor do four local outlier signals establish which SNP caused SuSiE's numerical non-convergence. The original GWAS effect values are rounded in published source, so avoid over-interpreting near-zero values.

Original pQTL allele/frequency verified at [39122708873](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39122708873). Original eGFRcrea allele/frequency verified at [39124577872](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124577872). UMOD T frequency pQTL 0.00141033, original eGFRcrea GWAS 0.1168, external 1KG EUR LD 0.135189. Cause unresolved.

**Next scientific gate:** explicit LD-vs-Z compatibility diagnostics requiring effect-allele signed LD matrix and aligned SNP ordering. Examine sensitivity to outlier SNP exclusion **in a separate ephemeral diagnostic run only**, preserving existing result files and frozen evidence. Report conditional quality and robust nonconvergence, not causality.
