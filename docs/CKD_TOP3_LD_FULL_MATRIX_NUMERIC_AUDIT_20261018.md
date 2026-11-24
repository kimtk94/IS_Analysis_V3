# CKD Top3 full-matrix numeric LD audit — 2026-10-18

**Scope:** Numerical integrity of existing 1KG EUR signed reference-LD matrices; strictly read-only. **Not** study-matched LD, full LD-positive-semidefinite (PSD) QC, or causal validation.

**Verified execution:** [GitHub Actions 39043729502](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39043729502), completed **success**. Source inspector: [scripts/audit_ckd_top3_ld_numeric_readonly.py](../scripts/audit_ckd_top3_ld_numeric_readonly.py).

| Metric | SDCCAG8 | INHBC | UMOD |
|---|---:|---:|---:|
| Matrix dimension | 3854 x 3854 | 2679 x 2679 | 3090 x 3090 |
| Maximum absolute asymmetry | 0.0 | 0.0 | 0.0 |
| Maximum |diag-1| | 0.0 | 0.0 | 0.0 |
| Maximum absolute correlation | 1.0 | 1.0 | 1.0 |
| Signed-matrix maximum asymmetry | 0.0 | 0.0 | 0.0 |
| Negative orientation-sign metadata | 3013 | 2052 | 2598 |
| Positive orientation-sign metadata | 841 | 627 | 492 |
| Result | FULL_SYMMETRY_SIGN_METADATA_PASS_PSD_NOT_TESTED | FULL_SYMMETRY_SIGN_METADATA_PASS_PSD_NOT_TESTED | FULL_SYMMETRY_SIGN_METADATA_PASS_PSD_NOT_TESTED |

The script reads all entries in each pre-existing dense float32 LD matrix and applies `LD_signed = LD * signs[:,None] * signs[None,:]`. This tests structural numeric consistency after the metadata sign transformation. It does **not** verify that effect allele vs reference allele coding was independently correct for every participant/source, only that the recorded sign vector is syntactically valid and preserves symmetry/diagonal.

### Outstanding scientific QC gates

1. **PSD/eigenspectrum** not checked. Eigenvalue analysis may be required if numerical pathologies could explain SuSiE instability; it must not be inferred from symmetric diag-one LD alone.
2. **Full signed beta alignment** has not been independently checked against raw molecular summary data and genotype dosage encodings.
3. **Summary vs reference-LD compatibility** remains an assumption. 1KG EUR `n=503` is only external reference LD and can diverge from GWAS/proteomic sample LD.
4. **ABF vs SuSiE SNP universes differ**, as recorded in [metadata audit](./CKD_TOP3_ORIGINAL_SNP_LD_METADATA_AUDIT_20261018.md). Do not interpret agreement/disagreement in PP.H4 as a pure algorithm comparison without SNP-coverage qualification.
5. **INHBC discordance** and **UMOD nonconvergence** remain substantive unresolved findings; no method flags are silently promoted to final scientific PASS.

**Paper-ready summary:** Existing LD matrices satisfy tested symmetry, finite-range, unit-diagonal and metadata-sign congruence checks. PSD, molecular LD matching and GWAS–pQTL effect-orientation claims remain pending; therefore Figure 3 is a descriptive analysis figure rather than experimentally validated causal molecular inference.
