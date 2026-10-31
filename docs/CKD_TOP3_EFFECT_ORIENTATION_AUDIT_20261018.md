# CKD TOP3 Stage2B→Stage2C effect orientation consistency — 2026-10-18

**Run:** [GitHub Actions 39111270087](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39111270087), SUCCESS.

This read-only server audit compared existing harmonized Stage2B GWAS–pQTL variant records with Stage2C SuSiE input records and 1000 Genomes EUR reference LD metadata. No input or LD files were changed.

| Criterion | SDCCAG8 | INHBC | UMOD |
|---|---:|---:|---:|
| LD-ordered variants checked | 3,854 | 2,679 | 3,090 |
| Position errors | 0 | 0 | 0 |
| Stage2B effect vs LD metadata allele errors | 0 | 0 | 0 |
| Effect/other allele vs reference REF/ALT errors | 0 | 0 | 0 |
| Reference-major allele sign derivation errors | 0 | 0 | 0 |
| Stage2B→Stage2C pQTL beta changes | 0 | 0 | 0 |
| Stage2B→Stage2C aligned eGFR beta changes | 0 | 0 | 0 |

**Pass definition:** SOURCE_CONSISTENCY_PASS_RAW_UNVERIFIED. The effect allele was checked against reference REF/ALT and the reference-major allele using recorded reference ALT frequencies. The signed LD orientation flags matched these derivations. Aligned pQTL and eGFR betas were unchanged in the Stage2B→Stage2C subset.

**Limitations:** This does *not* independently validate raw outcome/pQTL consortium summary files or recover genotype-level allele dosages. Nor does it establish that 1KG EUR n=503 LD is appropriate for each original study. Technical consistency should not be promoted to independent external replication or a causal gene verdict. In particular, INHBC ABF–SuSiE disagreement and UMOD SuSiE nonconvergence remain open.

Relevant script: [scripts/audit_ckd_top3_effect_orientation.py](../scripts/audit_ckd_top3_effect_orientation.py). Previous [LD PSD analysis](./CKD_TOP3_PSD_METHOD_AND_DECISION_GATE_20261018.md) found minimum eigenvalues only slightly negative within numerical tolerance, which does not explain or eliminate model-specific failure.

**Next scientific gate:** assess phenotype/pQTL source allele harmonization against *original* input and evaluate LD–summary-statistic consistency; do not rerun models or overwrite live outputs until a justified reason is documented.
