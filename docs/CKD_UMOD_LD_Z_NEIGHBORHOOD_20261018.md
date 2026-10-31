# CKD UMOD signed-reference-LD versus Z neighborhood diagnostic — 2026-10-18

**Execution:** [GitHub Actions 39126192065](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39126192065), SUCCESS; no source/output modifications. Script: `scripts/audit_ckd_umod_ld_z_neighborhood.py`. Workflow: `.github/workflows/ckd-umod-ld-z-neighborhood.yml`.

Matrix: 3,090 x 3,090 external 1KG EUR reference LD, sign-transformed using existing `effect_vs_ld_major_sign`, aligned to Stage2C SNP rows.

| Anchor | GWAS Z | pQTL Z | |r| >= .8 neighbors | GWAS Z residual RMS | pQTL Z residual RMS |
|---|---:|---:|---:|---:|---:|
| rs11860566 (MAF-discordant) | -0.4 | +1.969926 | 81 | 0.206596 | 1.885733 |
| rs77924615 (max |GWAS Z|) | +31.784108 | -102.722975 | 0 | NA | NA |
| rs12922822 (max |pQTL Z|) | +21.817364 | -103.698773 | 27 | 4.592946 | 1.159458 |

**Algorithm:** For each anchor, extract its signed reference-LD row and select other SNPs with absolute LD correlation >=0.8. Compute residual RMS of `Z_neighbor - r_neighbor,anchor * Z_anchor`. This is a *descriptive local pattern screen*, not a DENTIST-S or formal conditional Z test, not a sampling distribution test, and not a direct SuSiE refit.

**Do not infer**: A large Z-score at a SNP with no other |r| >=0.8 neighbors is not automatically a genotyping or LD error. The reference sample (1KG EUR) may be small, and other causes include differences in cohort LD, multiple causal variants and low-frequency variants. A strong signed correlation itself is not proof of study-reference LD suitability.

**UMOD question:** The original pQTL rs11860566 T frequency=.00141033, original GWAS T frequency=.1168 and 1KG EUR T frequency=.135189 remain strikingly discordant. However rs11860566 has relatively weak local GWAS Z and is not the strongest signal in the 3,090-variant universe. No causal attribution for SuSiE 1000-iteration nonconvergence has been established.

Next evidence gate: check genomic-build metadata for source GENPOS vs ID-internal position and evaluate LD conditioning/isolated diagnostic reruns. For a valid sensitivity comparison, refit in ephemeral space with and without the SNP **using the same subset of LD rows/columns and corresponding Z-order**, and report iteration count/convergence rather than implying causal influence. No change to archived figure or published evidence until complete.
