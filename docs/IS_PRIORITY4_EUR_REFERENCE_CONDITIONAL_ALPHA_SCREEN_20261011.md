# IS priority four: EUR503 reference-lead conditional z sensitivity

**2026-10-11 KST.** Exploratory diagnostic only; **not** actual GWAS conditional association, independent causal signal discovery, fine-mapping or coloc.susie. No original PP.H4 or 2,225 candidate universe changed.

Using BBJ IS `z=beta_ALT/se` with exact GWAS–EBI cis-match and EUR503 signed ALT-dosage reference correlation `r(j,lead)`, we calculated the heuristic lead-conditioned residual `z_approx(j)=(z_j-r*z_lead)/sqrt(1-r²)`. At shrinkage alpha, off-diagonal `r` was replaced by `(1-alpha)r`. The lead SNP itself was excluded for every alpha. **This formula should not be treated as formally calibrated conditional GWAS**, because the reference LD sample differs from the BBJ EAS study and LD uncertainty, sample overlap, and effect estimate correlations were not fully modeled.

| Locus | Evaluable SNPs per alpha | alpha=0: >=5.45 | 0.01 | 0.05 | 0.1 | 0.2 |
|---|---:|---:|---:|---:|---:|---:|
| chr4 FGF5/C4orf22 | 2,824 | 0 | 0 | 0 | 0 | 0 |
| chr10 CALHM2/NEURL | 2,814 | 1 | 1 | 1 | 1 | 0 |

The absolute-z threshold of **5.45** is an exploratory screen roughly near the two-sided nominal p=5e-8 threshold under a standard normal model. The one chr10 candidate is **GRCh37 10:105694487:G:T**, with heuristic z=5.498 at alpha=0.05 (EUR reference lead signed r=0.1092). It fails the same criterion at alpha=0.2 (maximum chr10 absolute z=5.424). Thus it is **shrinkage-sensitive**, and cannot be promoted as a secondary independent association.

The earlier 20 fixed 250kb EUR reference LD tiles and shrinkage spectral QC remain intact. A one-lead approximation is not block-wise multivariate conditional regression. Fine-mapping posterior inclusion probabilities and credible sets have **not** been computed. No SuSiE-coloc was run, and no original GTEx QTL donor LD matrix is authenticated.

Outputs: `/srv/is-analysis/results/is/audits/is_priority4_eur_ref_lead_conditional_20261011_v2/`, comprising per-SNP/per-alpha diagnostic z scores, scenario QC and JSON manifest.

Code: `scripts/is/audit_is_priority4_eur_ld_conditional_shrinkage.py`, test `tests/test_is_priority4_eur_ld_conditional_shrinkage.py`.

Next gate: investigate chr10 10:105694487:G:T against original BBJ GWAS conditional analysis or appropriately ancestry-matched EAS LD, and check whether it is a locus-specific independent signal; if cohort-correct matched LD and full QTL test universe can be obtained, formal conditional fine-mapping becomes eligible. Keep EBI QTL source separate from inherited historical QTL and note QTL effective sample N is not independently authenticated.
