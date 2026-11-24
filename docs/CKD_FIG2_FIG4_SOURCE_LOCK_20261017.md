# CKD Figures 2 and 4: source-locked manuscript panels

**Source:** `/srv/is-analysis/results/ckd/integrated_evidence/CKD_9GENE_EVIDENCE_MATRIX.tsv`, SHA256 `32fa0475461b206c48295833acfc3075bfb24d89c322b5f9be554b1883a2f3d3`.

- **Figure 2:** EAS eGFR *single-anchor Wald MR* effect estimate, 95% CI calculated as beta ±1.96×SE. These are source trait/protein units and are not equivalent to a causal drug intervention, multivariable MR or longitudinal KoGES validation.
- **Figure 4:** count of kidney eQTL source hits in 686-person meta, 356 tubule and 303 glomerular resources, accompanied by candidate's expression-dominant kidney cell label. Counts across datasets differ in power/coverage and must not be ranked as direct tissue-effect strengths. The underlying resources may contain *significant-only* tables: 0 hits means **no hits documented in that resource**, not biological null.
- Both figures use precisely 9 genes, unmodified aggregate original-source values, and require the pinned SHA256. Invalid SE/P or malformed counts stop plotting. Missing fields remain not assessed, not zero.
- Outputs `CKD_FIG2_EAS_ANCHOR_WALD.svg`, `CKD_FIG4_KIDNEY_TISSUE_QTL.svg`, and `CKD_FIG2_FIG4_SOURCE.tsv`.
- Script: [plot_ckd_manuscript_fig2_fig4.py](../scripts/plot_ckd_manuscript_fig2_fig4.py).
- CI: [synthetic reproducibility](https://github.com/kimtk94/IS_Analysis_V3/actions/workflows/ckd-manuscript-fig3-ci.yml) and [read-only server rendering](https://github.com/kimtk94/IS_Analysis_V3/actions/workflows/thesis-server-readonly-inventory.yml).
- Full science gate still **unverified** until harmonized allele effects, sample overlap, tissue datasets and prior choice are independently audited.

These manuscript figures are descriptive representations of existing analyses, not newly replicated estimates.
