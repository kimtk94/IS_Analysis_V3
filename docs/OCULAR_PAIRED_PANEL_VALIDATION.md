# Paired validation of the frozen ocular response panels

## Study and prespecified analysis

Use the seven-day, 100 nM dexamethasone study PRJNA729873 (SRP319815, Kathirvel et al., doi:10.1038/s41598-022-12443-7). ENA provides eight matched DEX and ethanol vehicle cell-line labels. The study includes four glucocorticoid responders and four nonresponders. Use all pairs together. The unit of inference is the matched donor-derived cell line. Do not treat genes or sequencing runs as independent replicates.

The primary panel is the exact 241-gene candidate archive from PR #333. The secondary panel is its 16 Reactome ECM organization genes. For both, fix the direction at the sign of the original 48-hour donor model. The three-strain GSE301517 results play no part in gene selection or weighting.

Use DEE2 STAR gene counts and match by stable Ensembl gene ID. Scale each library to counts per million using the sum of all assigned Ensembl gene counts and take log2(CPM + 1). Within each pair, multiply each gene's DEX minus vehicle log-expression by its frozen direction and average equally across mapped genes, including zero-count genes. Report the original and mapped panel sizes. A positive score means that the average change follows the 48-hour direction.

Before examining expression results, specify the primary population as complete pairs with DEE2 PASS or WARN for both arms. Report the per-arm QC flags and any pair exclusions. Summarize the median signed score and the count of positive, negative and tied pairs. Use the two-sided exact sign test across non-tied donor pairs for each panel and control the family of two tests with Holm adjustment. An exact zero is a tie. Report the complete-pair results including failed QC runs as a sensitivity analysis when their count data are usable. Do not select responder subgroups or genes by the validation outcomes.

This analysis estimates a panel-level response in cultured trabecular meshwork. Interpretation is specific to the source cohort, processing method and CPM normalization. It supplies no new per-gene FDR, causal time-course effect or glaucoma clinical endpoint.

## Reproduction

Download the frozen panel ZIP from GitHub Actions artifact 14912728421, the ENA read-run report and the DEE2 hsapiens project bundle SRP319815_NA.zip. The study-specific command checks the bundle digest and the sorted ENA field projection before joining the two sources by run and experiment. The output contains panel-level statistics and the availability of each pair.

```bash
python3 scripts/ocular_paired_panel_validation.py --panel frozen.zip --ena ena.tsv --bundle SRP319815_NA.zip --out paired/panel_summary.json
```

The paired-validation job in the ocular workflow uses the same command and publishes the JSON as an artifact. Source records are available from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJNA729873) and [DEE2](https://dee2.io/huge/hsapiens/SRP319815_NA.zip).
