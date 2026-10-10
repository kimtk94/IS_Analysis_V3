# IS G0022: directly observed OASIS nominal eQTL through Japan Omics Bokeh payload

**2026-10-10 KST. Status: Partial source-backed marginal sc-eQTL evidence, NOT valid locus-wide colocalization.**

The Japan Omics Browser OASIS gene-centric pages embed genuine Bokeh ColumnDataSource arrays (some columns base64 little-endian IEEE754). HTML tables may be empty while the plot serializes SNP-level molecular-QTL data. A source-preserving Python extractor decoded variant_id_hg38, rsid, gene_name, pval_nominal, effect_size, effect_size_SE, maf, cell_type and cluster_res without bulk 39GB E-GEAD-1054 download.

## Real queried source rows

| Gene | cell | SNP association rows | p>0.05 records | rs11066015 p |
|---|---|---:|---:|---:|
| ALDH2 | Mono-L1 | 2254 | 1435 | 0.000913 |
| ALDH2 | B_Activated-L2 | 2266 | 2017 | 0.845 |
| BRAP | Mono-L1 | 2140 | 1623 | 0.303 |
| BRAP | B_Activated-L2 | 2162 | 2150 | 0.81 |
| RPH3A | Mono-L1 | 2968 | 2235 | 7.62e-7 |
| RPH3A | B_Activated-L2 | 0 *browser Bokeh data not provided* | NOT ASSESSED | NOT ASSESSED |

**Important:** rs11066015 is only one of the four GWAS 95% CS variants. The other 3 high-LD variants were not retrieved from these pages. This can reflect source MAF filtering, selection, genotype imputation, plot subsetting or ID matching: absence on a browser plot is NOT proof of missing genotypes or biological null. Some OASIS gene-cell web pages do not embed QTL arrays. The output retains status separately rather than claiming true 0 associations. Different immune-cell type stratifications (Mono-L1 vs B_Activated-L2) should not be pooled or considered matched effective N.

Gene-dependent p-values cannot establish shared causal GWAS/QTL variant, relative gene causality, enzyme activity or drug action. Existing ImmuNexUT significant marginal gene associations and JCTF PIP-filtered molecular evidence are separate source evidence, not coherent coloc priors. For valid coloc, obtain unbiased all-tested cis SNP–gene rows in the EXACT per-gene cell group, harmonized alleles plus matched signed LD and GWAS coverage, check per-cell sampling and covariates, run fine-mapping and examine competing ALDH2/BRAP/RPH3A genes.

## Reproducibility

Code: scripts/is/audit_g0022_oasis_browser_qtl.py; tests/test_is_g0022_oasis_browser_qtl.py.
Output: /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_browser_gene_v1/G0022_OASIS_BROWSER_AUDIT.json
All six public original HTML pages snapshotted with SHA256 in output JSON. No private human genotypes, canonical AIS/CKD files or trading server services altered.

Official browser source:
https://japan-omics.jp/gene/OASIS-Mono-L1?input_value=ENSG00000111275
https://japan-omics.jp/gene/OASIS-Mono-L1?input_value=ENSG00000089234
https://japan-omics.jp/gene/OASIS-Mono-L1?input_value=ENSG00000089169

Scientific status: SOURCE_BOKEH_ASSOCIATIONS_PARTIAL; FOUR_CS_SOURCE_COMPLETENESS_UNVERIFIED; FULL_COLOC_NOT_PERFORMED.
