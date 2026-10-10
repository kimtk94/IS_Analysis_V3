# G0022 Japanese OASIS — EAS AIS allele identity and SNP coverage audit

**Date:** 2026-10-10, Asia/Seoul.
**Gate:** five real OASIS gene/cell SNP plots reviewed; rs11066015 identity matches (5/5); three remaining AIS credible-set variants not observed in these plots; OASIS signed beta harmonization and full GWAS–QTL coloc NOT PERFORMED.

## 1. Direct input provenance

GIGASTROKE East Asian acute ischemic stroke GCST90104545, GRCh37: G0022_FOUR_CS_SIX_GWAS_EFFECTS.tsv from existing audit script audit_g0022_rs671_stroke_subtypes.py, itself checking matched GRCh37 source variant ID, reference alleles, original GWAS effect allele, sign flip and normalized ALT beta/frequency. HECTD4 rs77768175 is a relevant edge case: its original source effect allele is REF A and validated output was transformed to ALT G. Requiring the original GWAS effect allele to always be ALT would be a scientific/engineering error; our input validation accepts either valid original allele but requires source-audited ALT-normalized beta.

The four candidate IDs use **prior source mappings**, not an independently re-run liftOver:
- rs11066015: GRCh37 12:112168009:G:A; prior GRCh38 chr12:111730205:G:A; AIS CS PIP 0.27247
- rs671: GRCh37 12:112241766:G:A; prior GRCh38 chr12:111803962:G:A; AIS CS PIP 0.37484
- rs11066132: GRCh37 12:112468206:C:T; prior GRCh38 chr12:112030402:C:T; AIS CS PIP 0.29776
- rs77768175: GRCh37 12:112736118:A:G; prior GRCh38 chr12:112298314:A:G; AIS CS PIP about 0.03778

## 2. Real OASIS browser source results

Public Japan Omics gene pages embed Bokeh SNP-level data including non-significant nominal eQTLs. Five cached HTML pages yielded these results:

| gene | OASIS cell type | plotted SNPs | SNPs with p>0.05 | rs11066015 eQTL p | QTL beta and SE, still unharmonized |
|---|---|---:|---:|---:|---|
| ALDH2 | Monocyte L1 | 2254 | 1435 | 0.000913 | -0.222 ± 0.0661 |
| ALDH2 | Activated B L2 | 2266 | 2017 | 0.845 | 0.0211 ± 0.108 |
| BRAP | Monocyte L1 | 2140 | 1623 | 0.303 | -0.0831 ± 0.0804 |
| BRAP | Activated B L2 | 2162 | 2150 | 0.810 | 0.0284 ± 0.118 |
| RPH3A | Monocyte L1 | 2968 | 2235 | 7.62e-7 | -0.427 ± 0.0838 |
| RPH3A | Activated B L2 | NO BOKEH DATA | NOT ASSESSED | NOT ASSESSED | NOT ASSESSED |

All five available plots span roughly 2 Mb on chromosome 12, consistent with the OASIS cis window of ±1 Mb from each gene TSS described in the paper. This verifies displayed SNPs with nonsignificant association are included, but not that **all tested SNPs** were plotted or that variant QC is equivalent across genes and cell types.

All five available plots contain exact rs11066015 / GRCh38 chr12:111730205:G:A with the same REF/ALT pair as the AIS GRCh37 SNP. AIS ALT A beta=-0.149, SE=0.0174, ALT EAF=0.2319, p=9.515e-18. OASIS MAF for rs11066015 was 0.297 in Monocyte L1 and 0.283 in activated B L2. **OASIS MAF must not be treated as oriented ALT EAF without checking minor/ALT status and source ancestry composition.**

## 3. Exact browser missingness is NOT biological absence

Each other GWAS AIS 95% credible-set SNP is located **within** the observed 2Mb plotted coordinate range on all five source pages but is not found by exact rsID:
- rs671, prior GRCh38 position 111803962: closest plotted SNP **214bp away**
- rs11066132, prior GRCh38 position 112030402: closest plotted SNP **211bp away**
- rs77768175, prior GRCh38 position 112298314: closest plotted SNP **35bp away**

These exact nearest-SNP distances were repeated across all five pages. This makes an out-of-cis-window explanation implausible for the browser omissions, while original source variant/genotype QC versus web visualization subsetting remains undetermined. Source study uses MAF selection; do NOT declare rs671 unavailable in the full E-GEAD-1054 nominal dataset without actually examining the relevant original TAR member.

## 4. Publication and effect-allele gating

- Original AIS GWAS allele normalization and ALT effect validated in preexisting G0022 pipeline.
- rs11066015 rsID / chrom / GRCh38 position / displayed G>A SNP allele identity observed in five source plots. Prior build liftOver map is NOT independently recalculated during this audit.
- The actual **OASIS beta effect-allele semantics were not independently verified** from original study matrix documentation; neither same G>A SNP string nor the negative beta proves signs can be compared.
- The OASIS Monocyte p-value being lower for RPH3A than ALDH2 does not show RPH3A is more likely to cause stroke; tissue-dependent power, expression and LD differ. A formal heterogeneity test has not been computed.
- No valid GWAS/OASIS full-cis SNP universe, source signed LD, conditional eQTL fine-map, coloc.susie, Mendelian randomization or causal mediation yet.
- Three missing Bokeh entries and RPH3A B L2 no embedded Bokeh payload are not evidence of null biology.

## 5. Reproducibility and next move

New source-preserving code:
- scripts/is/audit_g0022_oasis_gwas_coverage_and_alleles.py
- tests/test_is_g0022_oasis_gwas_coverage_and_alleles.py

Read-only source snapshots:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_browser_gene_v1/OASIS_{GENE}_{CELL}_original.html

Machine-readable output folder:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_ais_allele_coverage_v1
  - G0022_OASIS_FOUR_CS_ALLELE_COVERAGE.tsv — 24 gene-cell–SNP comparisons, distinct observed, not-plotted and no-source states, nearest SNP distances
  - G0022_OASIS_GENE_CELL_SOURCE_COVERAGE.tsv — 6 SHA256 source snapshots
  - G0022_OASIS_ALLELE_COVERAGE_SUMMARY.json — fail-closed full-coloc and effect orientation gates

Next highest-value task: obtain targeted full nominal E-GEAD-1054 cis rows including non-significant SNP–gene observations for ALDH2, BRAP, RPH3A and other competing candidates; confirm all four CS alleles and sample/MAF/INFO, source effect allele semantics and matched signed LD, then run proper multi-signal QTL fine-mapping/coloc. Preserve broader 2,225-gene IS candidate universe.

Official study: https://www.nature.com/articles/s41588-025-02266-3
Original source: https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-1000/E-GEAD-1054/
Browser example: https://japan-omics.jp/gene/OASIS-Mono-L1?input_value=ENSG00000111275

**No 39GB source download, original GWAS modifications or unrelated production-service modifications were made.**
