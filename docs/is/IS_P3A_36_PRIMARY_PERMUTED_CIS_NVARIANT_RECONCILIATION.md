# IS P3A — independent original eQTL Catalogue permutation vs nominal full cis tested variants

**Date:** 2026-11-23 KST. **Scope:** MASTER DEGREE ischemic stroke P3A. **36 real gene×tissue QTL groups, 20 unique Ensembl genes, ONE chr10 G0015 GWAS locus.** Coloc remains 0/36 HOLD.

## Real primary original-source evidence

The [eQTL Catalogue official methods](https://www.ebi.ac.uk/eqtl/Methods/) document testing variants ±1Mb from the Ensembl gene start. The [qtlmap source usage](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) defines the window by original phenotype_pos metadata. Crucially, the [official original Column specification](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/07918a267afdc1c84b5444bea9bd902e7201728f/tabix/Columns.md) defines **permutation-file n_variants** as the number of actual variants tested within the cis region of a molecular trait.

[GitHub Actions #45036525615](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45036525615): **SUCCESS**, independently downloaded original public full permutation files for QTD000131, QTD000171 and QTD000434 and queried their exact original nominal associations on GRCh38 chr10:101837472–105353833 using one bounded indexed query per dataset. Independent earlier nominal source regional SHA-256 and original Ensembl105 phenotype positions, 36 source identities, per-SNP QTL allele-N evidence all pinned and checked. Focused 14/14 tests passed and full repository regression **364 executed, 363 PASS and 1 Rscript SKIP**. No full 11 GB nominal file downloads.

## Direct comparison across all 36 assays

| Original QTL source | Assays | Physical original nominal source rows (including repeated rsID) | De-duplicated original nominal cis tested variant×gene assay records | Original permutation n_variants exact-match groups | De-duplicated SNP records | De-duplicated INDEL records |
|---|---:|---:|---:|---:|---:|---:|
| QTD000131 GTEx artery aorta | 10 | 57,780 | **55,808** | **10/10** | 51,864 | 3,944 |
| QTD000171 GTEx brain cortex | 11 | 58,823 | **56,748** | **11/11** | 52,838 | 3,910 |
| QTD000434 ROSMAP brain DLPFC | 15 | 100,822 | **89,844** | **15/15** | 68,269 | 21,575 |
| **Total** | **36** | **217,425** | **202,400** | **36/36** | **172,971** | **29,429** |

15,025 physical rows duplicated the same original variant across multiple rsID representations and were deduplicated by (position, REF, ALT, molecular_trait_id) within each gene×tissue group. The 202,400 de-duplicated associations are **not** 202,400 distinct genomic SNPs; many variants are tested against multiple genes and tissues. The previously accepted strict nonpalindromic canonical SNP-only subset contains **148,344** group-assay rows; the other genuine original variant classes are retained as source metadata (palindromic SNP and INDEL). Do not conflate 29,429 unique original gene-assay INDEL tests with **44,454 original physical INDEL rows before rsID deduplication**.

Every original permutation n_variants agrees **exactly** with the corresponding full nominal cis unique all-variant count (36/36). All original nominal molecular_trait_id values for these 36 gene-expression phenotypes match the Ensembl gene IDs. No original association is outside the pinned reference Ensembl105 ±1Mb phenotype window; the three actual streamed nominal source SHA-256 hashes independently equal previously archived source evidence. This supports complete **publisher-reported per-gene tested cis variant counts** (not complete source VCF membership reconstruction beyond published summary stats).

### Publisher original permutation file provenance

| Source | Original permuted gzip bytes | SHA-256 of full published .permuted.tsv.gz | Original permutation phenotype rows |
|---|---:|---|---:|
| QTD000131 | 785,582 | 6366a7272124f512bb3a5e6ca5a55f62329492ecd9b5c5921755d9a93e061d10 | 21,888 |
| QTD000171 | 882,596 | fc76a8b0b3de5156767bf6871ddff2b3f3af7f382c31135d415a48ca13e7a735 | 24,675 |
| QTD000434 | 841,205 | adf13cbe9f784cd56a832700b710ffd6abe76a962a38918667a1d80c51413091 | 23,076 |

Original phenotype reference: Ensembl105 eQTL Catalogue published MD5 b999b21dcffeecbfe0bf886a5a91997b. The original nominal Tabix region SHA digests are pinned in the 36-group output JSON and match the independently archived original union source checks.

## Formal scientific HOLD still required

1. **QTL exposure regression model N:** original QTL genotype allele-number AN/2 was observed for all 148,344 strict nonpal SNP gene-assay rows; GTEx aorta N ranges 366–387, cortex 186–205, ROSMAP 560. The molecular expression phenotype and genotype sample intersection plus covariate model may reduce effective regression sample size; AN/2 cannot automatically be treated as model-specific effective N. One blanket dataset N across every original SNP is not justified.
2. **GIGASTROKE stroke GWAS model N / independence:** original 2022 Nature Supplementary Table1 verified discovery EAS AIS 19,032 cases / 237,242 controls and EUR AIS 62,100 cases / 1,234,808 controls. The exact per-SNP meta-analysis N assumptions for the publisher summary statistic file remain unapproved; BBJ contributes 17,651 cases and 159,610 controls to EAS discovery and therefore is not an independent replication source.
3. **Allele/LD and multi-variant modeling:** original reciprocal signed SNP matches EAS 89,597 and EUR 128,282 gene-assay records; 8 EUR ROSMAP discordant allele records, palindromic variants and original INDEL harmonization remain separately held. Valid locus definition, sample-overlap assumptions and sensitivity require explicit review before formal calculations.
4. **Coloc:** ABF formal eligible **0/36**, executed **0**. SuSiE formal eligible **0/36**, executed **0**; no posterior H3/H4 generated and no biological causal gene claims.

**Next:** audit actual QTL regression sample inclusion metadata and prespecify accepted ABF N/case fraction use-case; then score all 36 genuine groups together. No arbitrary N_eff imputation or false BBJ-independent replication.

Code: scripts/is/audit_p3a_36_primary_permuted_cis_nvariants.py; test: tests/test_is_p3a_36_primary_permuted_cis_nvariants.py; workflow: .github/workflows/is-p3a-36-primary-permuted-nvariants.yml. The artifact contains primary-public source SHA256 and 36-group test counts, not SNP β, SE, or private genotype/patient data.
