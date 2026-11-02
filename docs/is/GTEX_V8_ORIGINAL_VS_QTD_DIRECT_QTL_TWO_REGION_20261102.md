# IS P3A — Official GTEx v8 imported FastQTL versus EBI-reprocessed QTD000161 direct source test

Verified 2026-11-02 from actual GitHub Actions [#41446592395](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41446592395). **13/13 deterministic tests PASS; source test is actual primary-data remote tabix, not a mock dataset**.

## Comparators (not interchangeable)
- Official archived imported GTEx v8 gene eQTL, tissue Brain_Cerebellar_Hemisphere, URL `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz`, source checksum listed in public `md5sum.txt`: `c07a1731eb96aafbd7aed131cb97515a`.
- Official EBI eQTL Catalogue uniformly reprocessed GTEx QTS000015 / QTD000161 gene expression, URL `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000161/QTD000161.all.tsv.gz`. **2023 release-6 raw source bytes are not yet authenticated to the current source URL**, since same dataset ID may be reused in 2026 r8 pre-release.
- These are molecular eQTL **not** the GIGASTROKE GWAS phenotypes. QTL source comparison is a scientific independent gate from prior verified 6,273/6,273 GIGASTROKE GWAS source-to-Drive numeric match.

## Technical source authentication
- The official original-imported file's MD5 was confirmed in the official public manifest, but the **entire 3.4 GB source file was not downloaded or rehashed**, so manifest correspondence alone is insufficient to independently verify whole original bytes.
- First remote tabix run held on bundled pysam libcurl SSL CA error 77 ([#41444884247](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41444884247)).
- Next run with OS tabix held because the distribution binary lacked both HTTPS and FTP protocol handling ([#41445212004](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41445212004)).
- Explicit `CURL_CA_BUNDLE` and `SSL_CERT_FILE` successfully enabled pysam index reading for both remote files; however header-only Tabix comments yielded no declared field row ([#41445907723](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41445907723)).
- Bounded first BGZF compressed-block retrieval (HTTPS Range, no full download) recovered the actual first-row headers from both sources. Original has `variant,r2,pvalue,molecular_trait_object_id,molecular_trait_id,maf,gene_id,median_tpm,beta,se,an,ac,chromosome,position,ref,alt,type,rsid`. EBI reprocessed source also includes `molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se,...`. Actual interval extraction then succeeded.

## Actual primary-data results

Two deliberately **small** GRCh38 windows, 5kb each, were queried by indexed range from both sources. The all-gene results are:

| Comparison | chr4:79,759,900–79,764,900 | chr10:103,079,000–103,084,000 |
| --- | ---: | ---: |
| Original imported source rows | 171 | 624 |
| Reprocessed QTD source rows | 128 | 572 |
| Original distinct valid gene–SNP associations | 162 | 576 |
| Reprocessed distinct valid gene–SNP associations | 112 | 484 |
| **Exact common gene ID + GRCh38 position + REF/ALT** | **112** | **462** |
| Beta sign concordance | **74/112 = 66.07%** | **345/462 = 74.68%** |
| Near-identical beta estimates | 0 | 0 |
| Near-identical p-values | 0 | 2 |
| Median absolute beta difference | 0.064281985 | 0.0870759 |
| Median absolute SE difference | 0.01883545 | 0.0294511 |
| Median absolute difference in log10(p) | 0.210128629 | 0.367241763 |

**Combined 574 exact common gene–variant associations, 419 matching beta signs (73.0%).** This is evidence of different statistics in the distinct original and reprocessed sources. The high-level direction was also suggested by the legacy source comparison in the connected 2026-10-11 Drive reports, but that *legacy inherited “GTEx_V8” source* was not then proven to be GTEx official original FastQTL, so do not equate old and current comparison denominators.

**Scope caveat:** these 574 observations are **all genes in the two genomic windows**, NOT yet the four research candidates alone. The code has since been extended to disaggregate Ensembl stable IDs:
- FGF5: ENSG00000138675
- C4orf22 (current symbol CFAP299): ENSG00000197826
- CALHM2: ENSG00000138172
- NEURL (current symbol NEURL1): ENSG00000107954

See [CI comparison JSON](../../research/is/evidence/gtex_v8_vs_qtd000161_two_region_actual_20261102.json).

## Inference policy
- A source model beta/SE difference **does not mean either association is invalid**; differences can arise from normalization, covariates, cis SNP testing, genotype processing, sample inclusion and statistic definition. This audit did not isolate the cause.
- The 2 × 5 kb sample windows are **not** complete original cis variant universes for any of the four genes and cannot independently produce a locus-wide colocalization estimate or an unbiased source-wide agreement rate.
- Matching gene IDs and REF/ALT does not prove that original regression betas use exactly the same effect-allele conventions unless explicitly confirmed from each source's documented schema; source-specific allele encoding must remain an independent harmonization gate.
- These results DO NOT validate QTL cohort LD, source-specific historical N, GIGASTROKE AIS case fraction, independence of BBJ replication, MR, colocalization PP.H4, or causal target conclusions.
- The complete original GTEx v8 nominal associations data distribution via requester-pays was not accessed or charged. The publicly hosted imported source was accessed by bounded free Tabix.

## Immediate next science acceptance task
Run gene-specific filters and broaden to a **predefined cis tested variant universe in each original-vs-reprocessed source**, with proper release and full-source provenance. Repeat comparisons under exact variant+gene keys, then use selected single-source QTL with its documented N, LD, effect-allele, and study-level assumptions before colocalization. Until then overall scientific status **HOLD**, while bounded source-difference detection is **PASS**.
