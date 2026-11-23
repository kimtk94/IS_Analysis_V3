# IS P3A — Original eQTL Catalogue Coloc AN/MAF route and ROSMAP source preview

**MASTER DEGREE / ischemic stroke; 2026-11-23 KST.** Prior PR #321 confirmed ROSMAP QTD000434 nominal N=560 in R6/R7/R8. Exact original per-gene regression N and normalized phenotype sdY remain scientifically unproven. All 36 original groups HOLD.

## New original-publisher Coloc method evidence

The [original eQTL Catalogue Coloc tutorial](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/tabix_use_case.md) and independent [official batch Coloc R implementation](https://github.com/eQTL-Catalogue/colocalisation/blob/master/bin/run_coloc_batch.R) construct a *quantitative QTL* dataset using association beta, SE squared as varbeta, study-group MAF, and the first row's allele number AN divided by 2 as the supplied N; direct measured sdY is not an explicit input. Coloc can then estimate a quantitative sdY from variance of coefficients, MAF and N: see [original coloc source](https://github.com/chr1swallace/coloc/blob/main/R/claudia.R).

The original [publisher QTL column definitions](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md) establish: AN is total allele count; for autosomal variants AN/2 is *usually* the sample size; MAF is minor allele frequency in the QTL mapping context; beta is regression effect and SE its standard error; ALT is the effect allele (GRCh38 position and reference).

### Crucial interpretation

This identifies a **legitimate alternative to obtaining the original normalized expression matrix**, but does NOT establish the exact sample size of any ROSMAP gene-specific association. Variation in AN across variants, missing genotypes, conditional genotype/sample intersection, and source release mismatches must be checked. The original example selects the **first** AN/2 and is not a blanket certification of this parameter for our datasets.

The official tutorial's GWAS example is *mean platelet volume*, a **quantitative** outcome. Our ischemic stroke outcome is **case-control**, requiring binary GWAS type cc, verified case fraction (cases divided by cases+controls), exact beta/log-OR scaling and allele harmonization. Do not label original ROSMAP sdY as measured, set sdY=1, invent effective N, or calculate formal posterior H4 for an unapproved pair.

## Source retrieval design

Workflow fetches and hashes original publisher documentation and real QTD000434 output URL from the official eQTL FTP manifest. It requests only a **bounded first 128 KiB HTTP byte range** from the official QTD000434 gene-expression full-cis .all.tsv.gz source. Requires actual HTTP 206 (range honored), trusted official HTTPS host, gzip signature, all actual publisher QTL columns, and numeric source AN/MAF/beta/SE. It emits only aggregate preview counts and AN/2 plus MAF ranges, no original per-SNP beta or sample identifier. Preview rows are **not** evidence for completeness of the 15 original chr10 ROSMAP genes.

- Script: scripts/is/audit_p3a_official_coloc_ana_qtl_source_probe.py
- Tests: tests/test_is_p3a_official_coloc_ana_qtl_source_probe.py
- CI: .github/workflows/is-p3a-official-coloc-ana-qtl-source-probe.yml
- Artifact: 36-group scientific HOLD decisions plus tiny publisher-source preview metadata.

## Scientific next gates

1. Retrieve original **indexed full-cis SNP-level QTL source** for all 15 authentic ROSMAP genes (and 21 GTEx assays), with locus build, complete declared window, measured beta/SE/MAF/AN, allele and genotype missingness audit.
2. Determine whether observed AN/2 is a defensible quantitative regression N for these 15 genes; estimate sdY from varbeta/MAF/N only as a **method-labeled estimator** with sensitivity and diagnostics, not as original normalized source sdY.
3. Independently verify GIGASTROKE EAS/EUR AIS original case-control beta scale, case fraction, sample counts, replication independence and complete intersected alleles, then perform eligible formal coloc.abf and sensitivity, never earlier.

**Formal Coloc remains 0/36** until original data and statistical evidence gates pass. Technical CI success is not biological causality.
