# IS priority-four genes: verified EUR503 reference regional VCF acquisition

**Date:** 2026-10-11 KST. **Classification:** REFERENCE_GENO_COVERAGE_PASS_ONLY; no GTEx in-study LD, no multi-signal coloc, no causal gene verification.

## Download and source verification

Source: UCSC-hosted [1000 Genomes Phase 3 (20130502)](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/), hg19/GRCh37 `ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz` and corresponding chr10 file. Their `.tbi` index URLs responded HTTP 200. The EBI FTP index HEAD route attempted earlier returned HTTP 404; this **does not imply the study data are unavailable**, as alternate official mirrors are accessible.

We used `bcftools view -r REGION -S EUR503.samples -Oz` over HTTP random-access indexed VCF, **not** full chromosome download. Source 503 EUR sample names came from the existing project's verified `EUR503.samples`. The subset outputs have indexed compressed VCF and **503 unique samples**.

| Locus | Genes | GRCh37 extracted region | VCF records (bcftools index -n) | Approx size |
|---|---|---|---:|---|
| BBJ_IS_L001 | FGF5, C4orf22 | chr4:80,581,087–81,784,150 | 33,067 | 1.9 MiB |
| BBJ_IS_L002 | CALHM2, NEURL | chr10:104,739,152–106,059,690 | 34,744 | 2.0 MiB |

A `bcftools query` count of REF/ALT split allele keys is slightly greater than VCF row count due to multi-ALT records (chr4 33,245 keys, chr10 34,896 keys); do not confuse allele count with independent sites.

## Original SNP overlap

The original 30 archival ABF PASS results include **10 assays in these four genes** (6 chr4, 4 chr10). Using source SNP `variant_id=chr:pos:REF:ALT` in **GRCh37**, without risk of GRCh38 input confusion, the extracted regional VCF matched **100% of 10 assay SNPs** (1,694/1,694 where 1,694 denominator, 1,583/1,583 where 1,583 denominator; original assay-specific counts preserved). **No REF/ALT-swapped source keys** were detected. Genotype ref allele reference FASTA and sequence-level consistency, including all multi-allelic conventions, should still be audited before signed LD.

All source VCF SHA256 values and the detailed 10-assay SNP fraction are stored at:
`/srv/is-analysis/results/is/audits/is_priority4_eur503_allele_overlap_20261011_v1/` (TSV and JSON).

Extracted VCFs under:
`/srv/is-analysis/results/is/audits/is_priority4_1000g_eur_region_20261011_v1/`.

Executable independent overlap script:
`scripts/is/audit_is_priority4_eur503_region_overlap.py`.

## Inference boundary

These data give an ancestry-appropriate public **EUR reference panel** for *sensitivity*, not GTEx donor's true in-study LD. Both GWAS (BBJ EAS) and QTL (GTEx, mostly EUR) sides require independent LD/allele harmonization. No `coloc.susie` or validated multi-signal claim is made. Full IS search universe remains 2,225 genes.

**Next:** convert restricted EUR regional VCF to exact signed ALT-dosage matrix, compute QC missingness and allele frequencies, then compare to GTEx molecular QTL MAF and previously extracted BBJ EAS reference; pilot descriptive signed cross-ancestry LD only. Maintain full locus SNP variant order and validation gates.
