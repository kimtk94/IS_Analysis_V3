# IS P3A — 15 ROSMAP gene cis-QTL original single-tabix regional audit

Date: 2026-11-23 KST. Scope: ischemic stroke research only; all 36 original P3A assays (15 ROSMAP QTD000434, 21 GTEx, 20 distinct genes, chr10 G0015). **All formal coloc.abf and coloc.susie remain HOLD.**

## Scientific purpose

The previous original [QTD000434 200-row preview](https://github.com/kimtk94/IS_Analysis_V3/pull/322) showed that actual source AN/2=560 can be observed in the original QTL file, matching the official eQTL Catalogue method example. Those 200 arbitrary source rows **are not the 15 true chr10 candidate genes**. This workflow directly tests all 15 candidate gene IDs together.

## Primary source and exact gene-window design

1. Official [eQTL Catalogue original phenotype metadata version v23.04.1 (Ensembl 105)](https://zenodo.org/records/7808390), published April 7, 2023. The publisher provides gene_counts_Ensembl_105_phenotype_metadata.tsv.gz (published MD5 **b999b21dcffeecbfe0bf886a5a91997b**, approximately 1.9 MB) with GRCh38 chromosome and exact phenotype_pos (the actual default cis-centre used by qtlmap). Download once and verify MD5. It is a historical **2023 source**, not certified identical to the Release 7 normalized expression or exact run configuration.
2. Extract the **15 exact original gene IDs** from the versioned 36-assay manifest, require one and only one phenotype annotation each, chromosome 10, and phenotype_pos within the source annotation gene coordinates.
3. Build **one union genomic interval** from minimum(15 phenotype_pos) minus 1,000,000 bases through maximum(15 phenotype_pos) plus 1,000,000 bases. The source [qtlmap usage](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) documents 1 Mb as the default half-window; the exact Release 7 qtlmap run option **has not been independently recovered** and therefore full-cis completeness remains scientifically provisional.
4. Exactly **one** remote indexed tabix query of the real 4 GB source, restricted to chr10 union interval, with no retry:
   https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz
   No downloading the entire file or performing 15 separate requests. [Official EMBL-EBI guidance](https://www.ebi.ac.uk/eqtl/Data_access/) warns frequent incomplete requests may trigger a firewall/DoS block, so keep requests minimal.
5. Stream all returned rows (up to a conservative resource guard) and **retain only in-memory variant QC** for the 15 target genes. Match molecular_trait_id, gene_id, GRCh38 position, REF/ALT, variant ID, numeric AN/MAF/beta/SE, p-value; safely collapse repeated rsIDs by exact variant and verify consistent effects; report different ALT alleles at the same site as multiallelic positions.
6. Save **aggregate results only**, per-gene source row count, unique variant count, duplicate count, multiallelic positions, AN range and mode, MAF range, reference-only legacy variant count mismatch. Do not emit original per-SNP beta values, sample IDs, genotypes, or rows.
7. Scientific acceptance is separate from technical pass. Report how many out of 15 genes have any original SNP-level observations; this does not prove that each Release 7 source cis region (which may have used run-specific parameters) is complete. A missing gene is NOT imputed from previous legacy results. All 36 QTL×stroke groups remain HOLD until GWAS binary beta/source checks and QTL-N/sdY uncertainty are addressed.

## Automatic regression gates

- Script: scripts/is/audit_p3a_rosmap_15_gene_single_tabix_qtl.py
- Tests: tests/test_is_p3a_rosmap_15_gene_single_tabix_qtl.py (synthetic biological and negative cases, no source requests)
- Workflow: .github/workflows/is-p3a-rosmap-15genes-single-tabix.yml
- Output artifact: aggregate 15 gene results, original publisher MD5 and SHA256, one combined query region, 36-group coloc HOLD status; zero SNP-level rows.
- Run at most one tabix process with single region; refuse more than 20 Mb union, unexpected source schema, wrong chromosome/alleles, non-finite numeric rows, inconsistent rsID duplicates, incorrect original phenotype metadata MD5, too many rows or gene variants.
- The eventual review must distinguish 2023 annotation compatibility from a proven Release 7 processing pipeline. **No arbitrary sdY=1 and no approved effective N without evidence.**

## Next step

If actual 15-gene source variant data are present, test per-gene AN stability and original-stat completeness against the observed archive, and separately validate GIGASTROKE AIS GWAS effect allele, case fractions, SE, study overlap and log(OR) units. Only then may formal coloc be conditionally eligible. The official publisher method permits a quantitative N/MAF-based *estimated sdY* route but does not eliminate uncertainty in its original regression N.

## Original TBI source transport implementation note (2026-11-23)

The initial direct remote HTSlib tabix HTTPS/FTP invocation failed to open the remote publisher source on GitHub-hosted runners. The validated alternative downloads **only the publisher 2.14 MB public Tabix TBI index** and **one contiguous HTTP 206 BGZF byte range**, computed offline from the index (chr10 101854259-105353833; 8,434,192 compressed bytes). Original TBI v1 has 22 contigs and a special pseudo-bin **37450**, which records metadata rather than ordinary virtual offsets. It must be excluded from BGZF chunk lists. The temporary index debug artifact used to diagnose this is removed from the final automated artifact; only aggregate per-gene scientific QC is retained.
