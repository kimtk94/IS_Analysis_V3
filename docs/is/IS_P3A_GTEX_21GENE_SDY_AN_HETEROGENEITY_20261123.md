# IS P3A — original GTEx 21 QTL sdY.est and heterogeneous AN/2 sensitivity

2026-11-23 KST. Extends original ROSMAP 15/15 QTL sdY.est source audit (PR #326) to all 21 remaining QTL gene×tissue pairs: 10 GTEx Aorta QTD000131 and 11 GTEx Cortex QTD000171. Genuine 36 gene×tissue groups, 20 genes, one chr10 locus. Formal coloc.abf/susie **0/36 HOLD**.

## Original publisher previously verified baseline

Original [source Actions #45038801844](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45038801844) independently scanned GTEx QTL all TSV rows in GRCh38 20-gene union chr10:101837472-105353833. Genuine canonical nonpalindromic SNP×gene-assay records: **44,693 Aorta**, **45,570 Cortex**, total 90,263 (not independent SNPs). Original whole union stream source SHA256:
- Aorta: **451ba8ba15493986f431b1751aff9ea512ab3dbc04ef41e38f95217d49cf3f53**, 241,786 QTL source rows.
- Cortex: **5964ba14a274f59247eebf88602da21706fdaea6c6e5df57990a3157f7854953**, 264,947 QTL source rows.

Original SNP AN/2 genotype allele-count proxies: Aorta **366–387**, Cortex **186–205**. All 21 gene assays have N heterogeneity; only approximately 98.86% / 99.04% SNP records use catalogue published nominal N=387 / 205. Therefore the published nominal cohort N is an **estimation scenario, not observed gene-expression regression N**. No guessed N_eff or sdY=1.

## New independent original source read and estimator

Read official eQTL Catalogue source files, one original public TBI index and one indexed HTTP 206 BGZF compressed range per study, never entire ~11GB eQTL source:
- https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000131/QTD000131.all.tsv.gz
- https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000171/QTD000171.all.tsv.gz

Verify original source region stream digest, original total region rows, all 21 source group IDs, per-gene SNP source counts, source MAF/SE, variant IDs, AN/2 bounds, duplicated associations and nonpalindromic alleles before estimation.

Exact original [coloc R sdY.est method](https://github.com/chr1swallace/coloc/blob/27776b1bb14e19d1542798ac7477fa378fc2fc0d/R/claudia.R): zero-intercept OLS of 2*N*MAF*(1-MAF) against 1/SE², then take square root of regression slope. Each source gene gets an estimate for nominal N assumption 387/205 (NOT certified regression sample N), and sensitivity diagnostics:
1. All original strict nonpal SNPs vs subset with variant AN/2 equal nominal N.
2. Source minimum observed AN/2 instead of nominal N, with exact square-root scaling.
3. Original MAF <0.1 vs >=0.1 regression estimate.
4. Sorted odd/even variant subset numerical sensitivity, **not independent** because LD.

A high estimator R² or split stability is not direct proof that original normalized expression sdY or regression N were measured. Do not auto-approve any coloc based solely on inferred N.

## Scientific HOLD and files

Source-only per-gene aggregate JSON report includes source SHA and counts; no original variant identifiers, individual SNP beta/SE/MAF rows or donor IDs. Original GIGASTROKE EAS+EUR beta=ln(OR) established across 725,707 original chr10 rows (PR #325); original coloc beta+SE case-control route without case fraction s established in PR #327. Next remaining work: normalized QTL model sample inclusion, complete cis GWAS-QTL GRCh37/38 harmonized variant coverage and prior sensitivity. BBJ is EAS GIGASTROKE discovery and not independent.

Files:
- scripts/is/audit_p3a_gtex_21gene_real_sdY_AN_sensitivity.py
- scripts/is/qtl_original_indexed_http_range.py
- tests/test_is_p3a_gtex_21gene_real_sdY_AN_sensitivity.py
- .github/workflows/is-p3a-gtex-21gene-real-sdy.yml
- docs/is/IS_P3A_GTEX_21GENE_SDY_AN_HETEROGENEITY_20261123.md

All **36/36 gene×tissue pairs formal Coloc HOLD**, none executed, no claimed H4 or causal signal.
