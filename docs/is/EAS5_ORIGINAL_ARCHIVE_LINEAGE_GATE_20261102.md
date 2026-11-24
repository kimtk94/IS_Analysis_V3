# IS P3A: Original EAS5 archive & distinct BBJ canonical lineage (2026-11-02)

## Confirmed historical evidence

- [2026-10-18 original BBJ ZIP-vs-EAS5 comparison](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124421779): 13,435,541 original GWAS rows scanned; 25,441 / 25,441 archived SNP candidate keys numerically consistent for alleles, beta and p under historical tolerances.
- Original EAS5 comparison input on the research server: /srv/is-analysis/results/is/audits/is_eas5_98_egene_fullcis_gwas_overlap_20261011_v3_fastaanchored/IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz
- The old comparison did NOT emit original EAS5 archive SHA-256. The server-only path is not automatically accessible to GitHub Actions. An exact-name search of connected Drive did not return the archive; this does not prove absence from every possible nested folder.
- GitHub Actions [41432321496](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41432321496): eight synthetic positive/negative source-gate tests PASSED, and the actual data availability gate returned HOLD_ARCHIVE_NOT_AVAILABLE_ON_CI_RUNNER, no source hash. Successful CI != completed scientific comparison.

## Additional distinct source in Drive

[2026-10-11 full-cis summary report](https://drive.google.com/file/d/1sIko_lf94xHcENl-CUSYT_RmYbHjZN4N/view) and [source manifest](https://drive.google.com/file/d/1DmacR2aa0HJqIOKoKRS2htL4LEzviQ11/view) document a DIFFERENT four-gene GTEx-QTD000161 reprocessed cis-QTL/BBJ/GIGASTROKE crosswalk:

- FGF5, C4orf22, CALHM2, NEURL; 20,689 gene-SNP records; 9,810 distinct GRCh37 SNPs.
- Processed BBJ canonical GWAS: 13,435,531 rows; GIGASTROKE EAS-AIS canonical: 6,725,216 rows.
- Historical original BBJ ZIP scan: 13,435,541 rows, 10 MORE than the separate processed canonical BBJ scan. Filtered/removed row identities and transformations remain unverified. The 10-row difference does not alone prove a biological discrepancy.
- The Drive manifest hashes the separate reprocessed EBI QTL TSV and hg38ToHg19 chain, NOT the original EAS5 archive nor processed BBJ canonical GWAS.
- Original GTEx v8 QTL lineage and study-specific LD still unresolved; do NOT infer shared causal variants from this overlap.

## CI-only implementation

- Script: scripts/is/audit_eas5_original_archive_gate.py
- Tests: tests/test_is_eas5_original_archive_gate.py (8 cases)
- Workflow: .github/workflows/is-eas5-original-archive-provenance.yml
- With no original file, emits scientific HOLD and an aggregate JSON artifact.
- If the EXACT authorized archive is mounted in a permitted runner, the script can calculate its own SHA-256, check compressed TSV, allele-matching annotation and original 25,441 variant count. Even a supplied SHA match is NOT independent source authentication and does not permit automatic N, MR or coloc promotion.

## Required next scientific gates

1. Obtain the original EAS5 compressed file via an authorized route. Compute SHA-256 from actual bytes, obtain a separate independent source manifest, and validate source-to-archive linkage.
2. Reconcile the source-vs-canonical BBJ 10-row difference with exact filters and original REF/ALT plus effect-allele definitions.
3. Re-run the original ZIP-vs-EAS5 crosscheck on pinned inputs, then separately identify original GIGASTROKE EAS-AIS release and its actual SE/N/case fraction.
4. Only after harmonization and appropriate LD assumptions consider per-locus coloc.abf. BBJ and GIGASTROKE EAS are not independently validated replication cohorts.

**Scientific status: HOLD**. No new patient-level data, new MR, or new coloc in this action.
