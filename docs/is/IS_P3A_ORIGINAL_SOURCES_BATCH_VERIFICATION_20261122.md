# IS P3A original-data batch provenance and scientific readiness — 2026-11-22

**Scope:** MASTER DEGREE ischemic stroke only. Sales, SFDC, CRM, clinical/patient-level data and unrelated research workflows are excluded.

## What was actually recovered and independently checked

- [2026-10-18 original server aggregate, Issue #193](https://github.com/kimtk94/IS_Analysis_V3/issues/193): 2,225 candidate gene IDs, 393 gene×QTL×region readiness groups, 36 source-flagged P3A follow-up groups, 25 with some basic numeric stats and 11 without; **0 accepted new coloc.abf groups**. The **36 individual natural keys** (`region_id × QTL_dataset × gene_id`) are **not** preserved in the GitHub repository or indexed Google Drive resources as of this audit. The original server-local `science/archived_fasta_abf_audit.json` is inaccessible while the server is offline. Do **not** create 36 synthetic gene names, infer an independent gene count from 36 groups, or label 36 loci as negative.
- [2026-11-16 original full GIGASTROKE EUR/EAS Actions run #43858821798](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43858821798): direct full-source gzip transfer/structure QC: EUR **7,482,032** SNP rows, EAS **6,725,216** SNP rows. Earlier local compressed-file SHA-256 values were computed, but **publisher original MD5 had not yet been verified**, so the prior pipeline rightly said structural QC PASS and source scientific HOLD.
- [2026-11-16 original full EBI source manifest](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43858821981): primary `GCST90104540` (EUR) and `GCST90104545` (EAS) official metadata provide GRCh37 release names and publisher MD5 values, but the N in metadata is **not** a verified SNP-specific N or case fraction for archived P3A.
- [2026-11-16 original ImmuNexUT source Audit #43858822068](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43858822068): the genuine human East-Asian eQTL source had sampled original rows and source MD5 of small metadata files. Actual cell-type-specific N was not established, association SE was not provided in that sampled original source header, variant AF not available, complete original 43GB release digest not verified; therefore coloc ABF and SuSiE inputs **HOLD**.
- [2026-11-21 GIGASTROKE publisher supplement audit #44900451784](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44900451784): official publisher supplement had BBJ AIS contribution; public BBJ 2021 source indexed but not downloaded/independently matched; independent EAS stroke replication is **NOT ESTABLISHED**.
- The separate 15 real archived BBJ gene-tissue groups / 225 exploratory ABF mathematical scenarios are **not** a replacement for the 36 newer EAS5 P3A group records. See [15-panel report](IS_15_GENE_TISSUE_BBJ_COLOC_EXPLORATORY_BATCH_REVIEW_20261122.md).

## New GitHub Action: two large source files in parallel, joined scientific gate

`.github/workflows/is-p3a-all-source-md5-review.yml`:

1. Offline tests make sure a status cannot be upgraded using N, source headers, or fictitious group identifiers alone.
2. EUR and EAS Actions matrix jobs independently stream the full original `GRCh37.tsv.gz` compressed data directly from the exact accession-specific GWAS Catalog publisher endpoint, verify **live publisher YAML's MD5**, independently calculate full-file **MD5, SHA-256, bytes**, and compare with original 2026-11-16 previously captured SHA and file size. Downloaded GWAS SNP rows are neither printed nor added to GitHub.
3. Final CI job downloads original full-source QC reports from pinned Actions run #43858821798, original ImmuNexUT QTL source science gate #43858822068 and publisher BBJ supplement source audit #44900451784. It combines them with the **actual** 36-group aggregate and generates a reviewed machine-readable report. Only aggregate source statuses and blockers are exported; no raw variant rows.
4. If original server JSON with authentic 36 source keys becomes available, the workflow can read `docs/is/evidence/IS_P3A_ARCHIVED_FULL_GROUP_AUDIT.json` and produce group-level blocker rows without inventing IDs. This optional source must be reviewed before committing; the original privacy-sensitive source paths, raw GWAS/QTL rows and patient-level data must be stripped.
5. **Never set eligible_for_coloc_abf=true** without original source-specific GWAS/QTL matched signed effects, exact alleles/build, trait type, case/control proportions, QTL N and complete cis region evidence; multi-signal SuSiE additionally requires a signed, SNP-order-matched ancestry-appropriate LD matrix.

## Scientific status checklist

| Gate | Status pending current publisher rerun |
|---|---|
| Official file URL/GRCh37 source metadata | Source verified by preceding CI |
| Full EUR/EAS SNP/gzip structural QC | Technical PASS from preceding CI |
| **Live full original publisher MD5 vs file bytes** | New matrix job (see current Actions outcome) |
| Archived EAS5 BBJ/EAS allele and effect reconciliation | HOLD: original 36 group source rows inaccessible |
| Per-group genuine GWAS effect SE/case fraction/SNP N | HOLD |
| Per-group complete eQTL effects/SE/AF/true donor N | HOLD |
| Source cohort-matched signed LD for multi-signal coloc | HOLD |
| Exact 36 P3A group identities restored | HOLD |
| New independently approved formal coloc.abf groups | 0/36 |
| New P3A coloc.abf or SuSiE computations | NOT_RUN |
| BBJ-excluded independent EAS stroke replication | HOLD |

**Acceptance meaning:** A technical PASS for a compressed-file checksum **does not establish MR causality, locus colocalization or clinical replication**. Report source integrity separately from statistical readiness. Even EAS GWAS Catalog accession and BBJ genotype data may overlap cohorts, requiring independent ascertainment for replication. Do not rerun only FGF5 while ignoring the entire candidate pool.

## Next minimal science work, once exact archive is recovered

**One batch pass**: import the *real* 36 group `region_id / QTL_dataset / gene_id` entries; for each, compute full region-wide variant coverage and independently signed matched variant counts; check original QTL N and GWAS Ncase/Ncontrol/SE/orientation; only source-approved groups enter `coloc.abf` R `matrix` jobs; publish every group including BLOCKED/NOT_TESTED; then compare PP.H3, PP.H4, H4/(H3+H4), prior sensitivity and correlated QTL tissue evidence in one gene-level matrix with explicit statistical denominators.
