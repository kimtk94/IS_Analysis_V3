# IS P3A – input completeness and evidence acceptance (2026-10-18)

## Objective

Advance the 25 source-stat-positive groups from exploratory flags to an individually audited `coloc.abf` input without inventing genomic metadata. The 36 groups and 68,612 source-stat-valid archived rows were verified in [#166](https://github.com/kimtk94/IS_Analysis_V3/issues/166). None are yet an independently approved coloc analysis input.

## Current archived table limitations

Read-only header audit of `IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz` confirms QTL effect + SE and BBJ/EAS effect + nominal p-value. It has NO explicit BBJ/EAS effect SE, no GWAS per-study N or case/control fractions, and no independent provenance document for source-specific allele coding. Where only two-sided Wald p and signed beta exist, an approximate SE might be derived only under a validated statistical model and with proper handling of p-value underflow; it is **not automatically an approved observed SE**.

### Required evidence by component

| Evidence gate | Current evidence | Required acceptance |
| --- | --- | --- |
| QTL effect | QTL ALT beta, SE, p in archive | Confirm cohort's effect allele orientation, effective per-tissue N and QTL summary-stat source provenance |
| BBJ ischemic stroke GWAS | BBJ ALT beta and p | Locate official summary stats/metadata with harmonized effect/SE or validated derivation, exact phenotype definition, Ncases/Ncontrols and trait type |
| EAS-AIS replication GWAS | EAS-AIS ALT beta and p | Same as BBJ, including effective sample counts, overlap and independent alleles |
| Genomic build and allele | Archived GRCh37 mapping and FASTA flags | Recheck independently against versioned FASTA, explicit REF/ALT/strand, exclusions for conflicting and ambiguous SNPs |
| Cis-window | 98 archived full-cis gene/dataset groups, 25 source-stat-positive followup groups | Confirm region definition from QTL source and full intersection; no cherry-picking only significant SNPs |
| ABF model | `possible_single_signal_ABF_followup` flag | Document one-variant-per-trait assumption, priors and sensitivity range, dataset types, MAF and case fraction where relevant |
| Output | No verified colocalization run | Real version-logged PP.H0–PP.H4, H3/H4, overlapping SNP counts, prior sensitivity and model warnings |

## Automatic fail-closed behavior

The restricted server script `harness_archived_evidence.py` audits the original archived header, comparing the presence of BBJ/EAS SE, case-fraction and N columns. It does not infer missing values or execute coloc. `harness_research_jobs.py` includes the audit in the validated aggregate manifest and the existing 5-minute source watcher / hourly full audit cycle.

Until authoritative datasets are independently approved, `coloc_abf_eligible=0` and `coloc_abf_executed=0` remain the correct statements. Published SNP exact-match flags and archived FASTA statuses are **historical source assertions**, not independently repeated genomic validations.

## Next implementation after source metadata resolution

Build a restricted P3A dataset manifest keyed by `region_id, QTL_dataset, gene_id, GWAS_source` with `trait_type`, `N`, `Ncases`, `Ncontrols`, `case_fraction`, `effect_allele`, `variant_build`, `summary_stats_sha256`, `reference_fasta_sha256`, `reviewed_by`, and `review_date`. Do not put tokens, individual genotypes, private patient information or unreviewed metadata into GitHub issues. Run synthetic `coloc.abf` tests before enabling any real clinical GWAS calculation. Prior sensitivity and single-signal caveats must be reported with result tables.

This document specifies research gates; it does not certify that official BBJ/EAS GWAS metadata or consented patient genotype files are present.
