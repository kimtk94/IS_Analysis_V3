# IS P3 analysis bifurcation — single-signal coloc.abf versus multi-signal coloc.susie
Date: 2026-10-17

## Source-supported starting point

The read-only [Harness header audit #155](https://github.com/kimtk94/IS_Analysis_V3/issues/155) verified three existing EAS5 inputs:

- `IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz`: ~20.96 MB compressed, with source-provided `BBJ_allele_exact_Fasta37`, `EAS_AIS_allele_exact_Fasta37`, signed ALT effect and standard error/nominal p fields.
- `IS_EAS5_GENE_TISSUE_GWAS_SNP_OVERLAP_QC.tsv`: source-reported `GRCh37_FASTA_VERIFIED_SNP_rows`, `FASTA_REF_CONFLICT_rows`, locus-level non-palindromic overlap and `eligible_for_coloc`.
- `IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv`: `possible_single_signal_ABF_followup`, source REF and ALT-direction approval flags.

The headers are verified, but counts and row-level scientific eligibility have not yet been audited in this new stage. Prior EAS5 exact-SNP booleans are NOT independent SNP identifiers. **Never interpret a reported source FASTA flag as a fresh independent FASTA verification.**

## Two separate paths

| Path | Analysis | Minimum required information | LD requirement | Scientific limitation |
| --- | --- | --- | --- | --- |
| **P3A** | `coloc.abf` single-causal-signal hypothesis per trait | Complete region-wide QTL and stroke GWAS summary statistics on harmonized intersecting variants, effect/SE or other compatible p-value inputs, trait types and appropriate sample metadata; AF/MAF and case fraction where required; verified alleles and priors | Not mathematically required by the coloc.abf likelihood | One causal variant **per trait** within the region assumption; do not claim resolution of secondary signals |
| **P3B** | `susie_rss` plus `coloc.susie` for multiple signals | Compatible summary stats, sample-size assumptions, matched signed LD matrices with SNP order matching the effect vector, allele and genome-build validation, convergence / credible-set diagnostics | Required for reliable SuSiE-RSS inference | Source-cohort matched LD is preferred; mismatched proxy LD is **sensitivity-only** |

### P3A guarded acceptance criteria

1. Compare *actual* QTL/GWAS variant identifiers and signed alleles, not merely boolean `*_exact_SNP` fields; require a justified GRCh37 reference and liftover audit.
2. Use the complete, predeclared region or cis-window SNP set rather than selecting only the lowest-p SNPs, as cherry-picking can inflate shared-signal posterior support.
3. Record region-specific gene/tissue/trait type (`quant`/`cc`), sample sizes, case fraction for case-control traits, MAF semantics and effect allele direction.
4. Run P3A as `EXPLORATORY_SINGLE_SIGNAL` when full compatible summary statistics exist but evidence for one causal variant is limited. Do not silently equate ABF H4 to gene causality.
5. Check `PP.H3`, `PP.H4`, `PP.H4/(PP.H3+PP.H4)` with appropriate prior sensitivity and number/coverage of overlapping SNPs. State when small overlap, sparse QTL and uncertain AF prevent reliable inference.
6. If source credible sets imply multiple signals, or sensitivity is unstable, route to P3B; do not declare an H4 causal gene without functional/evidence convergence.
7. Keep `coloc_executed=0` until genuine computation results with logged versions, input hashes and test-reproducible posterior fields exist.

### P3B acceptance criteria

The QTL source metadata for QTD000131 GTEx artery aorta (N=387), QTD000171 GTEx brain cortex (N=205), QTD000434 ROSMAP DLPFC (N=560) is pinned to [EMBL-EBI eQTL Catalogue metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/0bc275374c210b4feaaecf217c1d612472f7ec87/tabix/tabix_ftp_paths.tsv). These dataset sizes do NOT establish ancestry, LD matching or individual-variant allele validity. Both primary QTL credible-set and logBF archive URLs can be located from the pinned metadata, but their contents have not been independently audited.

Before P3B execution, independently verify source cohort ancestry, matched signed LD, genome build, SNP overlap, variant order and SuSiE convergence. Never imply that LD is available just because SuSiE archive results are linked.

## Engineering status

New offline module `scripts/is/harness_archived_evidence.py` and tests `tests/test_is_harness_archived_evidence.py` provide the **read-only audit** for archived FASTA flags, unique SNP strings, actual source-row sums and ABF followup flags. CI passes; server deployment is **not yet verified**. Production `research-refresh` now marks `MODULE_NOT_INSTALLED_NOT_AUDITED` rather than fabricating new findings or crashing when the module is missing.

Next acceptance gate: safely deploy that audit module with the existing controlled three-script installation, inspect its archive-derived QC totals, and only then design P3A restricted full-summary-statistics analysis for eligible loci.

## 2026-10-18 archived FASTA/ABF module: actual full data acceptance

The previously blocked module `scripts/is/harness_archived_evidence.py` is now included in the safe GitHub Actions installer, retaining existing SSH identity, locked atomic installation and original source-file immutability. The installation workflow completed SUCCESS. Server `research-refresh` Issue [#158](https://github.com/kimtk94/IS_Analysis_V3/issues/158) completed PASS in **55.523 s** and emitted `is-archived-fasta-abf-audit-v1` with integrity PASS.

Observed archived source data (counts are NOT new independent FASTA verification):
- 98 region-gene-QTL-dataset group rows; 321,538 archived gene-assay SNP observations.
- 25,866 distinct variant strings.
- 293,520 archived rows with BOTH BBJ and EAS reported FASTA exact-match flags true; 23,468 distinct variant strings and 66 distinct genes represented.
- Overlap summary covers 393 groups and reports 423,051 BOTH GWAS exact-match rows and 376,334 dual nonpalindromic rows; those counts are sourced from archived/derived reports rather than fresh allele verification.
- Of 393 readiness groups, `possible_single_signal_ABF_followup` is True in **36** and False in 357. The 36 are **P3A review candidates, not completed coloc analyses**.
- Source-reported QTL ancestry and molecular-QTL matched LD remain unverified in all 393 reviewed groups. All 393 source readiness rows report GRCh37 REF/effect direction assertion flags true; this does not override independent FASTA and allele provenance requirements.

Detailed server-local evidence JSON is stored in `/srv/is-analysis/results/is/github_actions_runs/evidence/science/archived_fasta_abf_audit.json`. Scientific gate decisions must use its `abf_flagged_groups` list and original source files, not only the aggregate 36 number.

**Next gate:** review each of the 36 for full region-wide QTL/GWAS SNP coverage, independent genome build/REF-ALT harmonization, sample sizes, trait model (ischemic stroke case-control), case fractions, MAF/EAF semantics, single-signal assumptions and priors before any coloc.abf computation. Maintain `coloc_executed=0` until validated model output genuinely exists.
