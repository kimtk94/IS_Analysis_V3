# Stroke GWAS replication-source audit: FinnGen versus GIGASTROKE

Audit date: 2026-11-21. Technical source QC PASS. Independent-replication claim NOT APPROVED. Primary LDSC effective-N HOLD.

## Audited official sources

* FinnGen R9 public manifest: https://storage.googleapis.com/finngen-public-data-r9/summary_stats/R9_manifest.tsv ; SHA256 13f3866a224e35471abc0290cba5ed151aad7a4e270ae70b3eacf3a92cf124a2 ; 2,272 entries
* FinnGen R10 public manifest: https://storage.googleapis.com/finngen-public-data-r10/summary_stats/R10_manifest.tsv ; SHA256 f5a97e7aad281ec1df05308d58b8a23587b2bc6568992be991c0a2221111ff63 ; 2,408 entries
* FinnGen R13 public manifest: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv ; SHA256 c9caafa9b98ee5ef050a705766bd1c451fc2154fb8e81b069022012a2f75ae20 ; 2,754 observed manifest entries
* Original GIGASTROKE corrected EUR 44-cohort Table 1: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/41586_2022_5492_MOESM2_ESM.xlsx ; 756,910 bytes; SHA256 7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4

CI verification: [initial manifest and source-label audit #44870079070](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44870079070) and [expanded exact original source case quantification #44870212469](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44870212469), both SUCCESS. Last source evidence [artifact #14693597353](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44870212469/artifacts/14693597353); 10/10 fixture tests PASS.

## Released GWAS endpoints are different from the clinical definition workbook

| FinnGen phenocode | R9 manifest | R10 manifest | R13 manifest | Clinical limitation |
| --- | --- | --- | --- | --- |
| I9_STR_EXH | ABSENT | ABSENT | ABSENT | Ischemic-leaning DF13 definition exists but has NO listed original GWAS result across these checked releases; includes I64 unspecified stroke |
| I9_STR | PRESENT | PRESENT | PRESENT | Includes I61 intracerebral hemorrhage; NOT clean AIS |
| I9_STR_EMBOLIC | PRESENT | PRESENT | PRESENT | ICD-10 I63.1; not necessarily TOAST cardioembolic subtype |
| I9_STR_EMBOLIC_WITH_FLIMMER | PRESENT | PRESENT | PRESENT | Requires AF for cases; controls explicitly require not AF |

| Reported phenotype-wide R9/R10/R13 cases | R9 cases/controls | R10 cases/controls | R13 cases/controls |
| --- | ---: | ---: | ---: |
| Embolic stroke | 1,373 / 342,673 | 1,485 / 374,631 | 1,919 / 450,257 |
| Embolic stroke with AF | 679 / 331,511 | 758 / 361,438 | 1,041 / 432,403 |

These are release-level endpoint counts, NOT per-variant effective sample sizes. The R9 and R10 manifest files have original-public-file SHA256 recorded here; the original publishers did not supply an independent published expected SHA for those files. Thus these SHA hashes are reproducibility anchors, not independent public checksum comparisons.

Official R13 manifest-linked compressed summary-statistics URLs for both embolic phenotypes responded HTTP 200 to an HTTP HEAD-only request (nominal lengths 797,290,613 and 790,024,136 bytes). No summary statistics file was downloaded, its header or variant data was NOT checked, and the FinnGen official [summary-data access procedure](https://www.finngen.fi/en/access_results) requires submitting a form.

## GIGASTROKE already used named FinnGen cohorts

GIGASTROKE corrected Supplementary Table 1 (EUR 44 cohorts) includes two explicit FinnGen-named rows:

| Original Table 1 row | Named GIGASTROKE cohort | Reported AIS cases | Reported controls |
| --- | --- | ---: | ---: |
| 46 | FINNGEN_HOSPITAL** | 8,097 | 162,569 |
| 47 | FINNGEN_Population** | 754 | 11,857 |
| Combined known AIS cases | | **8,851 out of 62,100 (14.25282%)** | Do not sum controls as proven unique |

FINRISK** (Table 1 row 27) is another Finland-related cohort candidate, but its sample counts are not included in the two FinnGen-labeled-row total.

The CE/LAA/SVD original case entries for these two FinnGen-named rows are NA, NOT 0. Do not use them to infer no participant overlap for these subtypes.

**Evidence classification:** overlapping named study origin is established. **Individual-level participant overlap between GIGASTROKE and FinnGen R9/R10/R13 is NOT measured.** The 8,851 cases are contributors to the historical GIGASTROKE original AIS source, not a verified count of individuals reused in FinnGen R13. Therefore FinnGen releases should NOT be described as independent external replication cohorts for this GIGASTROKE AIS analysis without source-backed participant-disjointness evidence. The original GIGASTROKE summary effects cannot be adjusted or made independent merely by subtracting this case-count proportion. CE comparisons are similarly subject to unquantified possible control/cohort overlap.

## Consequences and GO/HOLD

1. GIGASTROKE EUR AIS/LAA/CE/SVD whole-source integrity and allele QC: PASS. SNP-wise meta-analysis effective N, LDSC post-munge regression SNP set and independently calculated h2/intercept/rg: HOLD or NOT RUN.
2. FinnGen R9/R10/R13 clean AIS independent replication: **NOT ELIGIBLE** as checked, because I9_STR_EXH GWAS is absent and I9_STR is mixed hemorrhagic/ischemic; furthermore GIGASTROKE uses FinnGen-named cohorts.
3. FinnGen R13 embolic endpoints: **only an overlap-risk-qualified phenotype-specific direction/sensitivity candidate**, after published sample access, original file QC, GRCh38/37 harmonization, phenotype/power review, and overlap-risk disclosure; NOT TOAST-equivalent or independent replication.
4. Reassess genuinely nonoverlapping cohorts (e.g. independently validated East Asian stroke GWAS) for replication, taking ancestry/LD differences into account. Do not claim they exist in the needed format before source verification.
5. Cohortwise and pooled case/control harmonic effective N remain descriptions of published sample counts, not valid SNP-specific meta-analysis N estimates. No primary LDSC claim or new genetic-correlation value.
6. MR/colocalization can be assessed independently only under their own causal-assumption, cis-QTL, allele-harmonization, LD and replication gates. No MR/coloc performed here.

Implementation: scripts/ckm_stroke/audit_finngen_older_ais_manifest_overlap.py; tests/test_ckm_finngen_older_ais_manifest_overlap.py; .github/workflows/ckm-finngen-older-ais-manifest-overlap.yml. Draft PR #284 stays open and unmerged.
