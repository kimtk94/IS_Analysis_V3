# External subtype validation: source decision update (2026-11-28)

Status: HOLD_EXTERNAL_SUBTYPE_DATA. No external effect estimates received, no outreach sent, no individual-level data accessed.

## Prior work reused
PR #349 already documented the TPMI limitations and prepared a fixed 24-marker request and GENESIS-K outreach draft. Reuse the immutable source rather than create a second marker family:
- Source commit dea5733170cebab70946721e4d599ddff91c13bd
- docs/publication/external_validation_20261125/marker_request.tsv (Git blob 3035e008ddc29133c5f72008977c0dfe8a07551e)
- docs/publication/external_validation_20261125/OUTREACH_DRAFT.md
- docs/publication/IS_EXTERNAL_LAS_VALIDATION_FEASIBILITY_20261125.md
- docs/publication/IS_TPMI_EXTERNAL_VALIDATION_DECISION_20261125.md

## New verified evidence
- Actions 46126983365: GWAS Catalog v2 returns two Taiwan studies, GCST90866312 AIS and GCST90866311 early-onset AIS; neither provides full summary statistics through the catalog. Singapore PMID 42448752 returns an explicit total of zero. Neither is a usable independent LAS GWAS.
- Actions 46127093339: Taiwan PMC13239414 XML SHA256 737a6219e93ae1cab48424300531515fea5f3b0adefdb31f5023e7cd3b385302; explicit lack of TPMI subtype data; no supplementary-material node in this XML. eMethods access still unverified. https://pmc.ncbi.nlm.nih.gov/articles/PMC13239414/
- Legacy API audit found one literal match in the scanner itself and no direct v1 call in its bounded tracked-code scan. Dynamically assembled URLs, external dependencies and other branches not covered. No global code migration warranted from this finding.

## Candidate priorities and constraints
1. GENESIS-K: 5,328 enrolled patients across 30 centers Oct 2023–Nov 2025 (ESOC2026A1666, https://academic.oup.com/esj/article/11/Supplement_1/i452/8669192). Not 5,328 LAS cases. Public search results also surface KBDS/KEA records KSE102975 and KSE102978 (https://kbds.re.kr/KSE102975 and https://kbds.re.kr/KSE102978), describing a prospective stroke pharmacogenomics cohort. Exact record-to-cohort identity, access status, linked TOAST data, stroke-free controls and sample overlap not confirmed; direct page retrieval failed in this review. These are metadata leads only; do not retrieve patient genotypes.
2. Taipei Veterans General Hospital 716-case cohort: secondary candidate for authorized TOAST genotype-phenotype analyses, not a proven independent LAS case-control GWAS. Same requirements apply. No request sent.
3. European comparator: still unqualified. New release date or a generic ICD stroke phenotype does not establish independent matched LAS replication. GIGASTROKE itself explains the subtype limitation of large biobank replication (https://www.nature.com/articles/s41586-022-05165-3).

## Decision
Prioritize an approved aggregate-only request to GENESIS-K using all 24 original markers, then Taipei if appropriate. Require TOAST ascertainment, eligible controls, model scale, allele identity, EAF, per-marker N, QC, recruitment and overlap, uncertainty estimates and applicable governance. If only case-case data exist, a distinct protocol and endpoint are needed. Do not turn an East Asian association into replication of the EAS–EUR difference. Stop repeated public catalog searches for these same publications unless their release status changes. No novelty upgrade or scientific completion upgrade.
