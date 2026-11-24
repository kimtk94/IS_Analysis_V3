# Manuscript decision review — 2026-11-24
Status: **PROPOSED / human review required**. A focused literature audit, not a systematic review, registration, course submission, or completed scientific analysis. No research track is declared submission-ready.

## Scope and traceability
Nine tracks per [PR 333](https://github.com/kimtk94/IS_Analysis_V3/pull/333). Ocular/glaucoma is one track; Radiomics and harness are excluded. This documentation-only branch is based on main 2e5c4d36a91acd14ad8f31fade62d48ef5dd670a. Recent research evidence remains on its original branches: [IS #330](https://github.com/kimtk94/IS_Analysis_V3/pull/330) (merged), [#333](https://github.com/kimtk94/IS_Analysis_V3/pull/333), [#334](https://github.com/kimtk94/IS_Analysis_V3/pull/334), [#335](https://github.com/kimtk94/IS_Analysis_V3/pull/335) (drafts). No merge, code change, rerun, data download, private-data upload or automatic scientific promotion.

- [Novelty matrix](novelty_matrix.tsv): nearest competing work, what is already known, remaining candidate gap, and decision.
- [Decision protocols](decision_protocols.md): one decision-oriented question per track; endpoints, prerequisites, analysis, and stop rules.
- [Evidence register](evidence_register.tsv): claim-level evidence and next blocker. Update only against direct source/analysis readback.

## Portfolio disposition
| Track | Priority | Disposition | Next decisive deliverable |
|---|---|---|---|
| Ocular | 1 near-term | CONDITIONAL_GO | Unselected donor-level validation cohort and complete original paired score analysis |
| CKD | 2 genetic | CONDITIONAL_GO | Match independent EAS pQTL assay IDs, instruments, outcome phenotype and overlap; determine whether comparison is estimable |
| IS | Thesis core | CONDITIONAL_GO | One locus/one cell-context hypothesis with complete GWAS-QTL regional stats and credible LD; explicit cross-ancestry estimability |
| Muscle | Next exploratory | REDESIGN | Participant-level baseline-only features and continuous hypertrophy outcome availability; external cohort candidate |
| CKM–Stroke | IS support | LINK_TO_IS | Formal common-instrument subtype contrast only if distinct question and outcome sources available |
| HCC | Conditional | HOLD | Independent, paired RNA/protein/clinical cohort and locked clinical endpoint |
| Oral | Conditional | REDESIGN | Exact patient-primary/LN linkage and independent matched cohort; oral/HNSCC anatomical definition |
| Skin aging | Feasibility | HOLD | Full phenotype-specific GWAS and tissue-QTL coverage, independent endpoint feasibility |
| Metabolic | Controlled data | HOLD | Approved longitudinal genotype-phenotype intersection and recovery/relapse events |

## Common scientific gates
1. Novelty: compare closest original studies (including recent preprints), populations, accessions and claims; distinguish new question from new software.
2. Protocol: write estimand, unit of inference, endpoint, population, inclusion rules, smallest meaningful effect, uncertainty, multiplicity family and missing-data strategy. Freeze and version before **new** holdout outcomes are inspected; previously seen results are discovery. Human approval required; this is a draft, not preregistration.
3. Data: cohort/donor identity, overlap, legal access, version and checksums, attrition, covariates, scale/build/alleles, realistic power. Missing/unidentifiable is HOLD, not negative evidence.
4. Discovery: execute only fit-for-question analysis; report denominators, estimates, confidence intervals, diagnostics, nulls and failures.
5. Robustness: prespecified sensitivity models and bias checks; do not choose thresholds or subgroups after seeing P values.
6. Replication: independent data or explicitly label reused data as internal/conditional reproducibility; no pooling pseudoreplicates.
7. Interpretation: molecular localization is not causal mediation, association is not clinical prediction, and a difference in significance is not a significant difference.
8. Manuscript gate: trace every claim and figure to versioned outputs; human review of novelty, uncertainty, reporting and access/licensing. Publication may be justified by informative null/limits, not just significant results.

## Execution policy
GitHub Actions only for permitted public-data computations; new narrowly scoped manual-dispatch job after scientific gates justify it, not a repeated readiness-only matrix or harness. Existing CI unchanged. Controlled KoGES, dbGaP and private CKD inputs remain in approved environments; no credential, participant-level record or raw archive in Git. Do not assume public source download succeeds from catalogue listing. Preserve source URL/ID, accession, cohort, reference build, source hash, actual output artifact and its expiry; code and evidence lineage only under applicable data rules. Green CI != inferential PASS. Do not change open analysis PRs or merge without review.

Literature audit performed 2026-11-24 using targeted searches for each track, accessions, gene names, and competing methods; nearest papers below are verified entry points, not exhaustive proof of novelty. Recheck full texts, supplemental methods and citation networks prior to final protocol lock. In particular, no approved course submission text is produced here.
