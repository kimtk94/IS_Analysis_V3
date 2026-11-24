# Publication-directed portfolio (2026-11-24)

The user excluded Radiomics. The active portfolio has nine research tracks.
Harness recovery, scheduling, remote deployment and automated promotion are out of scope.
Use bounded, separate GitHub-hosted Actions for public input computations, not repeated readiness-only jobs. Do not send controlled individual data, credentials, or prior private Drive inputs to GitHub. Do not treat a green job as statistical evidence.

| Track | Next decision-oriented work | Advance criterion / limitation |
|---|---|---|
| Ocular / glaucoma | Compare full donor-paired and pool-paired effects against published independent TM bulk at the already specified effect-size strata. Use overlap denominators and missing FDR explicitly. | Manuscript on reproducibility and limitations; not evidence of clinical IOP effects or separation of treatment-correlated preparation. A genuinely unselected external test remains necessary for confirmatory claims. |
| CKD | Retrospective EUR-only selection, then all-family EAS anchor Bonferroni; private aggregate input remains outside Git. Follow with genuinely unused pQTL/outcome, model/assay equivalence and overlap assessment. | Separate significant EAS anchors from shared-signal evidence; no selection-reused direction as independent replication. |
| IS | Select a precise tissue- or subtype-specific question beyond published FGF5/ALDH2 pQTL MR. Run new coloc only with valid full regional SNP effects and stated QTL scale/N/prior. | Existing 15-assay sensitivity is exploratory, not formal validation of the 36-assay P3A. |
| HCC | Establish externally verified matched RNA/protein clinical cohort and a nonduplicative endpoint before refitting. | Same-cohort 48-patient reproduction and sample-label audit are not independent clinical validation. |
| CKM–stroke | Predefine LDL-C and SBP contrasts across LAA, CE, SVD; recover effect estimates/SE and sampling covariance, then test direct differences with a defined six-contrast family. | No comparison of significant versus nonsignificant as a heterogeneity test; no LDSC/MR until suitable sources and covariance are available. |
| MUSCLE | Run the actual 20 RT lead rsID-to-GTEx-v10-variant mapping, followed by skeletal-muscle significant eQTL/sQTL queries using the existing Stage 2A implementation. | A significant QTL at a suggestive RT lead is regulatory annotation, not GWAS replication or coloc. API errors/unresolved variants remain separate. |
| Skin aging | Try actual bounded public transfers for GCST90094903 / GCST90094904 (objective Korean wrinkles). Verify complete phenotype and allele schema before MR. | Do not substitute UKB perceived age or treat catalogue presence/headers as full-file validation. |
| Oral cancer | Read the actual public GSE103322 processed TPM matrix and inspect explicitly identified patient/site metadata; use official ExperimentHub EH5419 if GEO remains blocked. | Patient is the unit of inference; do not treat TPM as counts or infer clinical labels from names. Need a novel question and matched patient counts before inferential models. |
| Metabolic Resilience | Establish approved actual KoGES genotype/phenotype intersections and incident MetS events; keep predefined Cox primary and repeated MBI secondary. | No private cohort upload to Actions; no invented N, endpoint or positive result. |

Initial public jobs: Ocular three-model comparison, MUSCLE QTL search, Skin source access and Oral processed matrix structure. The CKD calculation CLI can run against the preserved input locally; no results are committed. Existing data sources, outputs and dates remain distinct. The new workflow has no schedules, server secrets, write token or harness dependency; it only runs when specifically dispatched or its own analysis files are changed in a PR. Results are generated as artifacts, not committed. Scientific decisions remain human-reviewed.

The existing main smoke suite requires the four immutable historical aggregate input documents from 0bf646fd6f29322d17da878253f37def24c01c74; the CI prepares only those paths. The new tests exercise no network or real data.
