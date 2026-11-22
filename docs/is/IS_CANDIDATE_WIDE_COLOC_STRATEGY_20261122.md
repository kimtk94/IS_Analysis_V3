# IS candidate-wide colocalization plan — 2026-11-22

**Scope:** MASTER DEGREE ischemic-stroke study. Exclude Sales-related research.  
**Decision:** Do NOT process FGF5, ALDH2, SH3PXD2A, COL4A2 as the only candidate pool. Run one *shared* candidate gate across the complete discovery universe first, then assay-level parallel statistical reruns only for scientifically qualified inputs.

## Verified denominators: gene != assay group != locus != independent replication

| Tier | Observed count | Actual meaning | Gate |
|---|---:|---|---|
| Broad discovery universe | 2,225 unique positional gene IDs | Discovery candidate reference, NOT 2,225 causal or coloc-eligible genes | Include as denominator and distinguish NOT_TESTED |
| Public QTL nominal source | 149 gene IDs | Nominal-QTL available, still without independently approved allele/LD | Prioritize original-source retrieval |
| Regional pilot not yet executed | 629 gene IDs | Regional analysis available only as design/pilot | Do not assign false negative |
| No rows in targeted slices | 91 gene IDs | Missing targeted source coverage, NOT null biological effect | Additional source query |
| Not in EAS 30 groups | 1,356 gene IDs | Outside particular EAS region tests | Preserve in discovery universe |
| P3 stage archived readiness | 393 gene×tissue×QTL/source groups | Repeated tests; grouped rows are not independent genes | Shared source metadata / error gate |
| P3A source-flagged follow-up | 36 *groups* | Possible ABF follow-up labels, not 36 independent genes | All included in batch screen |
| Source-stat-positive among 36 | 25 *groups* | Some source-flagged numeric QTL/GWAS rows, not validated complete cis summary | No fresh coloc allowed |
| Source-stat-unavailable among 36 | 11 *groups* | Inadequate archive statistics | Separate missing-data tier |
| Independently eligible new P3A coloc | **0** | Source GWAS SE/N, case fractions, authoritative alleles and complete cis unresolved | **HOLD** |

Primary pinned provenance: existing [Issue #193 aggregate evidence](https://github.com/kimtk94/IS_Analysis_V3/issues/193#issuecomment-6247948556); previous [P3A audit](IS_P3_ABF_VS_SUSIE_GATES_20261017.md), [GWAS identity protocol](IS_P3A_GWAS_SOURCE_MATCH_PROTOCOL_20261018.md). Baseline stored at [IS_CANDIDATE_WIDE_COLOC_BASELINE_20261122.json](evidence/IS_CANDIDATE_WIDE_COLOC_BASELINE_20261122.json).

## Execution layers

1. **Full universe:** process all 2,225 candidate IDs as a data-readiness atlas; do not require each candidate to be in a currently assayed cis region, and never label untested as negative.
2. **Assay-level gate:** for each existing 393 group, use natural composite key `region_id × QTL_dataset × gene_id` (and GWAS source once verified). Preserve original group IDs rather than inferred gene symbols.
3. **Parallel P3A batch:** source-stat-positive 25/36 groups are **first-class records**, along with the 11 blocked groups and the 357 not source-flagged. For actual `coloc.abf` use approved complete cis input, signed GWAS/QTL beta and SE or scientifically approved alternative, N and case fractions, valid AF/allele/build, trait model and prior sensitivity. A real computed PP.H0–PP.H4 is required for executed status.
4. **Parallel P3B batch:** run `susie_rss`/`coloc.susie` on *each qualified independent locus×tissue×QTL source* only when a signed ancestry/matched SNP-ordered LD matrix, donor N, credible-set checks and multi-signal diagnostics are approved; reference EUR/EAS proxy LD is sensitivity-only.
5. **Meta-evidence reduction:** aggregate by gene for presentation, but preserve the actual assay rows, tissue/dataset provenance, overlap SNP count, PP.H3, PP.H4, conditional ratio H4/(H3+H4), method, versions, N/p12/priors sensitivity, LD QC and explicit N/A values. Do not choose the maximum PP.H4 across dozens of assays without addressing selection/multiplicity, correlated tissues and discovery source overlap.
6. **External replication:** require independent stroke GWAS ascertainment; GIGASTROKE EAS discovery contains BBJ, so BBJ vs GIGASTROKE EAS is only descriptive comparison unless a BBJ-excluded set is sourced and validated.

## Four historical anchors are comparators, not the full universe

- **FGF5** archived inherited-source H4~0.779; separately reprocessed EBI×BBJ exploratory ABF H4~0.913 at assumed QTL N=175/p12=1e-5 (N/p12 sensitivity 0.370–0.991). Not independent replication.
- **ALDH2 / rs671** remains a coding/alcohol/BP/LD-pleiotropy *variant-centric* path, not automatically a coloc-or-MR causal protein/gene.
- **SH3PXD2A** and **COL4A2** are supporting cell/vascular hypotheses, not declared non-significant based on exploratory inherited H4.
- Existing [four-locus source-read audit](IS_FOUR_LOCUS_EVIDENCE_AUDIT_20261122.md) proves input table QC (7,969 rows), not recomputation of 646 historical coloc tests.

## Actually installed automation in this change

`scripts/is/candidate_wide_coloc_source_triage.py` + tests + `is-candidate-wide-coloc-triage.yml`:
- Validate *real existing* frozen aggregate 2,225 / 393 / 36 / 25 / 11 and the separate four-locus QC snapshot.
- Fail closed on fabricated approved N, false eligibility, group duplicate, missing source metadata or H4/replication status promotions.
- Optionally read the **actual original** `is-archived-fasta-abf-audit-v1` JSON including its 36 individual `abf_flagged_groups` and produce deterministically keyed per-group TSV; this is future-proof but requires the *original group-level audit JSON*, not aggregate issue text alone.
- The current Actions job **does not have original server-local 36-group JSON or private Drive genomic data**; it emits an **aggregate-only** source-verification output and does not invent 36 gene IDs. Synthetic multi-group tests are *software tests*, not biological data scans.
- The source audit JSON can later be committed in reviewed **metadata-only** form at `docs/is/evidence/IS_P3A_ARCHIVED_FULL_GROUP_AUDIT.json` (no private paths, individual SNP effects or donor data) to enable batch gene/dataset/region TSV generation automatically. No secret or remote server dependency is added.
- This change is deliberately limited to a **candidate-wide screening + scientific gate**, not an automatic run of coloc on unavailable/inadequate GWAS/QTL inputs.

**Next scientific bottleneck:** recover and independently verify original complete summary GWAS and molecular QTL effect/SE/case fraction/allele sources. Then run the approved groups together via bounded GitHub Actions matrix and combine a method- and cohort-aware evidence table, rather than validating FGF5 in isolation.
