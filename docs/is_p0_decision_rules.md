# IS P0 — Candidate evidence policy (proposal v0.1, 2026-10-10)

This document records **decision rules before the next screening**, not a retrospective
declaration that genes or loci are causal.

## Scope and immutable evidence boundaries

- Discovery: BBJ/GIGASTROKE EAS/Japanese, and EUR AIS **distance-merged provisional regions**.
  Do **not** label these 80 components independent GWAS loci without LD-clump/fine-map verification.
- Candidate generation: all 2,425 region–gene positional associations from 2,225
  unique gene IDs. Positional overlap is *not* proof of a target.
- Molecular evidence: the 646 successful GTEx v8 ABF tests were limited to four
  original BBJ loci. Do **not** use "646 tests" as the denominator for expanded
  genome-wide/EUR gene discovery.
- Orthogonal evidence: GTEx coloc/ABF and coloc/SuSiE share GWAS/QTL inputs;
  these are sensitivity/alternative models, **not independent replication**.
- Independent evidence includes genuinely distinct cohort QTL, disease-cell QTL,
  validated coding mechanism, and credible regulatory variant-to-gene experiments.
  scRNA expression localization alone is descriptive, not molecular causality.

## Precommitted assessment stages

| Gate | Criterion | Allowed interpretation |
|---|---|---|
| G0 source QC | Build, sample ancestry, phenotype, variant positions/alleles, variant overlap, reference LD, QTL provenance and per-gene/tissue test inventory verified | Eligible for molecular testing |
| G1 association | GWAS significance and region independence based on appropriate clumping/fine-mapping; preserve full GWAS/SNP and missingness ledgers | Locus associated with phenotype |
| G2 molecular | Report H0–H4, H4/H3, priors and p12 sensitivity, shared SNP N, tested tissues/genes, effect direction, QTL power; characterize multi-signal cases | Shared-variant hypothesis, strength graded not binary |
| G3 independent | Independent cohort/tissue QTL, appropriate cis-pQTL MR with pleiotropy checks, or functional enhancer/allele-specific evidence | Supported gene/mechanism hypothesis |
| G4 cell | Disease-relevant donor-level cell specific expression/regulation or QTL; compare relevant cell classes and absent cells | Cell-type hypothesis, not causal proof alone |
| G5 interpretation | Association, shared variant, mediated exposure, cell context and functional causality separated, with alternatives and negative checks recorded | Thesis-calibrated claims |

## Evidence status definitions

- SUPPORTED: multiple **independent** relevant evidence types converge and
  critical QC/ancestry/pleiotropy objections have been addressed. Claim must
  state precisely what is supported (locus, gene, mediation, or cell context).
- PRIORITIZED_UNVERIFIED: important GWAS/coding/molecular lead, but lacks
  independent evidence or required QC. All four historical "core" loci currently
  remain hypotheses unless G2/G3/G4 support is demonstrated.
- EXPLORATORY: plausible based on region, ABF or biology but insufficient
  confirmation. A high selected H4 alone does not promote a gene.
- INDETERMINATE_POWER_OR_COVERAGE: limited QTL sample size, failed SuSiE
  signal detection, missing matrix feature, missing variants, or mismatched
  ancestry reference. Do **not** describe as negative biological evidence.
- MECHANISM_NOT_SUPPORTED: an explicitly testable specified mechanism fails
  a predeclared, sufficiently powered and QC-passing falsification test.
  This does **not** imply a locus or the entire gene is biologically irrelevant.

## Explicit fail/hold logic

1. Do **not** convert H4 below 0.8 into a global rejection. H4 is a posterior
   sensitive to priors, signal strength, and region assumptions. Use H3 relative
   to H4, H4/H3, and p12 robustness; report the entire screen and a predeclared
   robustness grid. High H4 from post-hoc tissue selection remains exploratory.
2. If the QTL feature is **missing** or quality/LD prerequisites fail, set
   INDETERMINATE_POWER_OR_COVERAGE, not MECHANISM_NOT_SUPPORTED.
3. SuSiE no credible QTL signal is *not* a statistically established H3
   result and should not be interpreted as a falsification of shared causality.
4. To mark a **specific** macrophage mechanism unsupported, require an
   independently powered disease-relevant cell/QTL contrast with adequate
   markers, donor representation, region signal mapping and predefined effect
   precision threshold; descriptive healthy scRNA localization alone cannot.
5. FGF5 R3_4 adult temporal-lobe counts UNRESOLVED:
   first inspect original Seurat RNA features/symbol–Ensembl mapping and
   preprocessing filters; if absent, test other reference datasets. Until then
   MISSING_FEATURE, not zero expression or FGF5-negative.
6. For ALDH2 rs671, keep coding/mediated metabolic mechanisms separate.
   Do not interpret a single pleiotropic instrument as conclusive cis-MR.
7. No automatic causal-gene promotion based on a score without independent
   support and manual review. Record missingness and alternative explanations.

## Selection transparency and the expansion decision

- Preserve the original legacy core (FGF5, ALDH2, SH3PXD2A, COL4A2)
  and secondary (CALHM2, NEURL1, C4orf22, INA, COL4A1)
  as **historical hypothesis labels**, not fixed proof-based categories.
- CALHM2 ABF H4≈0.798 is numerically higher than FGF5 ABF H4≈0.779.
  FGF5 legacy priority reflects earlier locus/mechanism focus and is **not**
  justified by ABF rank alone; explicit signal, coding, pQTL, cell and
  cross-dataset comparisons are required before a thesis claim.
- Initial 646 ABF test rows cover selected 4-locus GTEx pairs; expanded
  discovery candidate universe is broader (EAS/Japanese + EUR) and must
  retain ancestry/phenotype distinctions.
- Continue **candidate inventory expansion** to close documented data gaps,
  but defer large-scale molecular tests on every expanded gene until
  G0 completion, independent signal QC, cost/power review and prioritization.
- Report all screened gene–tissue pairs, the QC denominator, filtered-out
  pairs, posterior results, missingness, and choice of selected examples.
  Coloc H4 posterior is **not** itself a multiplicity-adjusted discovery rate.

## Reproduction and non-destructive execution

Run read-only against data roots:

    python3 scripts/is/audit_is_p0_evidence.py \
      --root /srv/is-analysis \
      --out /srv/is-analysis/results/is/audits/p0_evidence_20261010

Output must be a new empty folder. Never overwrite input tables, phase
folders, or prior evidence without versioning. Server commands must not rely
on set -e. Commit source and tests only; never commit raw data or audit
outputs. This document is a decision policy, not a record of an executed
coloc rerun or of newly established causality.
