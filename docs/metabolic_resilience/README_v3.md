# METABOLIC RESILIENCE CODE — v3

## Operating model
- Drive: stage files + Colab runner backup.
- GitHub: canonical reviewed code and CI gate.
- /srv/is-analysis: data, execution, and result store.
- Colab: pipeline view and execution commands only.

## Prebuilt stages
- 05–33: discovery freeze, full pGWAS, coordinate/cis audit, marginal pQTL, 1000G EUR, clumping, frozen IVs.
- 40–47: outcome registry, harmonization, MR, FDR/overlap, Steiger.
- 50–55: coloc, multisignal decision, coding/assay artifact audit.
- 60–65: CKB/EAS registry, EUR candidate gate, cross-ancestry summary, EAS LD plan.
- 70–77: KoGES GRS, feasibility, variant coverage, analysis contract, longitudinal models, PA interaction, missingness, model diagnostics.
- 80–81: functional annotation and descriptive evidence matrix.
- 90–95: status, master gate, environment snapshot, Drive-Git drift, run plan, checksums.

## Locked rules
1. T2D is a separate disease-validation endpoint, not a fifth domain.
2. Confirmatory cis = GRCh37 protein-coding gene ±1 Mb; UKB-PPP coordinate from ID.
3. Primary pQTL: P<5e-8, F>10, INFO>=0.8; stringent sensitivity P<1.7e-11.
4. LD clumping r2<0.01; EUR uses 1000G Phase3 EUR503 v5b.
5. EAS primary analyses must not reuse EUR LD.
6. Weighted median k>=3; MR-Egger k>=10.
7. T2D coloc HOLD until case fraction is verified.
8. Multi-signal loci are routed to conditional/SuSiE sensitivity.
9. Protein-altering/coding/splice pQTLs trigger exclusion sensitivity.
10. KoGES primary = Cox incident MetS; secondary = repeated MBI LMM; PA interaction exploratory.
11. Functional annotation is supportive and cannot override MR/coloc evidence.
12. Final causal claims remain HOLD until master gate passes.
