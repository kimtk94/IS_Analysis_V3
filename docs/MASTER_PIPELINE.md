# MASTER OMICS CAUSAL PIPELINE v1

This document is the shared backbone for CKD and ischemic stroke (IS) analyses.

## Design principle
One engine, disease-specific configuration. Each stage has a fixed contract:
Input -> Method -> QC -> Output -> Figure -> Evidence.

## Stages

| Stage | Name | Core purpose | Minimum output |
|---|---|---|---|
| 0 | dataset_audit | audit source metadata/build/ancestry/schema | dataset_registry.tsv |
| 1 | exposure | build cis-pQTL/eQTL exposure set | exposure_instruments.tsv |
| 2 | instrument_qc | LD/QC/F-stat/allele harmonization | instruments_qc.tsv |
| 3 | primary_mr | causal screening | mr_primary.tsv |
| 4 | mr_robustness | WM/Egger/Q/Steiger/LOO | mr_sensitivity.tsv |
| 5 | coloc | shared causal signal test | coloc_summary.tsv |
| 6 | finemap | SuSiE/credible-set signal resolution | finemap_summary.tsv |
| 7 | cross_ancestry | EUR/EAS replication | ancestry_replication.tsv |
| 8 | cross_platform | Olink/SomaScan replication | platform_replication.tsv |
| 9 | transcriptomics | eQTL/SMR/HEIDI/RNA colocalization | transcript_evidence.tsv |
| 10 | phenotype_expansion | related outcome screen | phenotype_matrix.tsv |
| 11 | risk_factor | mediator/risk-factor MR | risk_factor_mr.tsv |
| 12 | disease_subtype | disease-specific subtype analysis | subtype_results.tsv |
| 13 | individual_validation | KoGES/individual-level validation | clinical_validation.tsv |
| 14 | tissue | tissue expression/QTL evidence | tissue_evidence.tsv |
| 15 | single_cell | cell-type localization | celltype_evidence.tsv |
| 16 | spatial | anatomical localization | spatial_evidence.tsv |
| 17 | phewas | benefit/risk/pleiotropy scan | phewas_summary.tsv |
| 18 | druggability | drug-gene/target direction | druggability.tsv |
| 19 | evidence_integration | final evidence tiering | evidence_matrix.tsv |

## Mandatory QC conventions
- Preserve gene/protein, chromosome, position, effect allele, other allele, beta, SE, EAF, P and N whenever available.
- Genome build and ancestry must be explicit.
- F-statistic > 10 is the default instrument-strength rule unless disease config overrides it.
- Colocalization reports PP.H0-H4; PP.H4 >= 0.80 is the default strong-evidence threshold.
- Fine-mapping must report credible-set count, size, purity and convergence.
- Cross-ancestry replication must never use mismatched LD reference panels.
- All multiple-testing rules must be written to result metadata.

## Disease-specific emphasis

### CKD
Primary chain:
plasma pQTL -> MR -> coloc -> SuSiE -> EAS replication -> transcript/kidney evidence -> KoGES longitudinal validation -> PheWAS/druggability.

Primary phenotypes:
CKD, eGFRcrea, eGFRcys, BUN, UACR, rapid decline, annualized eGFR slope, and CKD subtypes.

### Ischemic stroke
Primary chain:
plasma pQTL -> MR -> coloc -> stroke subtype -> risk-factor MR -> EAS replication -> brain/vascular cell evidence -> PheWAS/druggability.

Primary phenotypes:
any IS, cardioembolic stroke (CES), large-artery stroke (LAS), small-vessel stroke (SVS).

Risk-factor panel:
SBP, DBP, LDL-C, BMI, T2D, atrial fibrillation and smoking.

## Evidence tier
- Tier 1: MR + strong coloc + replication + orthogonal biological/clinical evidence.
- Tier 2: MR + coloc with partial replication or biological support.
- Tier 3: MR signal without sufficient locus-level or replication support.

## Figure contract
1. Study design
2. Proteome-wide discovery
3. Colocalization/fine-mapping
4. Cross-ancestry/phenotype replication
5. Individual-level validation
6. Tissue/single-cell/spatial evidence
7. PheWAS/druggability/final evidence matrix

## Implementation rule
New disease analyses should add or modify config only where possible. Shared statistical logic belongs in reusable scripts; disease-specific branching should be driven by config, not copied code.
