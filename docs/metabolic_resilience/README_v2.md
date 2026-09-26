# METABOLIC RESILIENCE CODE — v2

## Operating rule

- **Google Drive**: one script per analysis stage.
- **Colab**: pipeline structure + execution commands only.
- **Server**: `/srv/is-analysis` remains the analysis runtime and result store.
- Never place the whole analysis implementation inside the Colab notebook.

## Locked scientific rules

- Stage 2 shortlist: AOC1, OGN, TNFRSF6B, TFPI, SULT1A1, CTRL, IDUA, COMT.
- T2D is disease-endpoint validation, not a fifth metabolic domain.
- Full UKB-PPP pGWAS GRCh37 coordinate = position embedded in `ID`, not `GENPOS`.
- Confirmatory cis = Ensembl GRCh37 target protein-coding gene interval ±1 Mb.
- Exposure effect allele = `ALLELE1`; exposure effect = marginal `BETA`.
- Primary IV eligibility = `P < 5e-8`, `F > 10`, `INFO >= 0.8`.
- Stringent sensitivity = `P < 1.7e-11`, `F > 10`, `INFO >= 0.8`.
- 1000G reference = Phase 3 EUR503 GRCh37 v5b.
- Keep biallelic SNPs + indels.
- LD clumping = `r2 < 0.01`; default clump radius in Stage 3-C1 = 10,000 kb.
- MR-Egger only when `k >= 10`.
- Weighted median when `k >= 3`.
- Single IV -> Wald ratio.
- Coloc default priors: p1=1e-4, p2=1e-4, p12=1e-5; p12 sensitivity 1e-6 and 1e-4.
- T2D coloc is HOLD until the no-UKB case fraction is verified; never infer it.

## New code files

- `32_stage3c1_ld_clump.py`
- `33_stage3c2_freeze_instruments.py`
- `40_stage3d0_discover_outcomes.py`
- `41_stage3d1_extract_outcome_instruments.py`
- `42_stage3d2_harmonize.py`
- `43_stage3d3_mr.R`
- `44_stage3d4_run_confirmatory_mr.sh`
- `45_stage3d5_summarize_mr.py`
- `50_stage3e0_coloc_preflight.sh`
- `51_stage3e1_prepare_coloc_regions.py`
- `52_stage3e2_coloc_abf.R`
- `53_stage3e3_summarize_coloc.py`
- `60_stage4a_ckb_coverage_audit.py`
- `61_stage4b_eas_replication_scaffold.py`
- `70_stage5a_build_koges_grs_weights.py`
- `90_pipeline_status.py`
- `00_METABOLIC_RESILIENCE_RUNNER_v2.ipynb`

## Execution gate

Do not automatically advance past a gate when the prior stage returns `REVIEW`, `HOLD`, `WARN_LOW_NSNP`, or `FAIL`. Review the audit file first.


## Stage 4-5 continuation added in Git

- `62_stage4c_select_eur_candidates.py`
  - EUR-only gate before EAS replication.
  - MR FDR <0.05 + favorable direction + coloc PP.H4 >=0.80.
  - Multi-domain priority requires support in >=2 metabolic domains.
- `63_stage4d_eas_registry_validate.py`
  - Explicit CKB / KoGES / BBJ / TWB schema and sample-overlap registry.
  - Resource-specific schema is never guessed.
- `64_stage4e_cross_ancestry_summarize.py`
  - Separates concordant-significant, concordant-nonsignificant,
    opposite-direction, and not-testable EAS results.
- `71_stage5b_koges_feasibility_gate.py`
  - Implements the prespecified KoGES GO / MODIFY / NO-GO logic.
- `72_stage5c_validate_grs_variant_coverage.py`
  - Audits KoGES genotype coverage and dosage-allele orientation for external
    cis-pQTL GRS weights before participant-level scoring.
- `73_stage5d_write_analysis_contract.py`
  - Locks incident MetS as primary and repeated MBI as secondary before fitting.
