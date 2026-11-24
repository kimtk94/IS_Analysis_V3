# CKD MASTER end-to-end bridge plan

This runbook maps the existing CKD research outputs into the shared MASTER Stage 0-19 model without rewriting or discarding the validated CKD-specific pipeline.

## Principle

The existing CKD workflow remains the source of truth for completed scientific results. The MASTER layer is an orchestration and evidence-integration layer.

Therefore, a stage is classified as:

- **EXISTING**: a current CKD artifact directly represents the stage.
- **PARTIAL**: evidence exists, but the generic MASTER output still needs to be materialized.
- **PLANNED**: required input/result does not yet exist as a canonical CKD result.

An existing file does not automatically mean the conceptual MASTER stage is complete.

## Current bridge

| MASTER | Current CKD source | Bridge status |
|---|---|---|
| 0 Dataset audit | analysis-ready SHA256/public manifests | PARTIAL |
| 1 Exposure | stage1/instrument_qc.tsv.gz | EXISTING |
| 2 Instrument QC | stage1/instrument_qc.tsv.gz | EXISTING |
| 3 Primary MR | stage1/stage1_protein_summary.tsv | EXISTING |
| 4 MR robustness | new MASTER output required | PLANNED |
| 5 Coloc | stage2b_coloc/.../STAGE2B_COLOC_DEFAULT.tsv | EXISTING |
| 6 SuSiE | stage2c_susie/.../STAGE2C_SUSIE_DEFAULT.tsv | EXISTING |
| 7 Cross-ancestry | stage1 protein summary contains EUR/EAS support | PARTIAL |
| 8 Cross-platform | independent platform input required | PLANNED |
| 9 Transcriptomics | Stage3A kidney QTL evidence | PARTIAL |
| 10 Phenotype expansion | Stage1 renal phenotype screen | PARTIAL |
| 11 Risk-factor | new analysis | PLANNED |
| 12 CKD subtype | new analysis | PLANNED |
| 13 KoGES | public prototype / anchor preparation | PARTIAL |
| 14 Tissue | Stage3A kidney evidence | EXISTING |
| 15 Single-cell | Stage3B integrated evidence | EXISTING |
| 16 Spatial | new evidence layer | PLANNED |
| 17 PheWAS | new analysis | PLANNED |
| 18 Druggability | curated target/drug input required | PLANNED |
| 19 Final integration | generated after bridge | PLANNED |

## First server command

Run this after pulling the branch/PR:

```bash
cd /srv/is-analysis/IS_Analysis_V3

CKD_MASTER_ROOT=/srv/is-analysis/results/ckd/master \
  bash server/ckd_master_end_to_end_audit.sh
```

Expected outputs:

```text
/srv/is-analysis/results/ckd/master/
├── CKD_MASTER_STAGE_READINESS.tsv
└── CKD_MASTER_STAGE_READINESS.json
```

The audit is read-only and checks the actual server filesystem.

## Production sequence after the audit

1. Bridge Stage 0/1/2 metadata into canonical MASTER inputs without recalculating established CKD results.
2. Materialize Stage 4 robustness from the existing harmonized Stage1 multi-instrument tables where statistically appropriate.
3. Convert Stage1 EUR/EAS support to the generic Stage7 replication table.
4. Convert Stage3A kidney QTL evidence to Stage9 transcript evidence while keeping formal SMR/HEIDI separate from mere eQTL hit evidence.
5. Build the Stage10 phenotype matrix from existing CKD/eGFR/BUN/eGFRcys/UACR MR outputs.
6. Run controlled-access KoGES Stage13 only on the local server; export aggregate model summaries only.
7. Normalize Stage3A/3B into Stage14/15 localization; add Stage16 only when spatial evidence is truly available.
8. Add PheWAS and druggability inputs.
9. Run Stage19 final integration only when each upstream table has provenance and QC.

## Non-negotiable boundaries

- Public KoGES training data are workflow QA, not manuscript-level genotype inference.
- Stage3A kidney eQTL hit counts are not a substitute for formal SMR/HEIDI.
- Cross-platform evidence is not assumed if only Olink data exist.
- Spatial evidence is not inferred from cell-type localization.
- A missing planned stage remains missing; the bridge must not award evidence points for absent analyses.


## Materialize reusable existing evidence

After the readiness audit confirms the expected files, materialize only the reusable existing evidence:

```bash
cd /srv/is-analysis/IS_Analysis_V3

CKD_MASTER_ROOT=/srv/is-analysis/results/ckd/master \
  bash server/ckd_master_materialize_existing.sh
```

This creates, where prerequisites are available:

```text
stage07_cross_ancestry.tsv
stage10_phenotype_manifest.tsv
stage10_phenotype_long.tsv
stage10_phenotype_matrix.tsv
stage14_15_localization.tsv
stage19_candidates.tsv
stage19_evidence_manifest.tsv
```

Stage 4 robustness is deliberately off by default because the multiple UKB-PPP ST16 conditional signals can retain residual LD. For sensitivity-only materialization:

```bash
RUN_STAGE4_SENSITIVITY=1 \
CKD_MASTER_ROOT=/srv/is-analysis/results/ckd/master \
  bash server/ckd_master_materialize_existing.sh
```

A provisional Stage19 matrix can be produced for pipeline QA only:

```bash
RUN_PROVISIONAL_STAGE19=1 \
CKD_MASTER_ROOT=/srv/is-analysis/results/ckd/master \
  bash server/ckd_master_materialize_existing.sh
```

The provisional matrix is explicitly not the manuscript final tiering because Stage8/9/11/12/13/16/17/18 may still be absent.
