# MasterOmics canonicalization decision

Status: **ACTIVE**
Canonical execution engine: **MasterOmics (feat/masteromics-central-engine, PR #24)**
Reference implementation source: **MASTER pipeline v1 (feat/master-pipeline-v1, PR #28)**

## Decision

MasterOmics is the single execution backbone for CKD and ischemic stroke.

PR #28 is not merged as a second orchestration framework. Its scientifically useful
modules are migrated into MasterOmics one at a time behind the existing MasterOmics
contracts, checkpointing, provenance, resource registry, and CI.

This avoids three concurrent execution models:

1. legacy disease scripts
2. MASTER pipeline v1
3. MasterOmics

Legacy scripts remain frozen references until numerical parity or explicit replacement
is demonstrated. They are not deleted solely because a central adapter exists.

## Why MasterOmics owns execution

MasterOmics already provides:

- validated DAG execution and dependency contracts
- resumable checkpoints and input/code hashing
- atomic validated output promotion
- explicit unbound/data-not-validated states
- strict full-region vs MR-instrument separation
- signed ancestry/build-matched LD contracts
- shared MR, coloc, SuSiE, cohort and evidence interfaces
- pinned/resumable resource acquisition with provenance
- CKD frozen-result replay and raw-rebuild comparison
- GitHub synthetic/scientific contract tests

PR #28 has a broader translational stage vocabulary and several useful generic modules,
but running its orchestrator independently would recreate duplicate state, configs,
stage semantics, and result ownership.

## Stage mapping

| MASTER pipeline v1 concept | Canonical MasterOmics owner | Migration status |
|---|---|---|
| dataset_audit | acquisition + source_qc | covered |
| exposure | normalize + cis_regions | covered; source adapters remain explicit |
| instrument_qc | ld_qc + instruments + harmonize | covered |
| primary_mr | discovery_mr | covered |
| mr_robustness | sensitivity | covered |
| coloc | regional_coloc | covered |
| finemap | finemap | covered |
| cross_ancestry | replication | covered as contract; project adapter required |
| cross_platform | replication | **port candidate** |
| transcriptomics | annotation | **port candidate** |
| phenotype_expansion | evidence / project adapter | **port candidate** |
| risk_factor | evidence / project adapter | **port candidate; high priority for IS** |
| disease_subtype | replication/evidence project adapter | **port candidate; high priority for IS** |
| individual_validation | score + longitudinal + incident | covered for cohort contract; project adapter required |
| tissue | annotation | covered as contract |
| single_cell | annotation | covered as contract |
| spatial | annotation | **port candidate** |
| phewas | evidence | **port candidate** |
| druggability | evidence | **port candidate** |
| evidence_integration | evidence | covered as contract; richer scoring can be ported |
| manuscript/report | report | contract exists; report adapter pending |

## Migration rules

A PR #28 module is promoted into MasterOmics only when all of the following are true:

1. it consumes or produces an explicit MasterOmics artifact contract;
2. ancestry, build, allele direction, testing family and provenance are explicit;
3. it has synthetic/unit tests in the MasterOmics CI;
4. it does not duplicate an existing MasterOmics implementation;
5. production use remains blocked until source-specific data review is complete;
6. numerical differences from a frozen legacy result are reported, not silently accepted.

## Priority order

### P0 — IS production adapter

The current MasterOmics regression/rebuild path intentionally rejects ischemic-stroke
legacy schemas. The next centralization milestone is an explicit IS baseline adapter
that maps the existing BBJ/GIGASTROKE and functional-validation outputs into canonical
MasterOmics artifacts without changing the frozen IS results.

### P1 — IS disease-specific modules

Port in this order:

1. disease subtype
2. risk-factor MR
3. transcriptomics / cell-regulatory annotation
4. cross-platform replication
5. PheWAS / druggability

These should become project adapters or evidence modules, not a second runner.

### P1 — CKD completion

Keep the current CKD V1.1 frozen package as the numerical reference. Complete bounded
SuSiE production validation and controlled KoGES adapters when data are available.

### P2 — broader research programs

Metabolic Resilience, Muscle, and Skin Beauty should reuse the same contracts only after
CKD/IS centralization is stable. They should not force premature generalization of the
core disease engine.

## PR handling

- PR #24 remains the canonical central-engine PR.
- PR #28 remains a reference/migration source until its useful modules are ported.
- Do not merge both orchestration layers into main independently.
- Do not delete PR #28 code before a module-level migration inventory is complete.

## Source-of-truth boundary

Execution truth:

server data/results -> MasterOmics checkpoints/contracts

Presentation truth:

server checkpoints -> generated research manifest -> MasterOS/Dashboard

Static dashboard files may describe durable study design and figure scaffolding, but must
not own current status, live metrics, live findings, or next actions.
