# MasterOmics shared research structure

Current task: architecture/scaffold implementation. Data review is explicitly
DEFERRED. Initialization never claims scientific execution. CKD and ischemic stroke
use the same MasterOmics engine, checkpointing and artifact-contract machinery, but
now have project-specific scientifically valid stage catalogs. CKD remains pQTL/MR-first;
IS follows the frozen BBJ locus-first thesis design. Legacy source scripts and results
remain unchanged.

## Workspace

Each project has `project.json`, `config`, `raw`, `reference`, `interim`,
`standardized`, `instruments`, `regional`, `results`, `state`, `logs`, `reports`.
Config includes a source registry, analysis-policy placeholders, the stage DAG
and explicit adapter bindings. All paths are configurable; no credentials,
restricted cohorts or downloaded data belong in Git.

## CKD stage contracts and reuse map

| Stage | Artifact contract | Existing module / remaining boundary |
|---|---|---|
| acquisition | reviewed raw-source identity | science.acquire / local archive adapters |
| source_qc | integrity, access, ancestry/build | doctor / explicit source QC adapter |
| normalize | canonical summary statistics | science.normalize / source-specific adapters |
| cis_regions | full standardized cis region | science.region / validated coordinates/liftover |
| ld_qc | signed population-matched matrix | science.load_ld / reference-genotype QC |
| instruments | strong LD-selected instruments | science.instruments; legacy conditional signals kept distinct |
| harmonize | aligned SNPs + QC | science.harmonize; explicit ambiguity policy |
| discovery_mr | Wald/IVW + testing family | science.mr; no causal-proof claim |
| sensitivity | eligibility, Egger/median/Q | analysis.R MR + science.mr |
| regional_coloc | full overlap and prior sensitivity | analysis.R ABF; trait SD/case fraction |
| finemap | convergence/CS/pairwise evidence | analysis.R runsusie/coloc.susie |
| replication | independent validation definition | new project-specific evidence adapter needed |
| annotation | tissue/cell/eQTL evidence | science.annotate + CKD/IS reference adapters |
| score | approved allele-aligned individual score | science.score; genotype QC upstream |
| longitudinal | repeated phenotype model | analysis.R cohort; cohort-panel adapter needed |
| incident | validated incident endpoint model | analysis.R incident; endpoint adapter needed |
| evidence | all requested outputs + unresolved coverage | science.evidence / architecture coverage aggregation |
| report | evidence-backed research report | project report adapter needed |

## Ischemic-stroke locus-first catalog

| Stage | Purpose |
|---|---|
| acquisition | frozen BBJ/GIGASTROKE and functional-source identity |
| source_qc | checksum, ancestry/build and schema audit |
| normalize | canonical stroke GWAS representation |
| gwas_loci | BBJ L001-L004 discovery-locus contract |
| finemap | BBJ fine-mapping / credible sets |
| cross_ancestry | GIGASTROKE/EAS replication |
| molecular_coloc | eQTL/pQTL molecular convergence |
| mechanism | FGF5 / ALDH2 / SH3PXD2A / COL4A2 branches |
| celltype | frozen mouse/human cell-type evidence |
| human_annotation | GSE256493 author-validated vascular annotation |
| evidence | integrated locus-first evidence matrix |
| report | thesis/manuscript report contract |

The IS catalog intentionally does **not** require `discovery_mr`, `sensitivity`,
`score`, `longitudinal` or `incident`. Those belong to the CKD/pQTL-cohort
design and are not forced into the current stroke thesis merely for code reuse.
The shared engine is the reusable layer; the scientific DAG remains project-specific.

These are reusable contracts, not pre-bound claims of scientific completion.
Initialization leaves every enabled adapter unbound. In CKD, score/longitudinal/incident
are NOT_REQUESTED by default and may be enabled together. In the current IS thesis
catalog all 12 stages are required, including human annotation, so an unfinished
GSE256493 author-annotation layer keeps execution blocked rather than producing a
false complete state.

For CKD, full cis-region output remains separate from instrument output and coloc cannot
consume an instrument-only file. For IS, the equivalent guardrail is the frozen BBJ
locus definition: no cis-window/pQTL assumption is introduced merely to reuse code.
Across both projects ancestry/build, allele policy, replication definition and testing
family must be explicit before executable bindings are accepted.

## Commands

```bash
bash server/masteromics_structure.sh init --root /srv/is-analysis/masteromics_workspace
bash server/masteromics_structure.sh inspect /srv/is-analysis/masteromics_workspace/ckd/project.json
bash server/masteromics_structure.sh inspect /srv/is-analysis/masteromics_workspace/ischemic_stroke/project.json
```

Initialization creates only directories and metadata/configuration JSON, does
not inspect raw data, and refuses to overwrite existing project configs. It uses
only the Python standard library. Structure validation checks the exact catalog,
dependency order, required stages and optional-cohort prerequisites.

After data review and adapter implementation/binding:

```bash
python -m masteromics blueprint run /ABSOLUTE/PROJECT/project.json --plan
python -m masteromics blueprint run /ABSOLUTE/PROJECT/project.json --jobs 4
```

`--plan` prints structure/binding status without data validation or execution.
Run is blocked until every enabled stage has an explicit binding and scientific
policy. Binding schema: `argv` (list, no shell), `inputs` (paths), `outputs`
(engine TSV/JSON contracts). Every output must be owned by one stage, inside
the project workspace, and passed via its temporary `@out0`, `@out1`, … placeholder.
TSV contracts require columns and nonempty rows; JSON contracts require keys.
The existing central DAG engine supplies locking, hashing, environment capture,
validated atomic promotion, concurrency and resume. No placeholder success
files are generated to make an unimplemented stage pass.

## Completion truth

| State | Meaning |
|---|---|
| VALID structure | topology/contracts are structurally valid, not data-ready |
| UNBOUND_ADAPTER | implementation/config binding is missing; execution blocked |
| BOUND_NOT_DATA_VALIDATED | binding present, data/scientific readiness not established |
| NOT_REQUESTED | optional analysis excluded from current scope |
| NOT_EXECUTED | no scientific results established by the scaffold |

Existing CKD frozen-input comparison and raw-rebuild migration remain available
as separate commands; they are not silently run by this blueprint. IS full-source
migration, tissue evidence integration, approved cohort adapters and report
assembly remain explicit next tasks. Production-data review is deferred by user.

Tests use small generated fixtures for initialization/no-overwrite, dependency
gates, unbound execution refusal and binding contracts. Existing scientific
backend tests remain separate; CI success is not production data validation.

## 데이터 소스 연결

[DB 목록과 접근 계획](MASTEROMICS_DATABASES.md)을 기준으로 acquisition 입력을 선택한다. `python3 -S -m masteromics resources validate`와 `resources plan --ids ...`는 데이터 다운로드 없이 목록과 접근 조건을 검사한다.
