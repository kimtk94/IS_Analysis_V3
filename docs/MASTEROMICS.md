# MasterOmics 0.1 — central engine and CKD/IS migration

One code path now compiles any phenotype recipe into the same dependency graph.
Code stays in this Git repository; datasets and results live outside it. Legacy
scripts are unchanged. Initial presets specify CKD nine candidates and IS six
candidates, each with EUR discovery and EAS outcome analysis. These are candidate
reanalysis presets, not a completed proteome-wide screen.

## Server setup and execution

Use Python >=3.11. Create an environment and install `requirements-masteromics.lock`
before a run. R packages required: TwoSampleMR, coloc, susieR, jsonlite; for cohort:
lme4, survival. R packages are never installed by the runner. Actual R MR, coloc, SuSiE and cohort models passed synthetic integration CI
run 37104819367. Production CKD/IS/KoGES validation remains required.
CI records full installed R versions; production environment locking is still required.

Copy `.example.json` files to `ckd.json`, `ischemic_stroke.json`, `datasets.json`.
Fill every placeholder from observed data, including SHA256, trait SD/case fraction,
same-build gene cis intervals and ancestry-matched signed LD. Do not relabel a
GRCh37 outcome as GRCh38: validated liftover/reference allele normalization is an
upstream requirement. Dense per-locus LD is supported; genome-wide matrices are not.

```bash
cd /srv/is-analysis/IS_Analysis_V3
python3 -m masteromics inventory /srv/is-analysis/data /srv/is-analysis/input_inventory.json
python3 -m masteromics run projects/ckd.json projects/ischemic_stroke.json --registry projects/datasets.json --jobs 4 --plan
bash server/masteromics_run.sh projects/ckd.json projects/ischemic_stroke.json --registry projects/datasets.json --jobs 4
```

`--jobs 4` runs up to four independent stages concurrently within each project.
Projects run sequentially to limit memory. BLAS/R thread counts should be set during
server setup to prevent each job from consuming all CPUs. Never start two jobs for
the same output root. Locks prevent this. A failed project returns nonzero while
independent ready work may finish. Downstream stages are blocked.

## Dataset contract

One dataset is one trait/assay. Explicit mappings for `chr,pos,effect_allele,
other_allele,beta,se,pval,eaf,n` are required. Registry fields include ID, path,
checksum, optional download URL, ancestry, build, effect scale and trait type.
OR/log10P and archives require explicit upstream conversion; they are not guessed.
Current adapter handles delimited plain/gzip tables in memory. Stage full regional
inputs for very large studies; never use only significant SNPs for coloc. Multi-assay
files must be split or adapted before registration. Genome-wide streaming adapters
are not yet implemented.

LD consists of a square TSV with identical variant-key row and column order and a
variants TSV with `key,effect_allele` (the allele counted by signed correlation).
Key is `BUILD:CHR:POS:SORTED_ALLELE_PAIR`, e.g. `GRCh38:1:100:AG`. The matrix must be
finite, symmetric, PSD, diagonal one, and cover all requested SNPs. LD signs are
aligned to exposure alleles before fine-mapping. EUR exposure and EAS outcome use
separate panels. This cross-ancestry analysis is not independent exposure replication.
Gene intervals must document coordinate source and cis-window convention externally.

## Implemented analysis

Acquisition/checksum → full normalization → interval/significance/F filtering →
LD-based greedy clumping → conservative allele harmonization → Wald/fixed IVW/
multiplicative IVW → official TwoSampleMR median/Egger → regional coloc ABF with
three priors → signed-LD SuSiE/coloc.susie → annotation → final evidence.

Ambiguous palindromic SNPs are all dropped. Strand complements and indels are
intentionally unsupported in this first adapter. No silent fallback. Missing LD,
invalid builds, missing constant sample size, undefined quantitative trait SD,
no instruments, insufficient regional overlap, and nonconverged SuSiE stop the
relevant unit. No credible sets produce `UNRESOLVED_NO_CREDIBLE_SET`, not shared.
Current bootstrap methods use a fixed seed. MR is a statistical estimate, not proof
of a causal protein mechanism. Sample overlap, conditional-vs-marginal pQTL, LD
reference mismatch, assay-binding artifacts and phenotype transformations require
study-specific review. ST16 conditional independent signals do not replace full
marginal pGWAS for regional analyses.

Final outputs: `MASTER_MR_EVIDENCE.tsv` (BH within method across configured tests),
`MASTER_CANDIDATE_EVIDENCE.json` (all requested scientific output tables, source
checksums), `run_manifest.json`, `run_summary.json`, stage QC, logs and checkpoints.
NaN tests are excluded from BH. Declare the intended testing family when selecting
units; candidate-only BH does not represent proteome-wide FDR.

## Cohort boundary

Optional `cohort: {panel: ABSOLUTE_TSV, pcs: [PC1,PC2,...]}` executes an eGFR mixed
model with score×time, age, sex, PCs and random intercept/slope. Panel columns are
`id,time,egfr,score,age,sex,PC...`; repeated time rows/missing covariates fail.
A precomputed score panel or the dosage/weight scoring adapter below can be used.
Genotype QC, raw KoGES endpoint construction, eQTL integration, reverse MR,
Steiger and leave-one-out are not implemented in this release. Do not mark them
completed. Public training data without genotype cannot validate a genetic score.
A project requiring these can add external command stages using `dag`, with explicit
input hashes and TSV/JSON output contracts; a marker file alone is not evidence.

## Resume and validation

Same input/config/core/external script hashes and validated unchanged output → skip.
Input, configuration or code changes invalidate the stage and its downstream chain.
Each stage writes private temporary outputs and promotes only validated files.
A successful checkpoint is written after all outputs are promoted. No biological
empty outputs are currently promoted; they fail for review. Pipeline success means
all configured stages passed, not that unconfigured methods ran.

Python tests cover real synthetic effect 0.5, LD/build mismatches, allele swaps,
palindromic drops, output corruption, changed inputs/external code, header-only
rejection, dependency ordering and resume. R integration tests passed in GitHub CI: median/Egger effect 0.5, shared-locus
ABF and SuSiE posterior checks, full DAG/resume, LD order rejection, mixed slope
and Cox effect recovery. Production CKD/IS fits have not been run.

Development source: kimtk94/IS_Analysis_V3 commit 7b7311bf123c79128cfc936634dec72c705bd2a4.
Reviewed legacy MR, coloc/SuSiE entrypoints and workflow stubs; not all historical
code/backups/results or raw data. This new implementation remains separate from
legacy result interpretation.

### Additional cohort modules

`score DOSAGES WEIGHTS OUTPUT` computes an allele-aligned weighted genetic score.
Dosages must contain `id,key,dosage,effect_allele,other_allele,build`; weights contain
`key,beta,effect_allele,other_allele,build`. All weighted variants must be present for
every retained person. No silent imputation or participant dropping. Input genotype
QC, relatedness filtering and ancestry PCA must be completed upstream.

A cohort recipe can instead provide `dosages`, `weights`, `phenotypes`, `pcs` and
optional `incident_panel`; score and phenotype join then run automatically. Incident
panel columns: `id,followup,event,baseline_ckd,score,age,sex,PC...`. Endpoint definitions
and timing are externally prespecified; prevalent CKD is excluded. Cox output includes
HR and proportional-hazards diagnostics. Event derivation from raw KoGES and genotype
QC are still upstream adapters. These R cohort models passed synthetic CI; real KoGES fits remain unexecuted.

## Read-only server preflight

```bash
python3 -m masteromics doctor projects/ckd.example.json --registry projects/datasets.example.json --output /srv/is-analysis/masteromics_preflight.json
```

The example recipe intentionally reports BLOCKED_INPUTS until actual bindings are
provided. Doctor checks registry identity metadata, headers, cis coordinates, LD
paths and ancestry/build, pinned checksum syntax and R availability. It never
downloads or reads whole raw datasets. READY_FOR_SCIENTIFIC_GATES is configuration
readiness, not biological validation or a completed scientific run.

Latest verified CI: 6 Python tests plus 3 actual R integration tests and syntax/compile
checks. Commit 4d9b0c725377af33e244e6cbef8c57e25269af47.
