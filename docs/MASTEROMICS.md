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

## Server regression against existing CKD results

Run in an isolated checkout of `feat/masteromics-central-engine`. The adapter
reads existing results; it writes only a new, separate output directory. No
packages are installed and no legacy results are overwritten.

```bash
bash server/masteromics_regression.sh --genes SDCCAG8 GSTA3 --stages mr coloc
# Once the first run passes, compare all nine CKD candidates including SuSiE:
bash server/masteromics_regression.sh --stages mr coloc susie
```

Defaults: `/srv/is-analysis/results/ckd/stage1`, `stage2b_coloc`,
`stage2c_susie`; LD: `/srv/is-analysis/data/ckd/stage2c_ld/ld`.
Override with `--root`, `--stage1-root`, `--stage2b-root`, `--stage2c-root`,
`--ld-root`. Set `MASTEROMICS_PYTHON` to the scientific Python interpreter;
otherwise it checks the checkout, the existing `/srv/is-analysis/IS_Analysis_V3`
checkout and `/srv/is-analysis` for `.venv-ckd/bin/python`, then uses `python3`.
Rscript must resolve to the intended existing R environment. The wrapper respects
`CKD_R_LIB`/`R_LIBS_USER`, otherwise reuses `/srv/is-analysis/.Rlib` if present.
Set `MASTEROMICS_REGRESSION_OUT` to a fresh output directory when needed.

`REGRESSION_COMPARISON.tsv` contains baseline and new values, absolute deltas,
explicit tolerances, and PASS/DIFFERENCE for Wald/IVW beta, SE, P, IVW Q, SNP
counts, coloc H0–H4, SuSiE maximum H4 and credible-set counts.
Default acceptance is `abs(delta) <= atol + 1e-6 * abs(legacy)`: MR atol 1e-8,
posterior atol 1e-4; counts must match exactly. `REGRESSION_SUMMARY.json` records
input/baseline/LD/code hashes, coverage and errors. Exit 0 means all requested
comparisons passed. Missing genes, missing inputs, failed fits or unresolved
SuSiE produce nonzero exit and an INCOMPLETE report. DIFFERENCE needs review;
it is never automatically accepted. SuSiE holds a dense matrix in memory;
more than 8,000 variants requires explicit `--max-ld-variants` override.

The numerical replay freezes the original harmonized instrument set and full
regional inputs, and retains legacy median-rounded sample sizes and coloc
outcome sdY estimation. This proves numerical parity for prepared inputs,
not end-to-end equivalence from raw archives/genotypes. Central production
recipes still require an explicit phenotype sdY. `POLICY_DELTA.tsv` separately
records coloc changes after removing palindromic SNPs; these policy changes do
not count as same-input numerical failures. New clumping/instrument policies
must be reviewed separately from the legacy primary anchor Wald result.

The adapter currently supports the legacy **CKD** schemas, EUR regional/LD
analysis and EUR/EAS Stage1 MR. It rejects IS and EAS regional replay rather
than interpreting their files as CKD inputs. IS needs an explicit schema and
baseline adapter before server validation can be claimed. Real server CKD
results have not yet been tested by the assistant; synthetic legacy/new
parity is tested in CI.
