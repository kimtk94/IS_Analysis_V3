# Execution and deployment maintenance

This guide covers code execution, output identity, and operational reporting.
Scientific approval remains governed by the source-specific IS and CKD protocols.
The P3A evidence snapshots and their original source hash pins are unchanged by
these maintenance changes.

## Prepared output identity

The batch runner records `source_id`, `standardized_path`, `instrument_path`, and
`legacy_path` in both `execution_plan.tsv` and `batch_progress.tsv`. For `--run`,
each output root now uses:

```text
<output-root>/<ANCESTRY>/<batch_id>/<GENE>/<source_id>.tsv
```

`source_id` is derived from the manifest source identity. EUR/EAS and separate
assays therefore have separate files even when the gene symbol is the same.
Use the recorded paths when consuming outputs; an old `<GENE>.tsv` is not a
reliable reference to a particular ancestry or assay. Existing output files are
not automatically renamed, deleted, or treated as newly prepared results.

The direct R CLI accepts the original required arguments, plus explicit
`--batch-id` and `--source-id`. The runner supplies both. Calls without a source
ID retain the simple gene filename; instrument output is still separated by
ancestry. Batch execution is the preferred interface for multiple assays.

The preparation step recognizes the UKB-PPP `CHROM`, `GENPOS`, `ID`, `ALLELE0`,
`ALLELE1`, `A1FREQ`, `BETA`, `SE`, and `LOG10P` columns, including conversion from
negative-log10 p-values. It selects the target chromosome from tar archives and
decodes gzip members. Ambiguous archive members or missing required columns
stop preparation instead of selecting an arbitrary first file.

## Test modes

Install `requirements-test.txt` and base R to run all self-contained Python/R
regressions. `bash scripts/codex_smoke_test.sh --fixture-only` selects that suite.
It uses synthetic inputs and small evidence already tracked in the repository;
it does not authenticate a new publisher download or establish scientific
eligibility.

Four historical-source integration modules require a separate checkout at
`legacy/` from commit `0bf646fd6f29322d17da878253f37def24c01c74`:

- `test_is_24_original_lead_review.py`
- `test_is_gws_eas_24_original_leads.py`
- `test_is_p3a_36_group_recovery.py`
- `test_is_p3a_chr10_legacy_original_signature.py`

Use `--with-legacy` to run the full suite. A missing input fails that mode with
exit code 2 before running tests. The historical tests and production hash pins
remain intact. Existing source-audit workflows that prepare `legacy/` and call
smoke without an option continue to run the full suite. See the root README for
the pinned-checkout commands.

## Results and dashboard roots

The dashboard publisher accepts `--results-root`, the directory containing
`ckd/`, `is/`, `muscle/`, and other result subdirectories. Its default is
`${IS_ANALYSIS_ROOT:-/srv/is-analysis}/results`. The legacy `--source` option
means a work root whose `results/` subdirectory should be read.

Example local generation, without pushing the dashboard repository:

```bash
python3 server/publish_research_dashboard.py \
  --results-root /srv/is-analysis/results \
  --dashboard-repo /srv/is-analysis/brain_research_mr_scrna_seq
```

Missing or unreadable metrics produce an explicit input status and no fallback
scientific numbers. A generated manifest describes the files it could read;
it does not approve their biological interpretation. Inspect the manifest
before using the separate `--git-push` option.

## Queue installation and recovery

Follow [the pinned worker installation guide](../deploy/is-harness/README.md).
All worker modules and units must come from the same reviewed commit. The
installer validates that package before any timer is enabled.

The worker-update workflow separates validation from installation: pushes and
pull requests run tests only. A manual `operation=install` transfers the exact
commit's verified bundle, installs all modules and units, and reloads systemd
definitions without starting services or enabling timers. It holds while a
queue/AUTO_NEXT service, timer, related process, or execution lock is active.
The server research checkout is not updated as part of this home-directory
installation; its revision must be reconciled separately.

For diagnosis, select `audit=worker-status` in `IS Server Read-Only SSH Audit`.
The streamed check reports installed hashes/receipt, checkout revision, systemd
state, existing locks, process counts, bounded evidence freshness, and categorized
journal failures. It does not execute the worker, expose credentials or raw logs,
or write server state. Missing or unreadable observations remain explicit holds.

A stale queued or running request produces a visible hold instead of a green
"refresh succeeded" impression. A hold does not prove that a server process has
stopped. Check service status, installed revision, worker logs, execution locks,
and existing outputs before deciding whether a request should be retried.
The maintenance code does not automatically close or requeue old Issues.

For a stranded request such as #244, collect the read-only runtime report first,
then reconcile the request with any running process, execution lock, and existing
result/completion evidence. A queued request's age is based on its creation time;
editing or relabeling the same old Issue does not make it runnable under the
stale-request rule. Only after the actual execution outcome is known should an
operator resolve the old request and decide whether a new request is required.

## Cached fits and reference checks

SuSiE fit caches carry an input fingerprint covering ordered SNP/effect metadata,
the oriented LD matrix, sample sizes, model options, and R/package versions.
Unchanged inputs may reuse a converged fit. Legacy, corrupt, nonconverged, or
mismatched caches are recomputed; the first rerun after this migration can
therefore take longer. Cache metadata is an execution check, not a scientific
approval, and existing fit RDS files retain the normal fit-object format.

PLINK sample exclusion accepts IID-only or FID/IID lists and an explicit IID
header. Malformed multi-column data is rejected. FASTA lookups include reference
length in the cache key so a SNP and indel at one coordinate are validated
independently.

## Container configuration

Runtime environment files, keys, raw data, and analysis outputs are excluded
from the Docker build context. Supply credentials at runtime. A previously built
image is not changed by adding `.dockerignore`; rebuild the image when deploying
the reviewed version and handle old image copies according to the deployment's
credential-retention policy.
