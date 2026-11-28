# Reviewed source and execution cleanup — 2026-11-24

This change consolidates reviewed code on top of main
`2e5c4d36a91acd14ad8f31fade62d48ef5dd670a`. It prepares repeatable code checks and
preserves the provenance of existing evidence. It does not generate new research
results or change the interpretation of frozen analyses.

## Source provenance

| Component | Reviewed source | Change in this consolidation |
| --- | --- | --- |
| UKB preparation, CKD cache/ID parsing, dashboard input handling | [PR #331](https://github.com/kimtk94/IS_Analysis_V3/pull/331), commit `f0ce8cf8e30e4762622649be57ac4f85fc54ec4b` | Select the execution fixes and small regressions; integrate the numeric QC requirements reviewed in [PR #18](https://github.com/kimtk94/IS_Analysis_V3/pull/18). |
| Ocular P1/P4/P5 dependency chain | [PR #253](https://github.com/kimtk94/IS_Analysis_V3/pull/253), commit `0d39a9b41aa500d39327fdd8709f9ccb179f8502` | Restore the three missing Python modules and P5 R source; report the actual P5 classification cutoff. |
| CKM workflow source | [PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284), commit `4aff051003c6a65421516e69956c905498848a76` | Replace the moving research branch in 35 checkout steps across 31 workflows with this exact commit. |

The three restored ocular Python files are byte-identical to that reviewed
source. The P5 R change replaces the hard-coded `classification_score_threshold`
of 0.75 with `score_cutoff`; the statistical model is unchanged. The dependency
chain is `ocular_p4_fullcohort.py` → `ocular_p4_pseudobulk_pilot.py` →
`ocular_p1_integrity.py`.

Harness deployment and worker work were abandoned by project direction and are
outside this consolidation.

## UKB preparation and QC contract

The batch runner supplies the actual batch identifier and a deterministic source
identifier to the R preparation CLI. For each output root, batch outputs use:

```text
<ancestry>/<batch_id>/<gene>/<source_id>.tsv
```

Standardized data, instruments and legacy output each use their respective root.
Successful QC is next to the standardized file as `<source_id>.qc.tsv`. A direct
CLI call without `--source-id` retains the `<ancestry>/<batch_id>/<GENE>.tsv`
layout and `<GENE>.qc.tsv` sidecar. The execution plan and progress table expose
all output paths, `qc_path`, and `failed_qc_path`.

The preparation step selects the gene chromosome's unique archive member and
handles nested gzip compression. A missing/ambiguous chromosome member is an
error. It accepts original UKB-PPP headers (`CHROM`, `GENPOS`, `ID`, `ALLELE0`,
`ALLELE1`, `A1FREQ`, `LOG10P`) and selects ADD rows. It does not silently substitute
a different chromosome or assay.

Numeric rules:

- A supplied direct P column takes precedence over LOG10P. Invalid direct P is
  excluded rather than replaced from another column.
- Finite, nonnegative LOG10P is accepted. Zero produces P=1; finite numeric
  underflow to P=0 is counted. Negative, nonfinite and missing LOG10P is excluded.
- EAF is optional when the entire column is absent. If the column is present,
  missing/nonfinite values or values outside [0, 1] are invalid.
- Invalid coordinates/alleles, nonfinite effects, nonpositive or nonfinite SE,
  invalid P and nonfinite F statistics are excluded. F overflow is invalid even
  when beta and SE themselves are finite.

QC includes source path/member, identifiers, timestamp, bounded input limits,
row counts and reasons for exclusions. Its denominators satisfy:

```text
input_rows = non_additive_rows + off_gene_chromosome_rows + evaluated_rows
output_rows = evaluated_rows - invalid_any
```

Individual reason counts can overlap; `invalid_any` counts their union.
`PREPARED` describes the bounded input actually read. It is not evidence that the
whole source archive or downstream scientific analysis is complete.

An all-invalid attempt writes `<source_id>.failed.qc.tsv` (or
`<GENE>.failed.qc.tsv`) with `NO_VALID_ROWS` and exits unsuccessfully. It leaves any
previous successful QC and all three previous TSVs intact. Thus the normal QC
describes the most recent successful preparation, while the failed sidecar
describes the most recent all-invalid attempt. A later success preserves that
failure diagnostic; use the progress status and timestamps to identify the
current attempt. Parse/schema failures can occur before a row-level QC sidecar
can be produced.

## CKD and dashboard behavior

Related-ID parsing accepts IID-only lists, FID/IID pairs and headers explicitly
identifying IID. Inconsistent or ambiguous column layouts fail instead of
silently selecting a wrong ID column.

SuSiE cache fingerprints cover the actual ordered effect/variance inputs,
sample sizes, allele metadata, oriented LD matrix, model options and R/package
versions. Only a matching, converged fit is reusable. Unverified legacy caches,
input mismatches, unreadable caches and failed/nonconverged fits cause a cache
miss. This means a future authorized analysis may recompute an old fit; this
cleanup itself does not run a production SuSiE analysis. The base-R regression
uses a synthetic fitter, so it verifies cache behavior rather than the scientific
model's estimates.

Dashboard publication accepts `--results-root`; its default is
`${IS_ANALYSIS_ROOT:-/srv/is-analysis}/results`. The existing `--source` option
continues to refer to a work root containing `results/`. Missing/unreadable input
produces pending state and empty metrics. Historical constants are no longer
published as live results when input files are absent. The script's optional
publication action was not exercised as part of this cleanup.

## Workflow execution boundaries

The 31 CKM source workflows now use fixed source commit
`4aff051003c6a65421516e69956c905498848a76`. Updating a workflow no longer triggers
its original-data analysis on push. `workflow_dispatch` remains available, as do
the two existing owner-checked issue-comment command entry points. The lead
dispatcher's PR event runs only its offline tests; original source jobs require
their existing explicit inputs or owner command.

The ocular workflow runs small source-import, count-conservation, metadata and
R-summary checks on a PR. Its raw-UMI reconstruction and edgeR refit require
`workflow_dispatch` and a successful offline validation job.

**Remaining ocular prerequisite:** the manual comparison still retrieves its
baseline from GitHub run `41074586446`. Saving the archive separately preserves
the evidence but does not change that workflow's download source. After the
original artifact expires, a verified accessible baseline must be restored and
the retrieval step updated before dispatching the comparison. Do not re-run the
original matrices merely to discover that the baseline download has expired.

## Existing aggregate evidence preserved outside Git

The original artifact ZIPs below were preserved byte-for-byte on 2026-11-24.
Their sizes and SHA256 values were checked against GitHub's artifact metadata,
ZIP integrity was checked, and each internal file has a SHA256 entry in the
preservation bundle's manifest. No full source datasets or participant-level
matrices are included in that bundle. The archive is a separate project
deliverable, not a tracked repository result file.

| Aggregate artifact | Original run | Artifact ID | Original expiry (UTC) |
| --- | --- | --- | --- |
| `ocular-p5-edger-inference` | [41074586446](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446) | 12987425927 | 2026-11-29 14:58:33 |
| `ocular-p4-full-qc` | [41069563630](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41069563630) | 12984842522 | 2026-11-29 14:17:37 |
| `hcc-ng-stage24-patient-spearman-aggregate` | [41173779578](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173779578) | 13032152820 | 2026-11-30 05:46:19 |
| `ocular-gse313736-technical-confounding-identifiability` | [41430965326](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41430965326) | 13133928461 | 2027-01-01 10:42:49 |
| `is-EAS-GIGASTROKE-publisher-MD5-GWAS-EAF-completeness` | [43858822105](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43858822105) | 14224434625 | 2027-01-15 12:16:30 |

Bundle: `IS_Analysis_existing_results_20261124.zip` (941,795 bytes).

SHA256: `4493bd3fa4692fbe7968654eaa6920ce508763ce63f9f43d27ab26a7329ac664`.

Preservation does not upgrade the evidence: ocular treatment/library confounding
remains unresolved; the HCC correlation summaries remain descriptive; and the
IS source-completeness artifact remains a source audit. The IS EAF JSON existed
even though a subsequent report-formatting step failed in that workflow.

## Local and CI verification

```bash
python3 -m pip install -r requirements-test.txt
bash scripts/codex_smoke_test.sh
```

With no `legacy/` directory, the default smoke command runs the self-contained
suite. It uses only small fixtures. Install base R to execute the R-backed
preparation, synthetic cache and P5 metadata regressions; otherwise their skip
messages remain visible. The dedicated foundation CI installs Python test
dependencies and base R and checks `Rscript --version` before running the suite.
Ocular source tests report an explicit skip if NumPy is absent in an unrelated
historical workflow; foundation and ocular CI install NumPy and execute them.

Historical source verification remains a distinct integration gate. The existing
IS source workflows prepare a second checkout of this same repository at
`0bf646fd6f29322d17da878253f37def24c01c74` under `legacy/`. With that pinned checkout
available, use:

```bash
bash scripts/codex_smoke_test.sh --with-legacy
```

The explicit integration mode fails when its original inputs are missing; it
does not synthesize substitutes or report that source verification passed.
The default smoke command also selects the full suite when `legacy/` is present.
No scientific completion claim follows from a green code regression suite.
