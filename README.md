# IS Analysis V3

Reproducible UKB-PPP EUR/EAS pQTL preparation and downstream-analysis scaffolding.
The repository deliberately separates ephemeral code (`CODE_ROOT`) from persistent
Google Drive data (`WORK_ROOT`). See [the migration and execution guide](docs/v3_migration_guide.md).

## Validate a fresh checkout

Self-contained code regressions do not require Synapse, Google Drive, the old
research checkout, or production GWAS/QTL downloads:

```bash
python3 -m pip install -r requirements-test.txt
bash scripts/codex_smoke_test.sh --fixture-only
```

Install base R (`Rscript`) to include the preparation and cache regression tests.
Without R, Python reports those tests as skipped. The **Execution foundation
regressions** GitHub Actions workflow installs both Python dependencies and base
R so those tests execute in CI.

Four additional historical-source integration suites authenticate files from
the pinned old research commit. Run them only after preparing that checkout:

```bash
git fetch --no-tags origin 0bf646fd6f29322d17da878253f37def24c01c74
git worktree add --detach legacy 0bf646fd6f29322d17da878253f37def24c01c74
bash scripts/codex_smoke_test.sh --with-legacy
```

`--with-legacy` fails immediately if required files are absent; their original
hash checks remain enforced. With no flag, smoke runs the self-contained suite
in a fresh clone and the full suite when `legacy/` is present, preserving the
existing prepared source-audit CI gate. To run only the historical checks, use
`python3 scripts/run_repository_tests.py --suite legacy`.

See [execution and deployment maintenance](docs/execution_foundation.md) for
output-path migration, result-root configuration, and queue recovery checks.

## Quick start (Colab)

```bash
export CODE_ROOT=/content/IS_Analysis_V3
export WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V3
export SOURCE_WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V2
cd "$CODE_ROOT"
bash scripts/setup_ido1_test_drive.sh
source "$WORK_ROOT/ido1_test.env"
python3 -m pip install -r requirements-test.txt
bash scripts/codex_smoke_test.sh
```

### Required manifest and gene coordinates

The batch runner does not discover UKB-PPP URLs or gene coordinates by itself.
It requires these two reviewed TSV inputs:

* `ukb_ppp_download_manifest.tsv`: one row per gene and ancestry. The runner
  accepts either its compact columns (`gene`, `ancestry`, `source_url`,
  `size_bytes`) or the production columns (`gene_symbol`, `ancestry`, `url`,
  `expected_size_bytes`). `source_file`, `synapse_id`, `sha256`, and `md5` are
  used for deterministic naming, authenticated downloads, and validation.
* `gene_coordinates_hg38.tsv`: gene locations with the columns `gene`, `chr`,
  `start`, `end`, and `genome_build` (`GRCh38`).

If the reviewed files already exist in the V2 work root, copy them into the V3
layout and load the generated environment variables with:

```bash
export CODE_ROOT=/content/IS_Analysis_V3
export WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V3
export SOURCE_WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V2

bash "$CODE_ROOT/scripts/setup_ido1_test_drive.sh"
source "$WORK_ROOT/ido1_test.env"
```

Before downloading, verify that both files exist, have the expected headers,
and contain IDO1 for both EUR and EAS:

```bash
test -s "$WORK_ROOT/data/metadata/ukb_ppp_download_manifest.tsv"
test -s "$WORK_ROOT/data/reference/gene_coordinates_hg38.tsv"

head -n 3 "$WORK_ROOT/data/metadata/ukb_ppp_download_manifest.tsv"
head -n 3 "$WORK_ROOT/data/reference/gene_coordinates_hg38.tsv"
awk -F '\t' 'NR == 1 || toupper($1) == "IDO1"' \
  "$WORK_ROOT/data/metadata/ukb_ppp_download_manifest.tsv"
awk -F '\t' 'NR == 1 || toupper($1) == "IDO1"' \
  "$WORK_ROOT/data/reference/gene_coordinates_hg38.tsv"
```

The manifest should show two IDO1 rows (`EUR` and `EAS`), while the coordinate
file should show one IDO1 GRCh38 row. If the V2 files do not exist, obtain and
review these inputs from the project data owner; do not substitute the synthetic
files under `tests/fixtures` for a real download.

Production URLs of the form `https://www.synapse.org/Synapse:syn...` are landing
pages, not archive download URLs. Install and authenticate the Synapse CLI before
running a production download; the runner uses the manifest's `synapse_id`:

```bash
python3 -m pip install -r "$CODE_ROOT/requirements.txt"
synapse login
synapse get syn52363617 --downloadLocation /tmp/synapse_download_check
```

Do not put a Synapse password or authentication token in the manifest or Git.

Download and decompress a bounded test sample for the 15-gene batch containing
IDO1 without running the preparation step. IDO1 is capped at 500,000 lines and
each other gene is capped at 1,000 lines (both counts include the header):

```bash
python3 scripts/ukb_ppp_batch_manifest_runner_fast.py \
  --base "$WORK_ROOT/data/rawdata/pqtl/selected_targets" \
  --qc-dir "$IDO1_TEST_QC_DIR" --outdir "$IDO1_TEST_OUTDIR" \
  --standardized-dir "$IDO1_TEST_STANDARDIZED_DIR" \
  --instrument-dir "$IDO1_TEST_INSTRUMENT_DIR" \
  --download-manifest "$WORK_ROOT/data/metadata/ukb_ppp_download_manifest.tsv" \
  --gene-coordinate-file "$WORK_ROOT/data/reference/gene_coordinates_hg38.tsv" \
  --batch-size 15 --focus-gene IDO1 --focus-max-file-lines 500000 \
  --other-max-file-lines 1000 --test-data-only --stop-on-error
```

`--test-data-only` first checks the expected path under `--base` and reuses an
existing raw file, including one without manifest size or checksum metadata.
A missing file is downloaded; when validation metadata is available, an
existing file that fails validation is also replaced. Raw files may be
gzip-compressed; `--test-data-only` reads their decompressed contents
(equivalent to `zcat`) before applying the line limits.
If all expected test TSVs already exist, the batch is treated as sampled without
checking or downloading raw files. Otherwise, all raw files for the selected
batch are checked or downloaded before missing test samples are written;
existing test TSVs are preserved rather than overwritten.
The resulting plain TSV files are written under
`$IDO1_TEST_OUTDIR/test_data/{EUR,EAS}`. To retain only the downloaded archives
without creating these samples, use `--download-only` instead.

```bash
python3 -u "$CODE_ROOT/scripts/ukb_ppp_batch_manifest_runner_fast.py" \
  --base "$WORK_ROOT/data/rawdata/pqtl/selected_targets" \
  --qc-dir "$IDO1_TEST_QC_DIR" \
  --outdir "$IDO1_TEST_OUTDIR" \
  --standardized-dir "$IDO1_TEST_STANDARDIZED_DIR" \
  --instrument-dir "$IDO1_TEST_INSTRUMENT_DIR" \
  --download-manifest "$WORK_ROOT/data/metadata/ukb_ppp_download_manifest.tsv" \
  --gene-coordinate-file "$WORK_ROOT/data/reference/gene_coordinates_hg38.tsv" \
  --batch-size 15 --focus-gene IDO1 \
  --focus-max-file-lines 500000 --other-max-file-lines 1000 \
  --test-data-only --stop-on-error
```

Test samples are recreated under
`$IDO1_TEST_OUTDIR/test_data/{EUR,EAS}/<GENE>.tsv`; current task states are
rewritten to `batch_progress.tsv` and successful samples have status `sampled`.

Omit both `--download-only` and `--run` for a safe dry run. Use `--run` instead
to download the selected batch and run the R preparation step for each archive.

Prepared batch outputs are isolated by ancestry, batch, gene, and source. The
execution plan and progress TSV record their exact paths; downstream consumers
should use those paths rather than assuming a shared `<GENE>.tsv`. See the
[output migration guide](docs/execution_foundation.md#prepared-output-identity).

Runtime outputs and raw archives are ignored by Git. Only synthetic fixtures are
committed.
