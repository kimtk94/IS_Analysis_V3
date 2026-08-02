# IS Analysis V3

Reproducible UKB-PPP EUR/EAS pQTL preparation and downstream-analysis scaffolding.
The repository deliberately separates ephemeral code (`CODE_ROOT`) from persistent
Google Drive data (`WORK_ROOT`). See [the migration and execution guide](docs/v3_migration_guide.md).

## Quick start (Colab)

```bash
export CODE_ROOT=/content/IS_Analysis_V3
export WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V3
export SOURCE_WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V2
cd "$CODE_ROOT"
bash scripts/setup_ido1_test_drive.sh
source "$WORK_ROOT/ido1_test.env"
bash scripts/codex_smoke_test.sh
```

### Required manifest and gene coordinates

The batch runner does not discover UKB-PPP URLs or gene coordinates by itself.
It requires these two reviewed TSV inputs:

* `ukb_ppp_download_manifest.tsv`: one row per gene and ancestry, with the
  required columns `gene`, `ancestry`, and `source_url`. Optional `sha256` and
  `size_bytes` columns enable download validation.
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

Download the 15-gene batch containing IDO1 without running the preparation step:

```bash
python3 scripts/ukb_ppp_batch_manifest_runner_fast.py \
  --base "$WORK_ROOT/data/rawdata/pqtl/selected_targets" \
  --qc-dir "$IDO1_TEST_QC_DIR" --outdir "$IDO1_TEST_OUTDIR" \
  --standardized-dir "$IDO1_TEST_STANDARDIZED_DIR" \
  --instrument-dir "$IDO1_TEST_INSTRUMENT_DIR" \
  --download-manifest "$WORK_ROOT/data/metadata/ukb_ppp_download_manifest.tsv" \
  --gene-coordinate-file "$WORK_ROOT/data/reference/gene_coordinates_hg38.tsv" \
  --batch-size 15 --focus-gene IDO1 --focus-max-bytes 20000000 \
  --other-max-file-lines 1000 --download-only --stop-on-error
```

For a bounded test dataset instead of archive-only downloads, replace
`--download-only` with `--test-data-only`. The selected target gene is written
up to `--focus-max-bytes` (20 MB below), and every other gene in its batch is
written up to `--other-max-file-lines` (1,000 lines including the header):

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
  --focus-max-bytes 20000000 --other-max-file-lines 1000 \
  --test-data-only --stop-on-error
```

Test samples are recreated under
`$IDO1_TEST_OUTDIR/test_data/{EUR,EAS}/<GENE>.tsv`; current task states are
rewritten to `batch_progress.tsv` and successful samples have status `sampled`.

Omit both `--download-only` and `--run` for a safe dry run. Use `--run` instead
to download the selected batch and run the R preparation step for each archive.

Runtime outputs and raw archives are ignored by Git. Only synthetic fixtures are
committed.
