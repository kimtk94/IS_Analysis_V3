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

Omit both `--download-only` and `--run` for a safe dry run. Use `--run` instead
to download the selected batch and run the R preparation step for each archive.

Runtime outputs and raw archives are ignored by Git. Only synthetic fixtures are
committed.
