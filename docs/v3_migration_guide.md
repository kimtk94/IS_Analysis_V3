# V3 migration and execution guide

## Storage contract

| Variable | Purpose | Persistent |
|---|---|---|
| `CODE_ROOT` | Git clone and committed fixtures | No |
| `WORK_ROOT` | metadata, reference, raw staging, QC, and results | Yes |
| `SOURCE_WORK_ROOT` | read-only V2 metadata/reference source | Yes |

No raw UKB-PPP/GWAS summary statistics or derived results may be stored beneath
`CODE_ROOT`. `scripts/setup_ido1_test_drive.sh` creates the V3 Drive layout and
copies only missing manifest and coordinate files, without overwriting reviewed
V3 files.

## Pipeline stages

1. **Setup:** copy `ukb_ppp_download_manifest.tsv` and
   `gene_coordinates_hg38.tsv` from V2 when absent.
2. **Plan:** pair genes represented in both EUR and EAS, sort them, form stable
   15-gene batches, and optionally select the single batch containing IDO1.
3. **Prepare:** validate/stage selected sources and invoke
   `scripts/01_prepare_exposure_fast.R` per ancestry/source. Canonical summaries
   and filtered cis instrument candidates are separate outputs.
4. **Analyse:** run ancestry-matched LD/fine-mapping, causal checkpoints,
   GIGASTROKE adaptation, and brain-eQTL colocalization through the dedicated
   workflow entry points.

The runner always writes `batch_manifest.tsv`, `execution_plan.tsv`, and
`run_metadata.json`. `--download-only` downloads without preparing,
`--test-data-only` downloads and writes bounded complete-line test samples, and
`--run` downloads and prepares. Without any of these options, no archives are
downloaded or processed. Raw cleanup is intentionally not implemented in
focused testing.

## Provenance checklist

Record remote URL, commit SHA, code/work roots, input SHA-256 checksums, genome
build, batch/focus limits, statistical thresholds, timestamps, and cleanup state.
Legacy V2 QC may be copied under `results/qc/legacy_v2_reference`, but must never
be treated as completion state for a V3 batch.
