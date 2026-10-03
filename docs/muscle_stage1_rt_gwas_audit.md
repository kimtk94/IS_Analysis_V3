# MUSCLE Stage 1 — Resistance-training GWAS audit

## Purpose
Audit the genetic discovery layer for skeletal-muscle trainability before any MR, colocalization, or functional prioritization.

## Studies registered
- Yang et al. 2024, *Physiological Genomics*, PMID 38881426: 440 inactive adults, 12-week RT/HIIT, rectus femoris muscle-thickness response. Eleven RT and eight HIIT lead variants were reported at P < 1e-5. Full data require approved request to the Population Health Data Archive.
- Gu et al. 2026, *Journal of Cachexia, Sarcopenia and Muscle*, PMID 42455518: 187 young Asian adults, 12-week RT, DXA delta lean-body mass. Nine lead rsIDs are audited from the open article tables.

## Design rules
1. Keep reported and recomputed significance separate.
2. Tier A: P < 5e-8.
3. Tier B: 5e-8 <= P < 1e-5.
4. Tier C: P >= 1e-5 or other supporting evidence.
5. Missing P/alleles/build fields remain unresolved; never impute them from gene proximity.
6. Do not treat 2024 and 2026 cohorts as independent replication until participant overlap is explicitly excluded.
7. MR/coloc are gated on full, harmonizable summary statistics.

## CLI
```bash
python3 scripts/muscle/run_muscle_stage1_rt_gwas_audit.py \
  --outdir /path/to/results/muscle/stage1_rt_gwas_audit
```

Use `--offline` for deterministic smoke tests. Use `--xml-input` to test a saved Europe PMC XML file.

## Outputs
- `MUSCLE_STAGE1_STUDY_REGISTRY.tsv`
- `MUSCLE_STAGE1_2026_VARIANTS.tsv`
- `MUSCLE_STAGE1_2026_RAW_TABLE.tsv` when table parsing succeeds
- `MUSCLE_STAGE1_QC_FLAGS.tsv`
- `MUSCLE_STAGE1_SUMMARY.json`
- `MUSCLE_STAGE1_README.md`
- `raw/PMC13371584.xml` when downloaded
