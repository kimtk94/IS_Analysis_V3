# IS Broad Discovery — real genotype-based GWS clump screen (2026-10-09)

## Scope

After expanding 1KG EAS 504-sample genotype reference for seven GWS-containing coordinate components, we performed reference-LD clumping on each **source GWAS × region** having `P <= 5e-8` evidence. Cohort overlap across GIGASTROKE AS and subtype datasets means these comparisons must **not** be treated as independent replication.

- Input: `expanded_reference_v2/IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz`, GRCh37 genomic REF/ALT matched and effect-oriented.
- Signal index p threshold: 5e-8.
- Secondary variant p threshold: 1e-4.
- LD index threshold: `r² >= 0.1`.
- LD window: ±1000kb.
- LD genotype: EAS 1000 Genomes, n=504; four original BBJ region panels and three newly acquired panels, plus expanded G0015 and G0022.
- PLINK2 `--clump` completed successfully for 15 significant study-group combinations (others had no GWS variants).
- **18 index clumps** in total across those 15 combinations. This total double-counts signals replicated or shared between related GWAS reports.

## Apparent multiple index clumps

| GWAS source | Region | PLINK2 index clumps | Lead index IDs |
|---|---|---:|---|
| GIGASTROKE EAS AS GCST90104544 | G0022 (chr12) | 2 | `12:111629389:C:A`, `12:112930475:T:C` |
| GIGASTROKE EAS AIS GCST90104545 | G0022 (chr12) | 3 | `12:112241766:G:A`, `12:113031474:G:A`, `12:110675363:C:T` |

The remaining 13 significant study-region combinations each produced 1 index clump at the chosen threshold.

**Important:** These are **marker clumps**, not confirmed statistically independent causal effects, signals/credible sets, genes or proteins. Their count is highly sensitive to LD r² and clump radius, reference ancestry, sample size and differences between AIS/AS phenotypes.

## Scientific next step

G0022 should be prioritized for **conditional analyses and multi-signal SuSiE-RSS fine mapping** using the full, allele-harmonized GWAS locus and EAS reference signed LD. Prioritize after checking allele frequencies, LD matrix PSD, effective N, region boundaries, conditional association persistence and locus specificity. Do not use separate GIGASTROKE AS and AIS studies as if independently replicated cohorts.

## Code and output

```bash
python3 scripts/is/run_gws_ld_clump_screen.py
python3 scripts/is/run_gws_ld_clump_screen.py --execute
python3 -m unittest discover -s tests -p 'test_is_*' -v
```

Results:
- `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/clump_screen_v1/IS_GWS_LD_CLUMP_SUMMARY.json`
- `IS_GWS_CLUMP_EXECUTION_MANIFEST.tsv`
- `IS_GWS_LD_CLUMP_RESULTS.tsv`
- Individual `*.clumps`, PLINK2 logs and association inputs retained for audit.

No canonical CKD or historical IS files were overwritten.
