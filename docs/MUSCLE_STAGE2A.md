# MUSCLE Stage 2A — GTEx skeletal-muscle variant-to-gene mapping

## Goal
Take the 20 primary RT loci from Stage 1 v1.1 and independently test whether
each variant has significant regulatory evidence in human skeletal muscle.

## GTEx API rule
GTEx API v2 separates rsID lookup from association lookup:

1. `/dataset/variant` accepts `snpId` and resolves the dbSNP rsID to the
   dataset-specific GTEx `variantId`.
2. `/association/singleTissueEqtl` and
   `/association/singleTissueSqtl` are queried using the GTEx `variantId`.

Stage 2A v2.2 therefore never sends the source rsID directly as the association
`variantId`.

## Resources
- Ensembl REST variation endpoint: rsID -> current GRCh38 mapping.
- GTEx Portal API v2.
- Dataset: `gtex_v10`.
- Tissue: `Muscle_Skeletal`.
- Significant single-tissue eQTL and sQTL endpoints.

## Coordinate audit
The source papers report coordinates that may be based on GRCh37. Stage 2A stores:
- source chromosome/position
- Ensembl GRCh38 chromosome/position
- GTEx V10 `variantId` (b38)
- GTEx `b37VariantId` when available

Therefore a source-coordinate mismatch with GRCh38 is not automatically a
failure. `source_position_matches_gtex_b37` explicitly tests whether the
published coordinate corresponds to the GTEx b37 representation.

## Unresolved variants
A dbSNP rsID absent from the selected GTEx cohort is recorded as
`SKIPPED_NO_GTEX_VARIANT`. It is not counted as an API error and should later
be evaluated through ancestry-matched LD proxies in Stage 2B.

## Outputs
- `MUSCLE_STAGE2A_VARIANT_NORMALIZATION.tsv`
- `MUSCLE_STAGE2A_EQTL.tsv`
- `MUSCLE_STAGE2A_SQTL.tsv`
- `MUSCLE_STAGE2A_CANDIDATE_GENES.tsv`
- `MUSCLE_STAGE2A_SUMMARY.json`
- raw API JSON for provenance

## Interpretation
A zero significant GTEx association after a successful variantId query means no
significant precomputed association was returned for that exact variant in
skeletal muscle. It does not establish absence of regulatory activity.

## Next
Stage 2B expands each lead locus using 1000 Genomes EAS LD. This is especially
important for source variants not represented directly in GTEx and for
regulatory signals carried by high-LD proxies.
