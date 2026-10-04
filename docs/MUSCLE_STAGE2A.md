# MUSCLE Stage 2A — GTEx skeletal-muscle variant-to-gene mapping

## Goal
Take the 20 primary RT loci from Stage 1 v1.1 and independently test whether
each variant has significant regulatory evidence in human skeletal muscle.

## Public resources
- Ensembl REST variation endpoint for rsID -> GRCh38 normalization.
- GTEx Portal API v2.
- Dataset: `gtex_v10`.
- Tissue: `Muscle_Skeletal`.
- Significant single-tissue cis-eQTL and sQTL endpoints.

## Why rsID is the anchor
Source-paper genomic positions may reflect different genome-build contexts.
Stage 2A therefore re-resolves each locus from rsID rather than trusting the
published position as the primary join key.

## Outputs
- `MUSCLE_STAGE2A_VARIANT_NORMALIZATION.tsv`
- `MUSCLE_STAGE2A_EQTL.tsv`
- `MUSCLE_STAGE2A_SQTL.tsv`
- `MUSCLE_STAGE2A_CANDIDATE_GENES.tsv`
- `MUSCLE_STAGE2A_SUMMARY.json`
- raw API JSON for provenance/reproducibility

## Candidate score
A prioritization heuristic only:
- carry forward best Stage 1 evidence points
- +1 source-paper mapped gene
- +3 significant GTEx skeletal-muscle eQTL
- +2 significant GTEx skeletal-muscle sQTL

This is not a causal posterior.

## Interpretation
A zero GTEx result means no significant precomputed association was returned by
the selected API query. It does not establish absence of regulatory activity.

## Next
Stage 2B expands each lead locus using ancestry-matched 1000 Genomes EAS LD
before additional regulatory annotation.
