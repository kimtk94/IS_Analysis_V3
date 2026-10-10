# IS 30 GTEx QTL assays — existing EUR LD panels GRCh37→GRCh38 audit

**Date:** 2026-10-11 KST. **Read-only, existing data only.** No novel causal result.

## Source and build

The 54 1000 Genomes Phase 3 EUR reference regional PGEN/PVAR panels in the metabolic-resilience workspace have PVAR metadata documenting **GRCh37/b37/hs37d5**, not GRCh38. Original GTEx QTL variant positions use **GRCh38**. Therefore an earlier direct position-REF-ALT join yielding 0 SNPs was build-incompatible and was not adequate to conclude that EUR panels lack all IS SNPs.

Using the existing UCSC official hg19ToHg38 chain and a chromosome-specific best-scoring primary forward chain, the EUR panels were lifted to GRCh38 with exact displayed REF/ALT string retention (no reverse-complement transformation). We compared original gene–tissue SNP inputs for **30 numerically replayed ABF assays**, distinct from the untested original 616.

## Observed result

- Existing EUR PGEN/PVAR panels with readable variant tables: **54**.
- ABF assays with a directly allele-exact lifted SNP match in these 54 panels: **0/30**.
- ABF assays with at least 50% SNP coverage: **0/30**.
- Max per-assay matched SNPs across available region panels: **0**.
- EUR regional panel chromosomes: 1,2,4,5,7,9,11,15,16,17,20,22.
- 30 tested BBJ/GTEx assay chromosomes: 4 (16 assays), 10 (10), 12 (1), 13 (3).
- Only chromosome 4 occurs in both sets; existing EUR panels were developed for different regional loci.

No EUR genotype donor sample is asserted to be a GTEx study participant. Lifted coordinate/allele matching does not by itself establish reference FASTA, allele strand correctness, or an in-study signed LD matrix.

**Interpretation:** The existing *regional EUR cache* cannot currently supply useful shared-SNP LD for these IS assays. This is a scope limitation of cached regional loci, not proof that 1000 Genomes EUR reference genomes or public GTEx LD do not exist. Full chromosome extraction or region-focused EUR data provisioning remains a separate workload.

## Reproducibility

Script: `scripts/is/audit_is_eur_b37_liftover_coverage.py`.

Audited results: `/srv/is-analysis/results/is/audits/is_30_eur_b37_lifted_coverage_20261011_v1/`
- `IS_30_EUR_REFERENCE_B37_TO_B38_LIFTED_SNP_OVERLAP.tsv`
- `IS_30_EUR_B37_LIFTOVER_AUDIT.json` (input chain SHA256; 30-assay source count)

Scientific gate remains **NO_GTEX_IN_STUDY_LD, NO_VALIDATED_MULTISIGNAL_COLOC, NO_CAUSAL_GENE**. Preserve full 2,225-gene discovery scope. Do not treat untested status as lack of biological association.
