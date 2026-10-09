# IS GWS expanded EAS LD — V2 real-data QC (2026-10-09)

## Provenance and isolation

This V2 analysis adds expanded 1KG EAS GRCh37 reference panels for two originally BBJ-limited groups; the original V1 results remain unchanged.

- G0015 / chr10: `10:104163676-106220351`, 51,737 biallelic SNPs, 504 EAS samples
- G0022 / chr12: `12:109860321-113909176`, 101,913 biallelic SNPs, 504 EAS samples

The expanded panels are new `*.stable.{pgen,pvar,psam}` triplets. The pipeline selects expanded triplets when **complete**, otherwise retains legacy reference. A partial triplet raises an error.

Source results: `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/`.

V2 output isolated at:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/`.

## V1 vs V2 (real GWAS data)

Comparing the same 6 GWAS × 7 GWS groups (42 study-group pairs):

| Reference metric | V1 | V2 | Change |
| --- | ---: | ---: | ---: |
| GWAS/window variant records | 145,949 | 145,949 | 0 |
| Nonpalindromic allele oriented | 107,697 | **133,509** | **+25,812** |
| Palindromic requiring review | 1,465 | **1,760** | +295 identified, still not oriented |
| Missing reference positions | 36,787 | **10,680** | **−26,107** |

By expanded group:

| Group | V1 oriented | V2 oriented | Gain | V1 missing position | V2 missing position |
| --- | ---: | ---: | ---: | ---: | ---: |
| G0015 chr10 | 13,541 | **25,972** | +12,431 | 15,690 | 3,068 |
| G0022 chr12 | 21,661 | **35,042** | +13,381 | 15,557 | 2,072 |

**Unchanged five GWS groups:** V1 and V2 oriented counts are identical, by design. The comparison script rejects regression in expanded groups or unexplained changes in unchanged regions.

V2 orientation success: 133,509 / 145,949 = ~91.5%, compared with V1 ~73.8%. This is *record-level allele QC*, not 133,509 independent variants.

V2: 42/42 lead alleles oriented, 42/42 ALT-EAF differences versus EAS 1KG below or equal 0.10. All 42 pairs remain blocked from genuine multi-signal fine-mapping until LD matrix consistency, sumstat variant overlap, replication source and credible sets are checked.

## Reproducible commands

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
V2=/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2
python3 scripts/is/audit_gws_gwas_reference_overlap.py --root "$V2"
python3 scripts/is/compare_expanded_gws_ld_coverage.py
python3 scripts/is/build_gws_study_readiness.py --root "$V2"
python3 scripts/is/audit_gws_lead_eaf.py --root "$V2"
```

## Cross-ancestry next step

GCST90104540 EUR AIS full raw GRCh37 file must match the original GWAS Catalog MD5 before normalization. Do not treat an unfinished `.part` or synthetic test result as a real phenotype result.

EUR discovery is separate from the 30 EAS/Japanese coordinate components; it does not imply ancestry-independent causal loci and requires EUR ancestry-matched LD for confident fine-mapping.

## Important caution

- GWAS from GIGASTROKE AS and subtype files can share participants. 42 source-group comparisons are not 42 independent replications.
- All 869 positional gene-region rows and original nine anchors remain available; marker density improvements do not prove effects on any gene.
- Palindromic SNPs are kept as `PALINDROMIC_REVIEW`; never silently harmonize ambiguous strand orientation.
