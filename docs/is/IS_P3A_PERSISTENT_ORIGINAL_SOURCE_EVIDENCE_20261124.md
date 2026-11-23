# IS P3A — preserved original-source evidence and weekly independent integrity audit

Date: 2026-11-24 KST. Research: ischemic stroke original P3A chr10 locus, **36 authentic gene-by-tissue QTL assays**, 20 genes, 2 GIGASTROKE ancestry sources EAS/EUR. These are **72 ancestry-assay reporting cells**, not 72 independent causal signals.

## Actual preservation

The four independent previously publisher-source-verified aggregate JSON files were downloaded by GitHub Actions, checked byte-for-byte against their original audited SHA-256, rebuilt against all 36 source assay identities, and committed to the research branch with a manifest INDEX.json. **The originals are now Git tracked** in docs/is/evidence/verified_p3a_20261124/; retention of old Actions artifacts is no longer a prerequisite for checking the archived source claims.

| Original aggregate source report | Original Actions run | Exact JSON SHA-256 |
| --- | ---: | --- |
| QTL GWAS reciprocal FASTA reference and signed-allele SNP reconciliation | 45031834139 | 53b6d71229ed98ab79f404194fa80e94e95bd838fa119ba2961ef234f1308ba8 |
| Original 36-gene QTL permutation number-of-variants vs complete nominal cis | 45036525615 | 4a65a2e4ae9f260c9c154cb8f7a47af96ca95c1492185cd0a9739dbd178f7481 |
| Original 36-gene QTL AN/2 and SE/MAF proxy-N data | 45038801844 | 023750fde4a2a57c1604793dfa2b78ba65a29dbe06ff849e17382a559906eb2a |
| Original independently verified main 72-cell coverage | 45189787019 | 6fc64e2068b44e45a1bf16e1162e20062e4913cca9fc1d3b2e9b702cd118e344 |

This preserves the **original aggregated research findings**, not original per-SNP beta/SE or human genotypes. Original publication source cannot be reconstructed from this bundle alone, but the verified source chain, previous CI execution IDs, and exactly re-checkable 36x2 counts are now durable. Archive privacy checks reject SNP-level effect/donor keys, anomalous source row shapes, unexpected group identities and excessive per-source row arrays.

Original full-cis QTL-tested SNP/INDEL-by-assay entries 202,400 (not 202,400 independent genomic SNPs). Strict original nonpal canonical SNP-only entries 148,344. Direct allele-signed original GIGASTROKE GWAS matching among strict entries EAS 89,597 (60.398%), EUR 128,282 (86.476%). Original EUR ROSMAP allele mismatch 8 gene-assay rows remain HOLD. These are technical coverage diagnostics, not causal gene prioritizations.

## Ongoing automatic enforcement

GitHub workflow .github/workflows/is-p3a-durable-source-archive.yml has:

- One-time branch-only content write. GitHub Actions runs the integrity/privacy and repository regression tests, downloads only these four original preverified privacy-safe source reports and commits untouched JSON bytes and INDEX.json on the approved research branch. It does **not** write directly to main; merge remains subject to PR review.
- Weekly Monday 01:17 UTC / 10:17 KST main read-only validation and a check on pushes to the evidence code or directory. It verifies the **local Git-tracked** aggregate sources, their exact original SHA-256, the complete 36x2 recomputation and the tracked INDEX content (no git differences). It requires **no expiring Actions artifact, no fresh clinical data or FTP download**.
- Failure is fail-closed, not repaired by silently changing source expectations. The whole repository tests were also run during archive creation.

## Scientific distinction: why N/2 and sdY remain pending

Original eQTL Catalogue [association summary columns](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md) use AN/2 as the number of genotyped autosomal samples at the variant. Original [qtlmap usage](https://github.com/eQTL-Catalogue/qtlmap/blob/master/docs/usage.md) combines gene-expression phenotype, sample metadata, covariate adjustment and genotype information. Published genotype allele number does **not** independently establish the actual gene-specific regression sample intersection or the directly observed normalized expression phenotype standard deviation used by a given model.

Known source sample context: GTEx Aorta SNP AN/2 366–387, GTEx Cortex 186–205, ROSMAP original source 560. ROSMAP independently verified 2021 RNA/genotype QC-passing sample count 569, but Release7 nominal 560; no authenticated release-stage reason for the difference of 9 participants. Every actual 36 gene model sdY estimate uses explicit **provisional N**; numeric computation is not evidence that the original phenotype sdY was directly observed or measured.

The original signed beta+SE² case-control GWAS coloc pathway does not require case fraction s, as tested on native pinned coloc R source in merged PR #327. QTL-side N and phenotype scale, correct same-variant full cis input selection, ancestry and prior/LD sensitivity are distinct issues.

**Scientific status remains 0/36 formal coloc.abf and 0/36 coloc.susie approved/executed; no PP.H4, no independent BBJ replication, and no inferred N_eff or sdY=1.**

Source scripts: scripts/is/archive_p3a_verified_source_evidence.py and tests/test_is_p3a_archive_verified_source_evidence.py. The research issue is [#307](https://github.com/kimtk94/IS_Analysis_V3/issues/307).