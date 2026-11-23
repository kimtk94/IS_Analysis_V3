# IS P3A — 36 gene×tissue QTL assays × GIGASTROKE EAS/EUR source coverage matrix

**2026-11-23 KST — research QC only.** Exactly **36 genuine eQTL Catalogue gene×tissue assays** (10 GTEx Aorta, 11 GTEx Cortex, 15 ROSMAP), **20 distinct genes**, **one chr10 G0015 stroke locus**. There are **72 study×ancestry reporting cells**, NOT 72 independent GWAS signals or 72 causal tests.

## Why new reconciliation is needed

Previously verified reciprocal SNP harmonization [Actions #45031834139](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45031834139) reports GRCh38 QTL → GRCh37 GWAS → GRCh38 reciprocal liftOver, original FASTA reference concordance, and signed EAS/EUR allele orientation for **148,344** tested strict canonical nonpalindromic SNP×gene-tissue rows. Prior original nominal QTL versus official original permutation n_variants [Actions #45036525615](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45036525615) verified **202,400** unique within-gene-assay original full cis variant records, including SNPs and INDELs, across all 36 assays (36/36 original permutation tested variant count equals original nominal unique count). The difference **54,056** QTL gene-assay records includes palindromic/noncanonical/filtered SNP and original INDEL source variant classes, NOT all absent GWAS sites. The third actual QTL AN/2 regression model proxy report [Actions #45038801844](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45038801844) independently checked all 36 sample AN/2 distributions.

This new audit combines **the actual three original evidence files**, checks each against a *pinned SHA-256* before touching any group row, and recomputes exactly **36×2=72** gene×QTL tissue×GWAS ancestry summaries and the fully partitioned original source denominator. It does NOT redownload the original 260/283 MB GIGASTROKE GWAS, 3 multi-GB QTL files or UCSC reference FASTA: those were source-verified in the immutable upstream GitHub Actions runs.

The original trusted aggregate JSON SHA-256 and provenance are embedded in the Python audit file and GitHub Actions workflow. The evidence for original source 36 groups remains anchored to source-run IDs and exact unmodified bytes; silent source drift fails the run. Original artifact retention is 30 days; long-term recurring audits must eventually pin the privacy-safe evidence into an immutable repository snapshot or regenerate source-verified evidence, rather than relying forever on expiring Actions artifacts.

## Two genuinely different denominators

- **All original QTL-tested cis variants**: 202,400 gene-assay records = 172,971 SNP + 29,429 INDEL (after within-gene-assay repeated rsID dedup). This is **not** a count of distinct genomic SNPs.
- **Stringent nonpalindromic, single-base, source-valid QTL SNPs**: 148,344 gene-assay records, all with reciprocal GRCh37↔GRCh38 source reference allele matches.
- **Original GIGASTROKE GWAS signed allele match on the stringent QTL subset**: EAS AIS 89,597 records (**60.40% stringent SNP coverage**), EUR AIS 128,282 (**86.48% stringent SNP coverage**).
- **Matched under entire all-variant cis denominator**: EAS 89,597/202,400 (**44.27%**), EUR 128,282/202,400 (**63.38%**). This is **not** "the percentage of all original full-cis variants proven alignable"; it is the observed matched canonical nonpal SNP **share of a larger denominator that also contains untested-by-harmonization variant classes**.
- **Original GWAS missing position among stringent QTL SNPs**: EAS 58,747; EUR 20,054.
- **Conflicting original allele pair among stringent QTL SNPs**: EAS 0; EUR 8 **ROSMAP gene-assay records** across the 8 affected genes NEURL, CALHM2, INA, CNNM2, C10orf32, USMG5, AS3MT, RP11-724N1.1. These 8 gene-assay records may reflect **one shared genomic variant**, not eight independent loci. These are not matched and are not silently flipped.
- **Original ALT vs reference effect sign orientation**: original GWAS effect allele on QTL ALT is directly used; on QTL REF the original GWAS beta must be sign-flipped. Comparison is done after both FASTA and reciprocal liftOver, never by simple same-chromosome position equality.
- The **14,054 EAS** signed-matched records with |original QTL ALT AF − original GWAS ALT EAF| greater than 0.2 are a **nonexclusive descriptive subset of matched SNPs**; do not add to EAS missing counts. Population ancestry and meta-analysis discovery structure make frequency discrepancy a scientific review signal, not immediate allele mismatch proof.

All original indexed source variants, full original-file GWAS bytes/MD5/SHA-256 and historical Release7 QTL identities were verified in upstream source CI; the current job strictly derives and recomputes group-wise aggregate coverage, not a new independent full-source download. Overall "PASS" means the source evidence / arithmetic is intact, not full science-ready SNP membership for a formal Coloc.

## Eligibility and further steps

Current QTL quantitative trait model N and directly observed normalized expression sdY remain **unapproved** for all 36. The model-based estimator under AN/2 hypotheses was numerically calculated for all groups in [PR #326](https://github.com/kimtk94/IS_Analysis_V3/pull/326) and [PR #328](https://github.com/kimtk94/IS_Analysis_V3/pull/328). EAS/EUR GIGASTROKE **chr10 source beta matches ln(OR)** for 725,707 records in [PR #325](https://github.com/kimtk94/IS_Analysis_V3/pull/325); exact original coloc R-source direct beta+varbeta cc path does **not** require case fraction s in [PR #327](https://github.com/kimtk94/IS_Analysis_V3/pull/327).

All source QTL SNP ALT, original GWAS allele, and missingness are audited, but full **palindromic, indel, multiallelic** and QTL gene-specific source cis window/phenotype model/sample context have not been signed off for exact Coloc. A strict common SNP set must be selected before using identical SNP IDs in both model inputs. Dataset QTL and GWAS marginal effects need source-valid effect/SE scale and sensitivity to LD/priors and one-causal-signal assumption. **Formally eligible coloc.abf 0/36; executed 0/36. SuSiE eligible/executed 0/36; no PP.H4 values or independent BBJ replication claims.**

## Files

- Audit: scripts/is/audit_p3a_36x2_original_fullcis_gwas_coverage.py
- Tests: tests/test_is_p3a_36x2_original_fullcis_gwas_coverage.py
- GitHub Actions: .github/workflows/is-p3a-36x2-original-fullcis-coverage.yml
- Source manifest: docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json
- Output: artifacts/is/IS_P3A_36x2_original_cis_allele_coverage_HOLD.json (72 aggregate cells, no original SNP beta/SE/MAF or donor/genotype rows)
