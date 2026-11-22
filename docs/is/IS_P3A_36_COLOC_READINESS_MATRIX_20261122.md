# IS P3A — 36 gene×QTL group Coloc input readiness audit (2026-11-22)

**Purpose:** MASTER DEGREE / ischemic stroke. Read-only scientific input quality control. No Sales work, no MR or new H4 inference.

## Observed CI result

- [First complete source-joined execution #45023378225](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45023378225): **SUCCESS** after updating the QTL source recovery to the newest verified `main` artifact.
- Pinned input artifacts: recovered QTL group data [#45021251945](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45021251945), real EUR/EAS publisher chr10 scan [#45020605713](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45020605713), independent signed source reproducibility of 24 EAS SNPs [#45022609773](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45022609773).
- Synthetic fail-closed regression tests: 10/10 PASS.
- 36 authentic gene×assay keys; 20 gene IDs; one GWAS region `IS_XDATA_G0015` (chr10).
- QTL studies: QTD000131 10, QTD000171 11, QTD000434 15 groups.
- Official BH05 public eGene source matches: 25/36, with no BH05 row for the other 11. This is **not** equivalent to proving the previously reported 25/11 individual server source-stat membership.
- Source-only assay×variant row sums across those 25 official BH05 eGenes: 142,196 cis tested, 59,090 BBJ/EAS nonpal overlapping, 31,313 QTL–EAS AF delta ≤0.1. Counts are across group records, **not distinct SNPs**.
- Full verified publisher context: EUR 4,397 valid chr10 variant rows; EAS 3,913. These are common region counts, not groupwise joined SNP/SE inputs.
- Cross-region original 24/24 EAS nonpal SNP reproductions cover 16 chromosomes, but the 36 P3A QTL groups share chr10 where the archived region lead is G/C palindromic.

**Final scientific status: HOLD. New full P3A coloc.abf eligible/executed 0/36; coloc.susie eligible/executed 0/36.**

## Required scientific inputs for future eligibility

Each of the 36 natural keys (`region_id, QTL_dataset, gene_id`) must have original cohort-specific **full cis** association records, with variant ID, GRCh37 REF/ALT checked against a specified FASTA reference, effect/other allele, signed beta, measured SE, AF/MAF, SNP-level QC and source provenance. The full shared QTL/GWAS variant set must be deterministically established without selecting variants based only on GWAS significance.

- **Molecular QTL:** genuine tested cis SNPs, approved tissue/study donor N, cohort ancestry, effect/SE interpretation, variant-level source checks and per-assay source-stat membership.
- **GWAS:** official GIGASTROKE/BBJ accession-specific effect and SE for matched SNPs, exact AIS phenotype, population/ancestry, appropriate Ncase/Ncontrol and fraction (not aggregate totals from unrelated EAS cohorts), independence limitations.
- **Standard `coloc.abf`:** needs source-verified beta/SE (or independently validated alternative method), MAF and sample size/model information, full region overlap and priors. **No LD matrix is required by `coloc.abf` itself.** The single-causal-variant assumption is an explicit limitation.
- **`coloc.susie`:** additionally needs matched, SNP-order-consistent, allele-aligned signed LD with ancestry-appropriate fine-mapping QC, multiple-signal diagnostics and SuSiE settings.
- Bayesian prior sensitivity (p1/p2/p12), H3/H4 and warnings must accompany any future reported result. A valid H4 does not by itself prove the gene is causal.

A missing scientific input must remain HOLD; never infer N_eff, fabricate sample size, compute SE from nominal p without a validated model, or backfill cohort identity from unrelated GWAS metadata.

## Automation files

`scripts/is/audit_p3a_36_coloc_input_readiness.py` writes all 36 explicit statuses to a JSON + TSV and a concise markdown. `tests/test_is_p3a_36_coloc_input_readiness.py` includes ten failure-case controls. `.github/workflows/is-p3a-36-coloc-input-readiness.yml` retrieves fixed genuine CI evidence and runs the read-only join; no patient rows, per-SNP original data, credentials, or LD are uploaded.

The current gate is **a verified readiness matrix, not a full-cis variant harmonization algorithm**. Its source input deliberately contains only grouped original validation artifacts. Once genuine SNP-level and phenotype-level data are accessible, add a separate per-variant normalization/cis audit and run it against this 36-key registry; only its independently validated outputs may be candidates for a new ABF executor.

## Next research work

Obtain or re-access exact QTL Catalog source nominal full-cis files for QTD000131/171/434 (or the original legitimate archived source), their cohort donor N and allele reference, and exact GWAS cohort case/control counts. Prioritize recovering those *sources* rather than issuing gene-by-gene H4 estimates. Preserve all 36 keys and explicit blockers in subsequent CI runs.
