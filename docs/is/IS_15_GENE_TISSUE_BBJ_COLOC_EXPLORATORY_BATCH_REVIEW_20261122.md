# IS legacy 15 assay real-source batch colocalization arithmetic replay — 2026-11-22

**Research scope:** MASTER DEGREE ischemic stroke; Sales prospecting excluded.  
**Technical result:** Original Google Drive archived P01–P15 `stats.tsv` source files (15 separate gene×tissue assays, **11 distinct gene symbols**, **29,661 rows summed across assays**) were actually read in an authorized local analysis environment on 2026-11-22. A second independently implemented Python calculation of the Wakefield single-signal ABF formulas generated **225 model sensitivity scenarios** (15 assays × 3 hypothetical QTL scalar-N multipliers × 5 `p12` values), then source-derived **log-sum sufficient statistics** were saved in GitHub to enable CI mathematical recomputation without exposing variant rows.

**Scientific result:** **EXPLORATORY_HOLD**, NOT a clinical/causal gene validation. This is *not* an R coloc package execution, independent study replication, or a first-time analysis of the missing separate 36 P3A groups.

## Batch results: center assumed model

GWAS `beta/SE` interpreted as case-control log-effect; QTL `beta/SE` interpreted as quantitative trait. Prior effect SDs are 0.2 (GWAS) and 0.15×estimated sdY (QTL); `sdY.est` is calculated from scalar archived `qtl_n`, QTL MAF and var(beta). Default `p1=p2=1e-4`, `p12=1e-5`. `qtl_n` is an **archived assumption**, not an independently source-approved donor effective N. Same-region single-causal-signal assumption remains.

| Panel | Gene | Archived tissue | Source variants | Recalculated PP.H4 | PP.H4 at p12=10⁻⁶ | PP.H4 at p12=10⁻⁴ |
|---|---|---|---:|---:|---:|---:|
| P01 | FGF5 | Brain Cerebellar Hemisphere | 1,639 | 0.779004 | 0.26063 | 0.97241 |
| P02 | C4orf22 | Brain Cerebellar Hemisphere | 1,639 | 0.660664 | 0.16297 | 0.95115 |
| P03 | C4orf22 | Brain Cerebellum | 1,639 | 0.560987 | 0.11331 | 0.92742 |
| P04 | INA | Brain Anterior Cingulate Cortex BA24 | 1,556 | 0.536794 | 0.10385 | 0.92056 |
| P05 | NEURL1 | Brain Cerebellar Hemisphere | 1,557 | 0.677863 | 0.17384 | 0.95463 |
| P06 | CALHM2 | Brain Cerebellar Hemisphere | 1,557 | **0.797936** | 0.28310 | 0.97530 |
| P07 | RPL6 | Brain Cerebellar Hemisphere | 1,990 | 0.211835 | 0.02617 | 0.72883 |
| P08 | TRAFD1 | Artery Aorta | 2,206 | 0.180419 | 0.02154 | 0.68763 |
| P09 | COL4A2 | Brain Putamen | 2,977 | 0.331365 | 0.04722 | 0.83210 |
| P10 | COL4A2 | Brain Amygdala | 2,978 | 0.306811 | 0.04238 | 0.81570 |
| P11 | ALDH2 | Artery Tibial | 1,796 | **0.000057** | ~0.000006 | 0.00057 |
| P12 | ALDH2 | Artery Aorta | 1,796 | **0.000076** | ~0.000008 | 0.00076 |
| P13 | ALDH2 | Brain Cortex | 1,796 | **0.000078** | ~0.000008 | 0.00078 |
| P14 | SH3PXD2A | Artery Tibial | 1,557 | 0.351141 | 0.05134 | 0.84403 |
| P15 | COL4A1 | Brain Amygdala | 2,978 | 0.117435 | 0.01313 | 0.57093 |

**15 assay sum = 29,661 source gene-assay variant observations**, not 29,661 unique genome-wide variants (many assays reuse the same locus). All fifteen source tables satisfy bounded numerical QC; their exact source-byte SHA-256, counts, palindromic variants and large cross-study MAF-gap counts are recorded in [evidence/IS_LEGACY_15_GROUP_COLABF_SUFFICIENT_20261122.json](evidence/IS_LEGACY_15_GROUP_COLABF_SUFFICIENT_20261122.json). Raw SNP-level beta/SE/MAF was **not** uploaded to GitHub.

## Investigated discrepancies and meaningful candidate comparisons

- **CALHM2 (P06)** has the greatest H4 among these fifteen assayed group configurations (0.797936), marginally exceeding **FGF5 (P01)** at 0.779004. The prior `p12=1e-6` lowers them to 0.28310 and 0.26063, respectively. Neither supports an unconditional shared causal signal or a causal gene.
- **ALDH2** in the three *specific* P11–P13 tissues has PP.H4 below 0.0001 and H3 predominance. An earlier study chronicle quoted ALDH2 H4≈0.0790399 as a **historical aggregate/selected result**. It does not represent the same matched tissue model as the specific three re-read inputs unless a source-by-source index is supplied. Do not call the discrepancy a calculation failure or claim the entire ALDH2 locus is excluded: this is model/tissue-specific and requires archived original row provenance.
- The P01/P06/P14/P09 approximations closely reproduce the historical per-assay PP.H4 range (small differences could reflect R float/version, full source assembly or slightly different model settings). However the R coloc v5.2.3 package **was not executed here**; a future direct R rerun using the same source and dependencies remains a separate acceptance criterion.
- Correlated tissues, repeated loci and the same BBJ GWAS cannot be interpreted as **independent replication**. Ranking by max H4 alone is selection-biased; all assayed gene×tissue results and multiple-testing/source overlap must remain visible.

## Technical execution and reproducibility boundary

1. **Local original source stage**: the 15 real Drive files were read from distinct P01–P15 folders. Numeric beta/SE/MAF and archived per-row scalar N were inspected; source SHA256 computed; independent Python Wakefield ABF/logsum implementation ran all 225 scenarios. The full local exploratory report is a sandbox artifact retained separately; GitHub receives *only* three scalar logsum values per QTL-N factor, study-level audit metadata and source hashes, no private raw variant table.
2. **GitHub pure-Python mathematical replay**: `scripts/is/replay_legacy_15_abf_sufficient.py` calculates H0–H4 for each scenario from source-derived sufficient statistics. Three independent CI matrix shards process five groups each, upload outputs and a final job checks 15 distinct panels / 11 gene symbols / 29,661 source observations / 225 posterior scenarios / source hashes and fail-closed scientific status.
3. **Optional future original-source rerun**: `scripts/is/replay_legacy_15_abf_from_original_source.py` can re-read all 15 authorized, staged original `stats.tsv` source files on a machine with access, independently regenerate log-sums and compare them with pinned GitHub sufficient statistics; this is **not executed by public GitHub CI** because private Drive data are absent from the runner.
4. No actual 36-group P3A source-level data was imported. The separate [candidate-wide strategy](IS_CANDIDATE_WIDE_COLOC_STRATEGY_20261122.md) remains 2,225 positional gene IDs / 393 assay groups / 36 potential P3A groups / 25 with basic source stats / **0 independently eligible for fresh P3A**.

Mathematical reference: [official coloc R/claudia.R](https://github.com/chr1swallace/coloc/blob/main/R/claudia.R); in particular `sdY.est`, `approx.bf.estimates`, `combine.abf`. This implementation does **not** validate original study allele orientation, accession/build, donor N, matched LD, cross-cohort independence, multi-signal fine mapping or pleiotropy.

## Immediate next scientific gate

Recover the *original 36-group P3A audit JSON* from its former server or an approved copy, then complete authentic official BBJ/EAS outcome summary source matching, effect SE/case fractions, allele/FASTA build, QTL N and complete cis coverage **for all groups as a batch**. Only sufficiently validated source groups may proceed to a real R coloc/SuSiE Actions matrix. Invalid groups remain explicit HOLD rather than silently disappearing or receiving an inferred posterior.
