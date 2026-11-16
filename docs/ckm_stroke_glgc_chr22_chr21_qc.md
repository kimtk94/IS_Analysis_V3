# GLGC 2021 EUR: chr22/chr21 partial technical QC evidence

Evidence cut: 2026-11-16 (KST). **Technical data integrity only; not scientific GWAS P1 completion**.

## Confirmed source and protocol

- Public official GLGC 2021 EUR ancestry-specific LDL-C, HDL-C and inverse-normalized logTG summary statistics at `https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/`.
- GRCh37 `POS_b37`; REF is reference allele; ALT is effect allele in this GLGC summary. Columns are verified from the BGZF original file prefix via TLS-protected HTTP 206, and indexed per-chromosome data retrieved using `pysam.TabixFile`.
- SNP-only normalized records have `CHR POS EA NEA BETA SE P`. Multibase alleles, symbolic/ambiguous one-letter IUPAC alleles (e.g. R/D or R/I) are accounted for separately without declaring biological corruption; ambiguous palindromes `A/T`, `C/G` with pooled ALT AF 0.42–0.58 are counted for future frequency-aware harmonization.
- Canonical unordered CHR:POS:allele-pair duplicate check is **full over normalized SNPs within each chromosome**. Wald P vs beta/SE check is **sampled at every 101st variant, with absolute log10(P) difference >1 threshold**, and is NOT a genome-wide statistical validation.
- The complete original GLGC whole-genome file has NOT been downloaded, hashed or normalized; each result is one chromosome only. All LD/ancestry harmonization, cross-study sample-overlap investigation and MR/colocalization remain NOT_RUN.

## chr22 verified full chromosome — [Actions #43431051168](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43431051168)

Run source commit: `91de63c4916a4f078ca5490b35c530b1b1709b31`.
All five jobs SUCCESS: parser/contract test, three full-chromosome SNP QC jobs, and the chr22 aggregate. Aggregate `ckm-glgc-chr-aggregate` status PASS (technical), `scientific_status=HOLD`.

| EUR trait | Original rows | Retained SNPs | Multibase/non-SNP | Symbolic or ambiguous allele | Duplicate keys | Sampled P mismatches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LDL-C | 648,183 | 590,273 | 43,832 | 14,078 | 0 | 0/5,844 |
| HDL-C | 636,267 | 592,304 | 43,875 | 88 | 0 | 0/5,864 |
| logTG | 651,154 | 593,058 | 44,019 | 14,077 | 0 | 0/5,871 |
| **Total** | **1,935,604** | **1,775,635** | **131,726** | **28,243** | **0** | **0/17,579** |

## chr21 verified per-trait full chromosome — [Actions #43848671279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848671279)

Run source commit: `fea2ce570e5abcad7f1dc80e5287e8639baed514`.
Three actual chr21 source jobs independently completed SUCCESS with each source JSON `status=PASS`, `scientific_status=HOLD`, no row-cap truncation.
The **chr21 aggregate has NOT yet been generated**, because the rest of the chr21·20·19 batch was still queued as of this evidence cut.

| EUR trait | Original rows | Retained SNPs | Multibase/non-SNP | Symbolic or ambiguous allele | Duplicate keys | Sampled P mismatches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LDL-C | 640,388 | 581,991 | 44,910 | 13,487 | 0 | 0/5,762 |
| HDL-C | 629,325 | 584,370 | 44,955 | 0 | 0 | 0/5,785 |
| logTG | 642,765 | 584,192 | 45,087 | 13,486 | 0 | 0/5,784 |
| **Total** | **1,912,478** | **1,750,553** | **134,952** | **26,973** | **0** | **0/17,331** |

## Expansion completion gates

1. First batch run #43848671279 must finish all **9/9 chr21/chr20/chr19 × three traits** real-data jobs with technical PASS.
2. Independent `ckm-glgc-batch-summary-21`, `-20`, `-19` aggregates must each be `status=PASS`, no missing trait, no row-accounting error, and same committed code version across the original first batch.
3. [Sequential batch validator CI #43850181958](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43850181958) completed 12/12 offline tests PASS. The actual next batch chr18/chr17/chr16 remains NOT_DISPATCHED until conditions 1 and 2.
4. Do not equate all-autosome SNP technical QC with full-source integrity or downstream genetic association, MR or LDSC. Summary statistics `N` and `POOLED_ALT_AF` must be preserved or independently available for rigorous downstream harmonization.

Tracking: [Issue #285](https://github.com/kimtk94/IS_Analysis_V3/issues/285) · [Draft PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284).
