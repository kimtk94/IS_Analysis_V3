# G0022 — Japanese OASIS sc-eQTL full-nominal source feasibility

**Date:** 2026-10-10 KST
**Research status:** NOMINAL_SOURCE_PREFIX_VERIFIED; G0022_FOUR_CS_VARIANT_COVERAGE_NOT_ASSESSED; GWAS_QTL_COLOC_BLOCKED.

## 1. Original official source

The full-locus GIGASTROKE East Asian AIS SuSiE result includes four highly correlated variants (rs11066015 / rs671 / rs11066132 / rs77768175). Recent Japanese ImmuNexUT and JCTF sources offer ALDH2, BRAP and RPH3A *marginal* molecular evidence. Existing JCTF v2/v3 public SNP-gene QTL rows are threshold filtered, and cannot serve as unbiased all-tested cis source data for a valid shared causal-variant analysis.

OASIS E-GEAD-1054 is a **different Japanese single-cell immune QTL resource** (not ImmuNexUT E-GEAD-420). The official DDBJ IDF/SDRF states: 40 immune cell types, tensorQTL, GRCh38, 234 COVID-19 patients and healthy controls. The research paper reports 235 donors, a cohort/source-count difference needing reconciliation. Neither 234 nor 235 is verified per-cell-type or per-SNP genetic effective N.

Official index: https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-1000/E-GEAD-1054/
Study: https://www.nature.com/articles/s41588-025-02266-3
JOB browser: https://japan-omics.jp/variant?input_value=rs671

## 2. Actual nested ZIP/TAR/GZIP verification

We issued one HTTP 200 HEAD and a **bounded HTTP 206 Range GET for 262,144 bytes** of the original public DDBJ archive. The response was partially decompressed (without creating the full archive) using Python ZIP raw DEFLATE, TAR header and nested GZIP parsing.

| Attribute | Source-verified value |
|---|---|
| DDBJ E-GEAD-1054.processed.zip exact size | **39,271,327,981 bytes** |
| Official MD5 ZIP (NOT full-download re-checked) | 312f19beade7c98afdd19d8d3965c270 |
| ZIP first/only member listed in manifest | eQTL_summary_statistics.tar |
| Outer ZIP method | 8 (DEFLATE) |
| Official TAR extracted byte length | 39,259,422,720 |
| TAR first entry | B_Activated_PC15_MAF0.05.cis_nominal.txt.gz |
| First nested gzip compressed size | **970,351,338 bytes** |
| Downloaded archive body | **262,144 bytes** |
| Full original ZIP downloaded? | **No** |

ZIP's central file directory cannot directly index specific genes or cell types inside a single huge DEFLATE-compressed TAR; no convenient targeted HTTP byte range extraction of G0022 is supported by this packaging. A separately indexed file or staged full extraction is needed.

## 3. Actual original test records: non-significant variants INCLUDED

The source first member contains the exact columns:

    phenotype_id, gene, variant_id, tss_distance, pval_nominal,
    slope, slope_se, af, ma_samples, ma_count

The first **80 complete variant–gene rows** all had p>0.05 (first 5 were 0.3488858046, 0.7161006781, 0.5079023450, 0.5079023450, 0.6135907861). Original SNP ids use GRCh38 style chr1_55326_T_C; first observed gene is LINC01409 (chr1). This **proves this B_Activated source file includes unfiltered nominal, non-significant SNP-gene tests**, unlike JCTF released thresholded subsets. Inference about remaining cell-type files is not yet supported.

The original QTL variant identifier appears to include REF/ALT, but the beta effect-allele direction was not independently checked against official author documentation. The filename contains MAF0.05 filtering. Individual cell types may have different sample sizes and variant retention.

Crucial: The present audit **has not assessed whether rs671 or any of four G0022 credible-set SNPs is present** in E-GEAD-1054. This is NOT_ASSESSED (machine-readable null), not zero actual biological associations. No ALDH2/BRAP/RPH3A full cis tested locus extracted, and no coloc or mediation result computed.

## 4. Interpretation and priority

1. OASIS full-nominal Japanese eQTL statistics have better potential for valid EAS GWAS coloc than the JCTF threshold-filtered releases, but G0022-specific SNP and gene-assay coverage has not yet passed.
2. Highest-value next input is an authorized targeted export of full original chr12 G0022 cis associations (including non-significant variants) for ALDH2, BRAP, RPH3A and distal gene competitors, across tested cell types. Prefer per-cell indexed BGZF/tabix or Parquet.
3. JOB OASIS browser has gene and cell-type result pages, but browser results must not be treated as a complete full-cis denominator. Public rs671 variant page showed empty OASIS HTML table bodies in server-rendered HTML; this does NOT prove biological null, since dynamic loading/selection may exist.
4. If targeted source export is unavailable, a controlled resumable full ~39GB acquisition can be separately planned with MD5 validation. **It was not initiated during this task.**
5. When complete QTL source exists, check both presence and absence of ALL four high-LD credible-set SNPs and all tested comparison SNPs; check exact GRCh38-to-GRCh37 allele map, per-cell N/MAF, effects and SE, matched Japanese signed LD, and GWAS sample/cohort reuse. Then perform proper signal-level coloc and sensitivity checks.
6. The current locus is **ALDH2/BRAP/RPH3A multi-gene hypothesis**, not a proven stroke molecular mediator, pharmacological target, or causal SNP. Keep the broad 2,225-gene IS context.

## 5. Reproducibility

Code:
- scripts/is/audit_g0022_egead1054_nominal_archive.py
- tests/test_is_g0022_egead1054_nominal_archive.py

Real output path:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_egead1054_nominal_v1/G0022_EGEAD1054_SOURCE_NOMINAL_PROBE.json

Run:

    python3 scripts/is/audit_g0022_egead1054_nominal_archive.py --probe-remote --outdir /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_egead1054_nominal_v1

Verified source log gates:
- all original first-80 p>0.05
- HTTP Range enforced, body 262144 bytes
- G0022 four-variant coverage NOT_ASSESSED, null
- no valid new full-locus coloc; no causal target

No canonical GWAS, previous QTL, other projects or server production services modified.
