# GLGC 2021 EUR (LDL-C / HDL-C / logTG): 22 autosomes technically audited

**Evidence as of 2026-11-16 · technical QC PASS · scientific analysis HOLD.**

The original cohort's complete GRCh37-indexed chromosome-specific public GLGC 2021 EUR summary data were independently read from official BGZF/TBI sources by bounded GitHub Actions jobs for each trait and autosome. This is a technical, **SNP-only** audit. It does not prove source whole-file SHA-256 integrity, durable analysis-ready GWAS files, harmonization, LDSC, MR, colocalization, stroke mechanism or new biological evidence.

## Final independent metadata audit — original evidence, not extrapolated

- **Final GitHub Actions run:** [#43871252114](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43871252114) — **SUCCESS**, complete evidence artifact `ckm-glgc-all22-technical-audit-evidence` (22 per-chromosome aggregates + combined JSON).
- **Coverage:** 22/22 chromosomes PASS; 66/66 chromosome × trait technical checks PASS, zero missing/held chromosomes.
- **Original association rows counted:** **140,353,652**.
- **Retained SNP-only association rows, counted separately per trait:** **128,925,825** (this is NOT the number of distinct loci across traits).
- **Exclusions:** multibase/non-SNP **9,460,687**; symbolic/ambiguous allele **1,964,508**; other/zero-P candidate group **2,632**. Total excluded **11,427,827**, exact row accounting. Do not label excluded alleles biological errors.
- **Normalized canonical SNP duplicate keys:** **0** across individual chromosome × trait shard checks.
- **Ambiguous-AF palindromic source SNPs:** **384,159** records, which must undergo strand/frequency harmonization before MR.
- **Eight code revisions** were used for eight original source batches. The final manifest validates their source commit fingerprints and chromosome count; the audit reports `mixed_code_commits_review_required=true`. A source-identity and regression-parity review remains open.
- **Underflow accounting:** newer aggregate schema reports `underflow_p_zero_rows=2,263`; the initial 12 trait×chromosome aggregates lack that dedicated field. Original chr19 job logs confirm another **369** `P=0` rows; other initial three chromosomes' row accounting did not show other exclusions. Total **2,632** records require original `pvalue_neg_log10` rescue/interpretation before downstream significance filtering. Do not treat `P=0` as no association.
- **Scientific status:** `scientific_status=HOLD`, `whole_source_checksum_verified=false`, `ldsc_ready=NOT_ASSESSED`, `mr_ready=NOT_ASSESSED`, `source_snp_export_scientific_ready=false`. These statuses remain true after technical 22/22 PASS.

## Audited chromosome totals

| GRCh37 chromosome | Retained SNP association rows (3 EUR lipid traits) | Original source aggregate run | Technical |
|---|---:|---|---|
| chr22 | 1,775,635 | [#43431051168](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43431051168) | PASS |
| chr21 | 1,750,553 | [#43848671279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848671279) | PASS |
| chr20 | 2,930,879 | [#43848671279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848671279) | PASS |
| chr19 | 2,947,707 | [#43848671279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43848671279) | PASS |
| chr18 | 3,649,239 | [#43865163287](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43865163287) | PASS |
| chr17 | 3,659,008 | [#43865163287](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43865163287) | PASS |
| chr16 | 4,237,341 | [#43865163287](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43865163287) | PASS |
| chr15 | 3,815,651 | [#43866716532](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43866716532) | PASS |
| chr14 | 4,195,881 | [#43866716532](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43866716532) | PASS |
| chr13 | 4,586,725 | [#43866716532](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43866716532) | PASS |
| chr12 | 6,108,537 | [#43867261217](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43867261217) | PASS |
| chr11 | 6,340,948 | [#43867261217](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43867261217) | PASS |
| chr10 | 6,376,912 | [#43867261217](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43867261217) | PASS |
| chr9 | 5,603,486 | [#43867962230](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43867962230) | PASS |
| chr8 | 7,316,149 | [#43867962230](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43867962230) | PASS |
| chr7 | 7,521,660 | [#43867962230](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43867962230) | PASS |
| chr6 | 8,099,563 | [#43869303834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43869303834) | PASS |
| chr5 | 8,391,709 | [#43869303834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43869303834) | PASS |
| chr4 | 9,120,584 | [#43869303834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43869303834) | PASS |
| chr3 | 9,229,754 | [#43870196218](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43870196218) | PASS |
| chr2 | 11,091,087 | [#43870196218](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43870196218) | PASS |
| chr1 | 10,176,817 | [#43870196218](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43870196218) | PASS |
| **Total 22/22** | **128,925,825** | [Final metadata audit](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43871252114) | **PASS** |

## What is *not* yet established scientifically

1. **Source authenticity and durable rows:** source-side whole-genome compressed-file hashes, stable upstream version manifest and full non-truncated normalized GWAS association data containing source `N`, ALT-effect `EAF`, rsID, beta, SE and `NEG_LOG10P`. Existing jobs retained **QC JSON only**, not all normalized GWAS association records. The separate gated source-SNP export feature passed 11/11 offline unit tests but has not run on real full chromosomes.
2. **P=0 rescue and source scales:** resolve the 2,632 very-small/zero-P candidate records against original `pvalue_neg_log10` and sample size; do not silently discard the strongest associations. Confirm GLGC per-trait transformation, effective sample-size semantics and LD-score assumptions.
3. **Harmonization and LD:** enforce genome-build equivalence; reference/alternate effect sign, palindromic AF-based filtering, allele frequencies and EUR LD reference provenance for LDSC/MAGMA and independent instruments. Check ancestry and sample-overlap of CKM traits, pQTLs and the ischemic stroke outcome.
4. **Remaining CKM exposure sources:** P1 source selection, completeness and annotation for BMI, SBP and T2D (separate from the three GLGC lipid traits) remain distinct gates. Do not equate this lipid sub-component with overall CKM scientific P1 completion.
5. **Stroke biology and causal analysis:** downstream LDSC-SEG/MAGMA, pQTL/eQTL cis-MR, bidirectional/sensitivity tests, colocalization, ischemic stroke subtype definitions and East Asian replication remain NOT_RUN here.

## Evidence handling and next decision

Final CI generated a consolidated report `ckm_glgc_all22_technical_evidence.json` plus 22 source-linked aggregate JSON files in a GitHub Actions artifact. The resulting technical PASS should remain stable and **not** be repeatedly recomputed from scratch. The next effort should prioritize P=0 rescue and a durable, provenance-verified source SNP landing path, with **separate GO/HOLD criteria for biological inference**.

Tracking: [Research Issue #285](https://github.com/kimtk94/IS_Analysis_V3/issues/285) · [Draft PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284) · [Final audit CI](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43871252114).
