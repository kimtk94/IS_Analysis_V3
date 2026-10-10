# G0022: Source-exact OASIS cis-eQTL signals vs ALDH2 rs671 tagging in EAS/JPT LD

**Analysis timestamp:** October 11, 2026 (KST).
**Result classification:** G0022_TAG_LD_PROXIMITY_ONLY; NO_CONDITIONAL_INDEPENDENCE; NO_VALID_COLOC; NO_CAUSAL_MEDIATOR.

## 1. Why this source-only follow-up is useful

Existing work produced 8,921 GIGASTROKE EAS acute ischemic stroke (AIS) gene–cell–SNP QTL matched records, representing 2,866 distinct lifted GWAS SNPs. The public Japanese OASIS browser (per-gene source Bokeh embedded array) contains non-significant full-looking nominal records. But original donor cis universe and exactly three of four rs671-locus GWAS credible-set variants are not present in the downloaded browser plots. Therefore formal fine-mapping and colocalization are still blocked.

This work asks a smaller, useful question: **are the strongest nominal eQTL markers for ALDH2/BRAP/RPH3A even in the same high-LD block as stroke-associated rs671?** It uses only existing source-controlled 1000 Genomes reference allele counts, no additional inaccessible dataset.

## 2. Audited reference

- 1000 Genomes Phase 3 **EAS n=504**, split into JPT=104, CHB=103, CHS=105, CDX=93, KHV=99 using original panel sample IDs. **JPT is a subset of EAS**, not independent replication.
- PLINK exported 5,060 reference SNPs in reference_export.traw; source REF counted dosages were converted to **ALT dosage as 2-counted** after checking REF/ALT against exact variant ID.
- Only the requested matched source SNPs plus four 95%-CS variants were loaded for this calculation: **2,548 genotype records**. An exact reference allele check exists for every considered variant. For EAS or JPT correlation, pairwise nonmissing samples require n≥80 and nonzero variance.
- The original four credible-set r² checkpoints all reproduced correctly versus rs671:

| Known variant | EAS504 rs671 r² | JPT104 rs671 r² |
|---|---:|---:|
| rs671 | 1.0000 | 1.0000 |
| rs11066015 | **0.9875** | **1.0000** |
| rs11066132 | **0.9742** | **1.0000** |
| rs77768175 | **0.8696** | **0.9204** |

No in-cohort Japanese OASIS donor genotype was accessed. These are public panel LD reference correlations and may differ from the study's own LD.

## 3. Source-exact gene–cell associations and LD

| OASIS gene/cell | GWAS/QTL matched SNPs | Reference-LD-callable SNPs | Strongest marginal browser QTL among all matched SNPs | QTL p | EAS504 r² with rs671 | JPT104 r² |
|---|---:|---:|---|---:|---:|---:|
| ALDH2 activated B L2 | 1,747 | 1,564 | rs2283347 | 0.00616 | NA (reference not callable) | NA |
| ALDH2 Mono L1 | 1,736 | 1,553 | rs73412341 | 9.90e-5 | **0.0340** | **0.0176** |
| BRAP activated B L2 | 1,651 | 1,481 | rs933296 | 0.00891 | **0.0210** | **0.00575** |
| BRAP Mono L1 | 1,632 | 1,465 | rs12425405 | 8.06e-4 | NA (reference not callable) | NA |
| RPH3A Mono L1 | 2,155 | 1,919 | rs7966149 | **6.71e-41** | **0.0894** | **0.1052** |

The QTL lead SNP selection is nominal minimum p among OASIS site-matched regional variants, not a QTL conditional lead or fine-mapped PIP. If a lead is absent from the reference, it remains NA rather than being incorrectly assigned LD=0.

**Major interpretation:** strongest RPH3A monocyte eQTL signal (rs7966149, p~6.71e-41) has **low rs671 EAS r²=0.089**. Its extreme marginal eQTL significance therefore cannot be used as evidence that it shares the stroke/rs671 causal mechanism. A separate QTL signal, linkage, cell power, or different molecular mechanism may be involved. *No conditional independence or multiple-signal colocalization has been computed.*

## 4. Specifically inspect highly tagged rs671 proxy SNPs

Across the G0022 variants that are actually present in both OASIS plots and the 1000G reference panel, there are only **four** matched SNPs per populated gene-cell layer with EAS rs671 r²≥0.8. Three AIS GWAS top CS SNPs (rs671, rs11066132, rs77768175) are absent from the OASIS Bokeh page, even though they are present in reference genotypes; rs11066015 is present (EAS r²~0.988, JPT r²1.0).

| Gene-cell stratum | Number with r²≥0.8 | Of these, QTL p<0.05 | Min nominal QTL p in high-LD subset |
|---|---:|---:|---:|
| ALDH2 Mono L1 | 4 | 4 | **9.13e-4** |
| ALDH2 activated B L2 | 4 | 0 | 0.610 |
| BRAP Mono L1 | 4 | 0 | 0.299 |
| BRAP activated B L2 | 4 | 0 | 0.681 |
| RPH3A Mono L1 | 4 | 4 | **7.62e-7** |

These are not 4 independent SNP instruments or 4 separate causal findings; they are highly correlated proxies. The exact high-LD source set is thin relative to the entire locus and to all four GWAS CS markers. Strongest QTL p at a low-LD SNP and rs671-linked QTL at a weaker p can coexist, plausibly due to **two or more overlapping QTL signals**. Only conditional source QTL fine-mapping could separate this definitively.

No p-value threshold is used to declare no gene expression effect if another cell type shows no QTL association; absence/low significance may reflect low expression, smaller effective N or cell covariates.

## 5. Descriptive figure and reproducibility

Generated SVG **G0022_OASIS_RS671_EAS_LD_VS_QTL_P_DESCRIPTIVE.svg** contains five gene/cell facets of reference EAS r² (X) vs nominal eQTL -log10(p) (Y) for only SNPs callable in EAS/JPT. The five plots contain **7,982 gene-cell-SNP source matched observations**. Gold line marks r²=.8; red dots highlight strongest *among reference-callable SNPs* (may differ from all-browser lead in table). Original source data, not independent candidate generation.

Git scripts:
- scripts/is/audit_g0022_oasis_eas_jpt_tag_ld.py
- scripts/is/render_g0022_oasis_rs671_ld_qtl_evidence_svg.py
- tests/test_is_g0022_oasis_eas_jpt_tag_ld.py
- docs/figures/G0022_OASIS_RS671_EAS_LD_VS_QTL_P_DESCRIPTIVE.svg

Read-only input locations:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/reference_export.traw
- /srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_gwas_fullregional_overlap_v1/G0022_OASIS_EAS_AIS_REGIONAL_ALLELE_MATCHED_SNVS.tsv

New output:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_eas_jpt_rs671_reference_ld_v1/G0022_OASIS_EAS_JPT_RS671_REFERENCE_TAGGING_SUMMARY.json
- G0022_OASIS_GENE_CELL_RS671_LD_REFERENCE_SCORECARD.tsv
- G0022_OASIS_EAS_JPT_RS671_REFERENCE_TAGGING_BY_SNP.tsv
- G0022_OASIS_RS671_EAS_LD_VS_QTL_P_DESCRIPTIVE.svg

Each JSON records original checksum and audit/claim state.

## 6. Next valid, achievable work without difficult DB downloads

1. Search for credible conditional evidence that distinguishes **RPH3A's main molecular eQTL and rs671-associated QTL**. Existing publicly accessible OASIS plot data can support descriptive conditional readiness and LD tagging; a real multivariate conditional test requires *same-cohort LD, allele directions, effective sample size*, so it must remain exploratory if performed on EAS/JPT public reference.
2. Verify which OASIS QTL SNPs are missing due to publication-browser vs original genotype filtering through author-source QC scripts or available browser metadata; do not infer genotype absence merely from not-plotted marker.
3. Report separate molecular hypotheses ALDH2/BRAP/RPH3A and remaining candidates rather than collapsing all to a single causal gene. Retain original broad 2,225-gene IS universe and East Asian/Japanese transferability uncertainty.
4. Ensure the research web displays the new LD-vs-eQTL figure with mandatory limitations. This report is ready to be attached to IS research page; dashboard deployment not asserted in this step.

**Absolute prohibition:** do not turn a small reference r² for a leading molecular SNP into an assertion of independently causal cis-eQTL signal without a conditional test, nor a high r² into shared GWAS/QTL causal variant without coloc.

No original GWAS/QTL source, CKD canonical data, production services, restricted genotypes, or full 39GB OASIS ZIP were modified/downloaded.
