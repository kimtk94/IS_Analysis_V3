---
type: research_master
research: IS
priority: thesis_main
status: phase9f_human_reference_complete_p0_evidence_audit
next_stage: p0_evidence_gates_then_phase10_functional_celltype_atac
updated: 2026-10-10
generated: false
---

# IS MASTER — Ischemic Stroke

> Thesis-main research program. Current working structure: BBJ East Asian ischemic-stroke genetics → molecular convergence → cell-type functional validation.

## Research question

Can East Asian ischemic-stroke loci be connected to reproducible gene-, protein-, and cell-level mechanisms suitable for thesis-grade target prioritization?

## Current phase

- **Phase 9D literature benchmark:** COMPLETE
- **Tissue-label repair / convergence master:** COMPLETE
- **L003 statistical branch:** FROZEN
- **Bulk eQTL stage:** sufficient for prioritization
- **Phase 9F-E human vascular reference:** COMPLETE (2026-10-09), normal temporal-lobe reference localization only, **not disease-state DGE**.
- **Phase 11B coloc/SuSiE:** COMPUTATIONAL_COMPLETE_EXPLORATORY; paper-grade multi-signal colocalization remains **PENDING** validation.
- **Broad discovery V2:** 80 provisional ancestry-aware distance components, 2,425 positional gene–region connections (2,225 unique gene IDs), plus one NEURL1 legacy-anchor-only row.
- **P0 next stage:** Source-denominator, evidence-tier and FGF5 feature QC before costly molecular expansion.

## Evidence integrity update — 2026-10-10

- The legacy 4-locus GTEx ABF stage includes **646 PASS gene–tissue tests**, **43 tested genes** and **16 tissue datasets**; H4 ≥ 0.5 in 6 tests, ≥ 0.75 in 2, ≥ 0.8 in 0. These results are not multiple-testing-adjusted or genome-wide across the expanded universe.
- FGF5 best ABF H4≈0.7787 and SuSiE H4≈0.7522 share input evidence, so they are **not independent replications**. Both remain sensitive to molecular QTL power, priors, ancestry and LD provenance.
- **FGF5 human Phase 9F-E reference:** the processed Seurat RNA counts feature is UNRESOLVED in 17,349 RNA features / 80,515 cells; this cannot establish zero expression. The original ~3.3 GB RDS feature-space and an independent tissue reference require QC.
- **SH3PXD2A cell-type:** human R3_4 dataset has 33,630 Microglia and Macrophages. The gene's top detected cell class was Oligodendrocytes, not macrophages. Neither observation alone demonstrates cell-specific disease mechanism. Donor-level/disease-state tests are pending.
- **COL4A2:** smooth muscle localization is descriptively concordant with the vessel hypothesis, but ubiquitous basement-membrane biology means localization alone does not prove causal regulation.
- **Ancestry scope:** expanded EUR and EAS regions are not directly pooled independent loci; LD-clump/fine-map, effect-allele and source-cohort checks are required.
- **Status and policy:** all historical core genes are *prioritized hypotheses*, not established causal targets. See [IS P0 decision rules](is_p0_decision_rules.md) and the immutable audit under results/is/audits/p0_evidence_20261010_v3 (source results are intentionally not committed).

## Core framework

1. BBJ ischemic-stroke GWAS harmonization
2. Locus definition / ancestry-aware fine-mapping
3. Brain / artery eQTL ABF colocalization
4. SuSiE multi-signal molecular colocalization
5. pQTL / MR integration
6. Functional convergence
7. Human and stroke sc/snRNA + vascular ATAC
8. Final thesis target prioritization

## Primary loci and mechanism branches

| Locus | Core gene | Working mechanism | Key evidence | Current branch |
|---|---|---|---|---|
| BBJ_IS_L001 | **FGF5** | Regulatory / expression / protein convergence | ABF H4 ≈ 0.779; SuSiE H4 ≈ 0.752 | pQTL / MR / BP pathway |
| BBJ_IS_L003 | **ALDH2** | Ancestry-specific coding / protein / metabolic | rs671; BBJ P ≈ 1.99e-18; p.Glu504Lys | coding-protein-metabolic |
| BBJ_IS_L002 | **SH3PXD2A** | Cell-specific immune / vascular regulation | bulk ABF H4 ≈ 0.351; macrophage/monocyte rationale | cell-QTL / ATAC |
| BBJ_IS_L004 | **COL4A2** | Cerebrovascular structural / regulatory | bulk ABF H4 ≈ 0.331; mural/vascular-cell rationale | SMC / pericyte / endothelial ATAC-eQTL |

### Secondary / exploratory molecular screen

- **CALHM2** — ABF H4 ≈ 0.798
- **NEURL1** — ABF H4 ≈ 0.678
- **C4orf22** — ABF H4 ≈ 0.660
- **INA** — ABF H4 ≈ 0.536
- **COL4A1** — secondary gene within L004; weaker bulk support than COL4A2

These remain exploratory until orthogonal validation supports promotion to a core causal mechanism.

## ALDH2 / rs671 freeze

- Primary locus: **BBJ_IS_L003**
- Primary variant: **rs671**
- GRCh37 variant: **12:112241766:G:A**
- BBJ beta ≈ **-0.10975**
- BBJ P ≈ **1.99e-18**
- BBJ EAF ≈ **0.2469**
- JPT conditional GWS variants: **0**
- Low-LD conditional GWS variants: **0**
- Functional consequence: **missense_variant**
- Protein change: **p.Glu504Lys**
- Working interpretation: **rs671-dominant EAS signal with no robust independent secondary GWS signal in JPT-LD sensitivity**
- Paper-grade limitation: **cohort-matched BBJ LD unavailable**

## FGF5 convergence

- BBJ_IS_L001 core gene
- Best bulk ABF H4 ≈ **0.7787**
- Multi-signal SuSiE H4 ≈ **0.7522**
- Working interpretation: **convergent regulatory / expression / protein mechanism**
- Next layer: pQTL/MR and blood-pressure pathway integration

## Cell-type functional priorities

### SH3PXD2A

Compare immune/microglia, oligodendrocyte, fibroblast and vascular cell contexts without assuming macrophage enrichment. Healthy 9F-E contains 33,630 microglia/macrophages with lower SH3PXD2A expression than several other classes. Disease-state donor-level pseudobulk/cell-QTL remains untested. Bulk-tissue non-confirmation alone does not reject a disease-specific mechanism.

### COL4A2

Prioritize vascular smooth-muscle cells, pericytes, perivascular fibroblasts and endothelial contexts. COL4A1 is retained as a secondary L004 gene rather than merged into a single causal claim.

### ALDH2

Use single-cell expression primarily for biological localization. The main mechanism test is coding/protein/metabolic, not simple bulk-eQTL mediation.

### FGF5

First resolve FGF5 row-name/feature coverage in original human RDS. No inference of biological nonexpression can be made from the missing Phase 9F-E counts feature. Test independent cerebellar and vessel references, then pQTL/MR/BP and alternative mechanisms.

## Working thesis structure (hypothesis chapters, not causal claims)

1. **FGF5** — convergent GWAS-expression-protein pathway
2. **ALDH2 / rs671** — ancestry-specific coding/protein/metabolic mechanism
3. **SH3PXD2A** — cell-specific immune/vascular regulatory mechanism
4. **COL4A2** — cerebrovascular mural/stromal regulatory mechanism

## Manuscript figure plan

- **Figure 1** — Main thesis overview: BBJ EAS → loci → molecular convergence → cell-type validation
- **Figure 2** — EAS genetic architecture: four loci, fine-mapping, ALDH2 rs671 branch
- **Figure 3** — Molecular colocalization: FGF5 ABF/SuSiE plus exploratory molecular screen
- **Figure 4** — Mechanism branches: FGF5 / ALDH2 / SH3PXD2A / COL4A2
- **Figure 5** — Cell-type validation: sc/snRNA + vascular ATAC + final target map

## Next actions

- Reproduce audit across 80 provisional components and the 646 ABF tests, with source hashes, separate 2,425 positional pairs plus one anchor-only row, and missingness coverage.
- Register evidentiary statuses and predeclared decision gates from docs/is_p0_decision_rules.md.
- Audit FGF5 original RDS features in a memory-safe environment; validate alternate tissues rather than interpreting UNRESOLVED as zero expression.

- Phase 10 human and disease-state sc/snRNA integration
- Vascular ATAC / regulatory mapping for SH3PXD2A and COL4A2
- FGF5 pQTL/MR and blood-pressure pathway integration
- Maintain ALDH2 as a coding/protein/metabolic branch rather than forcing a bulk-eQTL model
- Document ancestry-matched LD limitations explicitly
- Keep secondary molecular candidates exploratory until orthogonal validation

## MasterOS links

- [[01_RESEARCH/IS/IS_MASTER|IS Master]]
- [[01_RESEARCH/CKD/CKD_MASTER|CKD Master]]
- [[01_RESEARCH/RESEARCH_INDEX|Research Index]]
