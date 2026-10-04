---
type: research_master
research: IS
priority: thesis_main
status: phase9d_complete
next_stage: phase10_functional_celltype_atac
updated: 2026-10-05
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
- **Next stage:** **Phase 10 functional cell-type / ATAC validation**
- Phase 11 molecular colocalization remains exploratory until cohort-matched molecular-QTL LD is available.

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

Prioritize macrophage, monocyte, microglia, endothelial and vascular-immune contexts. Bulk-tissue non-confirmation should not be interpreted as rejection of a cell-specific mechanism.

### COL4A2

Prioritize vascular smooth-muscle cells, pericytes, perivascular fibroblasts and endothelial contexts. COL4A1 is retained as a secondary L004 gene rather than merged into a single causal claim.

### ALDH2

Use single-cell expression primarily for biological localization. The main mechanism test is coding/protein/metabolic, not simple bulk-eQTL mediation.

### FGF5

Test vascular and brain regulatory context together with protein-pathway convergence.

## Working thesis structure

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
