# Ocular P9 — TM steroid-response functional assay preregistration checklist

**Status: protocol proposal and public aggregate evidence review, NOT performed wet-lab validation.** 2026-10-31.

## Why the next gate differs from P5–P8

The existing 500 nM BeamA single-cell pseudobulk (57 paired donors, 9 source pools) and 100 nM 48h bulk TM published study yield consistent RNA signals for selected ECM/GC genes. However, the former's treatment arms use condition-specific 10x libraries, and between-study dose differences are confounded by donor cohort, platform, cell-state selection and protocols. Internal leave-one-pool-out edgeR testing, selected-gene concordance, and exploratory pathway enrichment **cannot identify dose-specific causation, protein expression, aqueous outflow change or patient IOP susceptibility**.

## Tiered candidate panel (must be audited from actual source artifacts)

| Tier | Targets | Intended question | Threshold/interpretation |
|---|---|---|---|
| **Core: cross-study concordant RNA controls** | **MYOC, TGFBI, MMP2, ACTA2** | Are glucocorticoid and matrix/cell-state RNA effects reproducible in newly randomized paired cultures? | Prior two-study same sign, both BH-FDR<.05, each abs(log2FC)≥0.3. MYOC anchors glucocorticoid responsiveness, but cannot itself distinguish library batch bias. |
| **Context-specific contrasts** | **SDC1, EGR3, LPAR1, TFPI** | Do previously opposite 500/100 nM effects become genuine dose-specific differences under matched experimental conditions? | Previously opposite-sign and significant in two nonidentical studies. Do **not** infer reversal until tested in one donor/time/batch-balanced randomized design. |
| **Secondary ECM gene priorities** | **LUM, IGFBP3, VCAM1, NTM** | Is the cross-study concordant ECM/adhesion signal broader than original small candidate list? | Both-cohort gene-level FDR<.05, same sign, but post-hoc candidate selection. |
| **Primary-only glucocorticoid assay reference** | **FKBP5** | Does a known GR-responsive target work in the assay? | Primary public aggregate detected it; the P7/P8 frozen selected-gene files do not provide an independently verified 48h numeric value for FKBP5; don't invent it. |
| **Low-abundance exploratory genes** | **CHAC1, ERFE, PARM1** | Is an apparent large log2FC sustained above qPCR detection limit? | logCPM<2 in primary source; separate limit-of-quantitation and primer specificity, not in minimum four-gene clinical mechanism panel. |

**Alias QC**: `CTGF` is the older name for `CCN2`, not a second independently quantified gene. Avoid reporting CTGF plus CCN2 as two validation hits. `FN1`, `COL1A1`, `CCN2`, `THBS1` and `TIMP1` can be negative/ambiguity controls for transcription-vs-ECM physiology; their external gene-wise FDR is not significant in the 48h official source, so they must not be called cross-study validated.

## Minimum prospective wet-lab design

- **Biological unit**: independent human donor TM culture; split each donor into `vehicle`, `100 nM DEX`, `500 nM DEX` at a **common 48h timepoint**. Choose number of donors from a variance/power pilot, not an arbitrary advertised sample size.
- **Randomization**: distribute all three arms from every donor across processing order, culture plates, extraction batches and measurement order; record passage, viability, donor age if permitted, and matched vehicle solvent. Technical wells are not additional donor replicates. If dose remains perfectly aliased with library/preparation batch, pharmacological specificity is unidentifiable.
- **Primary RNA endpoint**: RT-qPCR MYOC/TGFBI/MMP2/ACTA2, stable reference genes empirically checked under DEX, standardized Ct exclusion/assay efficiency, qPCR analyst blind to arm. Predefine donor-blocked dose contrasts, correct multiplicity across core panel. In contrast to existing P7/P8 signs, this is a **new assay yet to be done**.
- **Secondary time test**: 0/100/500 nM × 1h/6h/48h as a feasible second phase, prioritizing contextual SDC1/EGR3/LPAR1/TFPI. Keep time and dose within the **same donors and protocols** before interpreting any direction interaction.
- **Independent mechanism endpoint**: protein/secreted MYOC, secreted MMP2 and gelatinase activity by zymography; TGFBI or ECM deposition/collagen/fibronectin immunostaining; if feasible stiffness, contractility or aqueous outflow surrogate. These are distinct readouts, not computationally validated by RNA-level concordance alone.
- **Contamination/technical QC**: mycoplasma testing, culture sterility/endotoxin assessment as relevant, viability, batch randomization evidence, no-template qPCR and melt-curve/primer specificity, identity controls. Do not infer Mycoplasma negative from 10x 3′ UMI counts.
- **Analysis barriers**: dose×time and donor variation; 100 vs 500 nM dose contrast specified *within one randomized experiment*. Do not label any association a clinically useful steroid responder predictor without patient steroid challenge and ΔIOP outcomes.

## Literature-backed precedent

- [Mehrotra et al., 2026, PubMed 41636427](https://pubmed.ncbi.nlm.nih.gov/41636427/) — 100 nM DEX TM/SCE time course (1h, 6h, 48h); external experimental context.
- [Mao et al., 2017, PubMed 27939447](https://pubmed.ncbi.nlm.nih.gov/27939447/) — cultured human TM DEX treatment, MYOC secretion up; MMP2 abundance/activity down by Western blot/zymography; treatment durations differ from the 48h TM transcriptomic study and cannot be treated as direct protein replication.
- [Delayed glucocorticoid MYOC response, PubMed 11726619](https://pubmed.ncbi.nlm.nih.gov/11726619/) — MYOC as GR-response assay marker by qPCR/immunoblot.
- [ARHGEF26-associated TM myofibroblast/ECM phenotype study (FASEB J., 2024)](https://faseb.onlinelibrary.wiley.com/doi/10.1096/fj.202400400RR) — experimental ECM/stiffness-related phenotype precedent; **not directly established by this project's selected-genome analysis**.

## Research GO/HOLD/STOP checkpoints

| Level | GO | HOLD or STOP |
|---|---|---|
| Computational P9 gene panel | Identical source accession/gene ID and numerical crosswalk; every core agrees in sign and BH-FDR; alias and low-abundance tagged | HOLD if source artifact versions drift or a core marker fails QC |
| Experimental P9 proof of pharmacological effect | Blinded paired donor/condition randomized assay and reproducible RNA effects, independent from 10x library prep | HOLD until actual new donors and dose contrasts tested |
| ECM mechanistic claim | Protein/enzymatic ECM, imaging or outflow functional evidence plus perturbation where feasible | HOLD until functional readouts and biologically motivated perturbation |
| Microbiome add-on | Verified patient/run crosswalk, accessible low-biomass blanks, confounder metadata and independently interpreted shotgun/16S | HOLD currently: ENA BioSample XML patient grouping unresolved |
| Patient steroid glaucoma prediction | Clinical steroid exposure plus paired ΔIOP/glaucoma outcomes, external clinical validation | **STOP for current cell and ocular-microbiome data** |

This document is a **preregistration checklist for future data collection**, not retrospective preregistration of the already post-hoc P5–P8 gene selection or proof of human tissue function. P9a/b/c CI must still pass on immutable artifact references before this plan is marked computationally validated.
