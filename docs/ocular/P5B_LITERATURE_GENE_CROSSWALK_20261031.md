# Ocular P5b pre-registered novelty and reproducibility crosswalk — 2026-10-31

## Why the novelty bar is high
The GSE313736 public cultured human TM single-cell eQTL/DEX deposition already reported 74 starting donors and 316,799 total singlets, DEX upregulated 63 and downregulated 17 genes at an abstract-defined high-fold-change threshold. Our reproducible *downloadable final* cell metadata has 59 donors (57 BeamA primary test eligible). Original gene list has **not been verified** from an individual supplementary table; similar total counts alone do not establish exact overlap.

The 2026 Singh et al. multiomic [preprint (PMCID PMC13015409)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015409/) (**not yet peer reviewed** in the cited PMC version) used 3 independently derived primary human TM strains and RNA+ATAC+chromatin/HiC; found **3,423 dex DEGs**, 1,527 higher and 1,894 lower, and 78 IOP/POAG prioritized genes after GWAS/chromatin linking. These are **different technologies, doses/criteria/biological donor counts** from our 57-donor BeamA single-cell pseudobulk, so the raw DEG counts 3,423 and 13,717 cannot be compared as if they used the same analysis universe.

## Exact named gene crosswalk: our P5 primary edgeR vs the Singh preprint text
**Source of our values:** [GitHub Actions P5 corrected primary](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446) [aggregate gene-level table artifact](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446/artifacts/12987425927) `p5_all_gene_results.tsv.gz`. Positive log2FC is induced by DEX.

| Gene | Our log2FC | Our BH-FDR | Earlier study's explicitly described TM DEX pattern | Match? | Novelty verdict |
|---|---:|---:|---|---|---|
| HIPK2 | +2.740051 | 8.731e-54 | +1.7 log2FC; gained super-enhancer at promoter/TSS | YES | **Known** DEX/chromatin candidate |
| FKBP5 | +3.984015 | 3.066e-47 | Increased expression, gained H3K27ac/lost H3K27me3 | YES | **Known** DEX positive control |
| ZBTB16 | +7.513283 | 4.060e-44 | Increased expression; B-to-A compartment and new enhancers | YES | **Known** DEX/chromatin response |
| ATP10A | +1.352570 | 1.147e-25 | Increased expression and enhanced regulatory-region activity | YES | **Known** DEX/POAG locus example |
| EDIL3 | -0.713442 | 3.845e-17 | -0.61 log2FC and super-enhancer loss | YES | **Known** negative-control direction |
| CEBPD | +0.185847 | 1.559e-5 | Listed as a regulatory/topology example; clear numerical expression direction not established from cited text | NOT SCORED | No novelty claim |
| JUN | -0.114965 | 0.1410 | Regulatory/topology example; clear numerical expression direction not established | NOT SCORED | No novelty claim |
| ALOX15B | Not in 16,050 tested genes | N/A | Increased expression/eRNA; do not call absence a discordant DE result | NOT ASSESSABLE | Not tested after filtering |

**Five confirmed examples concordant in signed direction (5/5)**. They are *literature-based positive controls* and not proof of novel target discovery or blinded cohort replication. All 5 were explicitly named in the source paper. These numbers should not be described as exact fold-change concordance because the datasets/model differ.

## Older external DEX RNA-seq/qPCR comparison (2022)
An independent 2022 study of short- and long-term DEX exposure in cultured primary human TM ([full text, PMC9117214](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117214/)) reports Table 2 RNA-seq and qPCR effects for named genes. These are *not* the same experimental dose/time/donor/cell-state setup, so only signed directions are cross-comparable:

| Gene | Our 57-donor BeamA DEX log2FC (48h) | Earlier 2022 RNA-seq 16h | Earlier 2022 RNA-seq 7d | Interpretation |
|---|---:|---:|---:|---|
| **PTK2B** | **+3.14596** | +2.34 | +1.8 | **Already reported** in TM DEX RNAseq and tested by qPCR; do not present as novel discovery |
| FKBP5 | +3.98401 | +2.91 | +4.47 | Known glucocorticoid positive control |
| ZBTB16 | +7.51328 | +7.70 | +6.47 | Known high-effect DEX response |

PTK2B is among our highest-ranked FDR genes, but the direct 2022 evidence substantially reduces novelty of its **expression increase**; any claim of new regulatory mechanism requires a separate experimental or genotype-dependent result. The 2022 paper also identifies FST, MMP1, PDPN, RGCC, SERPINB2, VEGFC; those are literature comparison candidates only, not automatically replicated without checking our table.

## MAOA and mitochondrial fission are also NOT novel DEX/steroid-linked TM mechanisms
A 2026 peer-reviewed experimental study ([PMID 41579974](https://pubmed.ncbi.nlm.nih.gov/41579974/), *Free Radical Biology and Medicine*, doi:10.1016/j.freeradbiomed.2026.01.037) already linked sustained endogenous steroid/cortisol exposure to TM fibrosis, mitochondrial fragmentation and the **MIEF1–MAOA–DRP1** pathway. It reported MAOA upregulation in human TM and gene knockdown/overexpression experiments supporting a functional mechanism.

In our corrected P5 pseudobulk, **MAOA log2FC +4.442772, BH-FDR 7.56e-52**; MIEF1 +0.184497 (FDR 1.16e-10), DNM1L/DRP1 +0.477118 (FDR 3.85e-26). Despite substantial response, **claiming discovery of MAOA-driven steroid TM fibrosis is a STOP** without a demonstrably new cell-state/donor-specific and independently validated mechanism. Cortisol/chronic stress and DEX 48h are **different exposures**, thus effect size magnitude should not be equated.

## Additional independent data
Mehrotra et al., IOVS 2026 ([PMID 41636427](https://pubmed.ncbi.nlm.nih.gov/41636427/)) already found FKBP5, FAM107A and ECM/adhesion DEX effects; 46 ENA runs of PRJNA1243143 do not mean 46 donors, and accession-level titles do not reliably encode treatment/time/strain. External exact directional RNAseq validation remains **HOLD** until sample-specific metadata are reconciled.

## Analysis gates to avoid overclaim
- Formal 9-pool leave-one-pool-out secondary edgeR refits and corrected anonymous MDS are being executed in `.github/workflows/ocular-p5b-pool-sensitivity.yml`. Robust primary and non-robust sensitivity fits must be labeled as different inferential procedures; rely mainly on effect-size/direction stability, not direct FDR count equivalence.
- High discovery fraction **13,717 / 16,050 = 85.46%** warrants technical batch confounding evaluation; DEX and vehicle are distinct 10x libraries in each pool, which donor blocking cannot independently identify.
- Original GSE313736 full gene-set match still pending. The preprint already describes HIPK2, ZBTB16, FKBP5, EDIL3, ATP10A at mechanistic depth, so reporting them as a *novel* main claim is a **STOP**.
- Potential novelty can instead involve **quantitatively validated between-donor steroid-response heterogeneity**, cell-state sensitivity, and independent clinical genotype/outcome-backed mechanism, if data allow. Experimental confirmation remains missing.
- This is *in-vitro* DEX exposure, **not** clinical steroid-induced delta-IOP/side-effect prediction.

Sources:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13015409/ (Singh et al., 2026 preprint, Fig. 2–3)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736 (GEO)
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446 (actual P5 result)
