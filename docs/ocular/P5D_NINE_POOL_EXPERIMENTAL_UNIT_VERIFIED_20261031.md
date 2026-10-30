# Ocular P5d: completed conservative 9-pool library-level edgeR inference — 2026-10-31

## Actual execution
- **[GitHub Actions #41090827791](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791)**: full PASS including 18 original matrix stream, corrected 57-donor edgeR rerun, conservative 9-pool unit edgeR, no individual donor count exports after deletion.
- [Verified artifact #12995175346](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791/artifacts/12995175346) provides `p5d_poolunit_model_summary.json`, `p5d_poolunit_all_gene_results.tsv.gz`, `p5d_18_pool_library_MDS.png`, `p5d_top100_pool_consistency.tsv`, corrected original donor MDS and full gene-level summary.
- Code: `scripts/ocular/ocular_p5d_poolunit.R`, CLI inputs local runner-private `paired_raw_counts.tsv.gz` and matching donor-to-pool design `samples.tsv`; outputs contain only gene-level statistics and aggregate numbers, no donor identifiers or individual-level expression.

## Model construction
- From the validated 57 paired donor pseudobulks, sum **raw integer counts** within each of the 9 original source pool × Vehicle/DEX conditions, yielding **18 pooled 10x-library-level profiles**. Each original pool has a separately prepared DEX vs Vehicle 10x library; 9 matched library pairs.
- Inference: edgeR TMM / `filterByExpr(group=treatment)` / robust NB QL model `~ pool + treatment`. Rank **10/10**, residual degrees of freedom **8** (vs prior `~ donor + treatment` rank 58/58 and residual df 56).
- This is a **conservative experimental-unit sensitivity** that does not pretend 57 donors are 57 independently prepared 10x libraries. It is not a fully independent replication study or proof that drug treatment and library preparation are separable.

## Actual verified numbers
| Criterion | 9-pool library-unit model |
|---|---:|
| Testable genes after its own expression filtering | **21,081** |
| BH-FDR<0.05 | **13,600** (64.51%) |
| BH-FDR<0.05 and absolute log2FC>=0.25 | **11,652** |
| Common with prior donor-level model (16,050 genes) | **16,049** |
| Pearson gene log2FC vs donor-paired model on common genes | **0.99624485** |
| Spearman log2FC vs prior model | **0.99543167** |
| Original donor-model FDR<0.05 among common genes | **13,716** |
| Genes still significant FDR<0.05 in both models | **11,079**, **80.774%** of prior significant common genes |
| Primary top100 direction retained | **100/100** |
| Primary top100 genes with same direction across all 9 pools | **100/100** |
| Pool unit model common dispersion | **0.02010249** |

Different `filterByExpr` behavior and different experimental units mean 13,600/21,081 cannot be compared as an identical denominator to the previous 13,717/16,050. The strict common-set overlap figures above are the appropriate comparison.

## MDS figure: a significant confounding warning
The pooled-library **MDS dimension 1 visually completely separates Vehicle (left) from DEX (right) for all 9 paired pools**. In the reproduced 57-donor MDS, Vehicle and DEX samples are also visually separated on axis 1. This is compatible with a large true DEX response **and/or** a consistent treatment-associated technical library-batch effect. Nothing in this design can unambiguously assign the separation to one cause because each source pool's two treatment arms are also two distinct 10x libraries.

**Main decision:** technical/statistical robustness of **the within-study differential expression direction** is strengthened. Causal separation of DEX biology from preparation/library artifacts, independent external reproduction, novelty and clinical steroid-induced IOP phenotype remain **HOLD**. Avoid treating 9 source pools as 9 external labs.

## Additional necessary gates
1. Finish the distinct P5b 9 leave-one-pool-out donor-level edgeR refits [run #41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185); its model uses robust=FALSE for runtime while primary is robust=TRUE, so prioritize log2FC direction and correlations.
2. Investigate the 85% donor-model DEG rate and library-wide treatment separation, including a no-effect negative control / independent replicate samples on separate sequencing prep.
3. Repeat full model with broader TM labels and without post-treatment classifier threshold selection. P5e [run #41091146021](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091146021) shows pooled BeamA fraction +2.406pp vs paired donor median −0.421pp, so selection matters.
4. Account for official 74→70→59 donor attrition and 316,799 vs 169,398 singlets [run #41091777826](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091777826).
5. External independent study PRJNA1243143 treatment/time/biological strain annotations remain unavailable from public ENA, so no strict replication yet.
6. No direct IOP data from cultured TM cells: no patient drug adverse-effect predictions.
