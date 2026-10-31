# Master's Thesis Manuscript Draft — human trabecular meshwork DEX response reproducibility
**Working version:** 2026-11-01, source-evidence-aligned draft (not submitted/published).  
**Proposed title:** *Cross-study reproducibility and cell-annotation sensitivity of glucocorticoid-associated transcriptional changes in human trabecular meshwork*.  
**Important:** The external Mehrotra 2026 journal article has an October 2026 erratum, content not yet checked; any result relying on its archived supplement is labeled **provisional pending erratum comparison**. The Singh 2026 chromatin study is a **preprint** at the version verified.

## Abstract — manuscript-ready *conservative* working draft

**Background.** Glucocorticoid exposure alters gene expression in trabecular meshwork (TM) cells, but differences in cell annotation, experimental design and processing can complicate transcriptomic inference. We examined the reproducibility of TM glucocorticoid-associated transcriptional responses across publicly accessible experiments and evaluated sensitivity to cell-selection decisions.

**Methods.** We reconstructed raw UMI pseudobulks from a 10x single-cell TM dataset (GSE313736), retaining 57 donor-paired dexamethasone (DEX) and vehicle samples for the primary BeamA analysis. Differential expression was analyzed using donor-blocked edgeR quasi-likelihood models. We repeated the analysis at two BeamA confidence thresholds and with an alternative scPred TM_Fibro annotation using the same 55 overlapping eligible donors. Gene-effect directions were compared with published 100 nM, 48-hour bulk TM results and normalized 500 nM, seven-day bulk TM data (GSE301517), accounting for independent TM strain rather than file count.

**Results.** In the original source, 13,717 of 16,050 tested genes met a BH-FDR threshold of 0.05, a very broad response that warrants explicit consideration of batch and selection effects. The 57-donor BeamA confidence reanalysis yielded a gene-level log2 fold-change correlation of 0.999 between thresholds 0.75 and 0.50. Within the fixed 55 donors, the BeamA versus scPred annotation comparison yielded a correlation of 0.950 among 13,975 genes testable in both models. MYOC, TGFBI and MMP2 retained effect direction and statistical significance in both within-source sensitivity analyses. These three transcripts also had the same change direction across three experimental sources, whereas ACTA2 was directionally discordant in the seven-day dataset. Across 194 commonly annotated Hallmark EMT genes, 124 showed the same signed effects in the three settings. The cross-study results are descriptive and the external journal-study contribution remains subject to the published erratum audit.

**Conclusions.** These public-data reanalyses support prioritization of selected glucocorticoid-associated TM transcripts and demonstrate numerical robustness to the tested cell-annotation choices. The results do not establish treatment-independent causality, ECM protein/function changes or an individual steroid-induced glaucoma risk prediction.

## 1. Introduction (evidence-proportional draft)

Elevated intraocular pressure (IOP) is a major modifiable risk factor for glaucomatous optic neuropathy, and the trabecular meshwork regulates a substantial fraction of aqueous humor outflow resistance. Corticosteroid exposure can cause ocular hypertension in susceptible individuals, prompting investigation of glucocorticoid-responsive pathways in human TM cells. Previous experimental work has documented MYOC induction, changes in matrix-related transcripts, cytoskeletal reorganization, and—in some systems—senescence-related phenotypes following DEX exposure.

RNA profiling of primary TM cultures has already been reported at different DEX exposure durations, including a 2022 study assessing shorter and seven-day conditions and a subsequent 2026 peer-reviewed study of 100 nM treatment at 1, 6 and 48 hours in TM and Schlemm's canal endothelial cultures. A separate 2026 preprint integrated DEX-responsive TM RNA profiles with chromatin topology and glaucoma-related genetic association signals. These studies make the detection of glucocorticoid-responsive transcript and ECM pathways an established biological observation rather than a novel discovery by itself.

However, computational comparisons of public datasets can be sensitive to source populations, donor/technical replication, cell-state annotation and treatment-linked sequencing batches. In the single-cell study GSE313736, vehicle and DEX cells were processed as condition-specific 10x libraries, and the chosen BeamA cell label was assigned after treatment. We therefore examined which transcript-effect signs are stable across three different TM experiments and whether a selected signal persists in the original dataset when the cell-confidence threshold or the TM annotation system changes. We present direction and methodological robustness as **candidate-prioritization evidence**, not as evidence that transcriptomic signatures cause glaucoma or clinically predict steroid-induced IOP elevation.

**Aim 1:** quantifiably reconstruct donor-level transcript effects, retaining explicit unit-of-analysis and provenance.  
**Aim 2:** evaluate the sensitivity of these effects to BeamA confidence and independent cell-label algorithms with a fixed shared donor set.  
**Aim 3:** describe signed consistency in two other public TM experimental systems without incorrectly pooling non-comparable effect-size estimators.

## 2. Materials and Methods

### 2.1 Public evidence sources, experimental doses and biological units

The discovery resource was GEO GSE313736, a 10x single-cell RNA-seq dataset from cultured human TM cells exposed to vehicle or 500 nM DEX at the sample-level stated 48-hour condition. Its GEO Series description reports 74 source donors and 316,799 singlet cells, whereas the present strictly reconstructed analysis uses 18 original 10x libraries from nine source pools and **57 donor-paired profiles** passing the fixed criteria defined below. The GEO summary and recovered count-matrix populations are not interchangeable denominators.

The second, independent external transcript source was Mehrotra et al., 2026, an *Investigative Ophthalmology & Visual Science* journal article describing primary TM cultures from nine QC-passed strains treated with vehicle or 100 nM DEX for 48 hours, with earlier 1- and 6-hour experiments. We used the published full-gene summary values rather than refitting that study's raw sample count data. **The article has an October 2026 erratum. Until the corrections and published supplementary table version are compared against the frozen numbers, external gene-level interpretation is provisional.**

The third, separately conducted study was the 2026 **Singh et al. preprint** underlying GEO GSE301517. Its GEO sample metadata indicated 21 RNA samples, comprising seven Resource samples, seven Control samples and seven DEX samples. Treatment/control profiles were obtained from three named human TM strains (TM62: 4+4; TM70: 2+2; TM72: 1+1) exposed to 500 nM DEX for seven days. The seven Resource samples were excluded from treatment-effect comparisons; the technical sample file count was never treated as 14 independent donor observations. A compressed processed gene-by-sample CPM matrix, not raw integer gene counts, formed the third-source transcript comparison.

### 2.2 Primary scRNA metadata inclusion, demultiplexing and cell selection

GSE313736 public 10x filtered matrix files were downloaded for 18 libraries (nine paired pool-level condition sets) with explicit gzip CRC/feature-order/coordinate-record checks. Only confident single-genotype droplet assignments mapped to valid donors were considered. Barcodes were cross-joined between classification metadata and the original 10x filtered count matrices. The primary analysis required at least **30 actually matched selected cells in both the vehicle and DEX arms** for each donor and a donor assignment to one pool. The primary BeamA classification used the deposited CellTypist label `BeamA` with a confidence threshold of at least 0.75, yielding 57 donor pairs. These labels are applied **after treatment**, a limitation intrinsic to the study.

For confidence sensitivity, the same metadata and original UMI matrices were reaggregated at BeamA confidence at least 0.50, retaining exactly the original 57 eligible donor pairs. For independent-annotation sensitivity, the scPred `TM_Fibro` annotation was used, and the analysis was restricted to the **55 donors that passed both selection definitions** with at least 30 matrix-joined cells in both conditions for each annotation. ScPred cell assignment and BeamA assignment can label some of the same cells; this sensitivity analysis does not create a new independent cohort.

The original GEO source sample-level annotation reports 48-hour treatment, but at least one individual GSM growth-protocol narrative contains a conflicting duration reference. This inconsistency should be cited in the source-provenance limitations rather than concealed; it cannot be resolved by the transcript data alone.

### 2.3 Raw integer count pseudobulk aggregation and differential-expression models

Within each selected donor-condition profile, matrix-market `gene × cell` integer UMI counts were summed into a donor-level raw-count vector. Paired vehicle and DEX samples were kept distinct, and source-pool assignments were checked. Biological replication for the single-cell discovery study was thus counted at the donor level, **not at the cell level**.

Gene expression modeling used edgeR negative-binomial quasi-likelihood procedures with `model.matrix(~ pair_id + treatment)`, TMM library normalization and robust dispersion/QL fitting (`estimateDisp(..., robust=TRUE)` and `glmQLFit(..., robust=TRUE)`). Expression filtering used `filterByExpr` with treatment grouping and original library-size information, and treatment-associated expression estimates were evaluated by `glmQLFTest`. Gene-level multiple testing used Benjamini–Hochberg false discovery rate (BH-FDR). The original 0.75-confidence analysis tested 16,050 genes, while 0.50-confidence reanalysis tested 16,728; differences in the testable universe were retained explicitly instead of treating filtered genes as zero or nonsignificant.

The fixed 55-donor dual-annotation comparison used the **same** paired edgeR design for BeamA and TM_Fibro, with 16,051 and 14,003 genes respectively passing expression filters. Gene-wise concordance calculations were restricted to common Ensembl IDs testable in both models, and any selected candidate filtered from one model was marked untestable rather than assigned a pseudo-statistic. No clinical steroid-challenge phenotype was present.

### 2.4 External gene identity mapping, sign consistency and third-study estimand

Genes were mapped by matched Ensembl IDs, with version stripping where documented. For the published 100 nM, 48-hour comparison, we used published gene-level log2 fold-change and adjusted significance values. For the 500 nM, seven-day study, treatment replicates were first averaged **within** each TM strain and condition; a separate descriptive `log2((mean CPM_DEX+1)/(mean CPM_vehicle+1))` was calculated for each of the three strains, and the three strain-specific values were equally averaged. The pseudocount was **one CPM**, and the three strains—not the seven files per arm—were the observational unit for descriptive direction assessment.

We did **not** apply edgeR or DESeq2 to the third study's processed CPM as if it were raw counts, did **not** report a third-study gene-wise FDR or claim a 48-hour versus seven-day causal time contrast, and did **not** pool its pseudocount log-ratio with the 48-hour published model-based log2FC into a numerical effect-size meta-analysis. Directional agreement was reported only for genes with mapped, finite values in the relevant source.

### 2.5 EMT gene set and exploratory outputs

An earlier pathway analysis motivated a post-hoc evaluation of 194 Ensembl-matched genes belonging to the MSigDB Hallmark EMT annotation. We reported the numbers of genes whose effect directions agreed between each pair of experiments and across all three. Since the genes are correlated, this gene-count statistic was interpreted as a **descriptive fraction**, without deriving an independent-gene replication p-value, and **without** equating the EMT annotation with actual increased matrix deposition or reduced aqueous humor outflow.

### 2.6 Reproducibility, provenance and data-access constraints

Original source identifiers, data file hashes, public data handling, parameter settings and final CI run/artifact identifiers were recorded in the GitHub study archive. Key frozen analytical runs include original confidence refit #41199532259, fixed dual-annotation reviewer #41207583296, three-experiment signed comparison #41186736618 and final figure generator #41208972011. Aggregate gene-effect results and anonymized figures were published as CI artifacts. Individual donor × gene raw count matrices, pseudonymous donor linkage files and original genotypes were not uploaded to public artifacts.

## 3. Results

### 3.1 Reconstructed primary donor-unit DEX-associated effects

In the 57 donor-paired BeamA ≥0.75 analysis, 16,050 genes passed filtering, of which **13,717 (85.5%)** met BH-FDR <0.05. Such an unusually large fraction must be evaluated cautiously because the original treatment and vehicle material were processed in different 10x libraries, and library-level artifacts and post-exposure annotation selection may align with condition. Thus FDR measures uncertainty conditional on the fitted model and does not isolate an unconfounded pharmacological effect.

### 3.2 Cross-study transcript directions, with exposure-duration and estimator caveats (Figure 1)

MYOC exhibited positive signed changes and TGFBI and MMP2 negative signed changes in the primary 500 nM, 48-hour BeamA analysis, the published 100 nM, 48-hour bulk study, and the 500 nM, seven-day TM study. In the respective source effect metrics, MYOC estimates were **+4.026, +1.860, +4.349**, TGFBI **−1.869, −0.843, −1.529** and MMP2 **−0.589, −0.357, −0.756**. The third-study figures are strain-averaged, CPM-based pseudocount log-ratios and therefore cannot be interpreted as comparable model-based effect magnitudes. ACTA2 was positive in both 48-hour sources (**+0.469, +0.374**) but negative in the seven-day source (**−0.399**), indicating experimental-context-dependent directional discordance rather than evidence of a causal time reversal.

The primary-source evidence for the 100 nM comparison must be re-evaluated if the published October 2026 erratum changes any source gene-level supplementary result. No equivalent third-study gene-wise inferential statistics were calculated.

### 3.3 Within-source confidence-threshold sensitivity (Figure 2A)

Reaggregation at BeamA confidence ≥0.50 preserved the original 57 paired donor identities but increased the number of matrix-joined selected cells. The original and relaxed filter shared 16,050 testable Ensembl genes and yielded a log2FC Pearson correlation of **0.99911597** and median shared absolute difference of **0.0112462**. MYOC, TGFBI and MMP2 retained effect directions and BH-FDR<0.05. An internally defined large-effect gene screen (`abs(log2FC)≥3`, BH-FDR<0.05) gave 83 versus 88 genes, of which 81 overlapped; the remainder were largely attributable to magnitude thresholds or changed expression-test universes. This **does not constitute an identity-level comparison with the original GEO authors' reported 80 high-effect genes**, whose individual locus list is not yet verified.

### 3.4 Same-donor alternative cell-annotation sensitivity (Figure 2B)

Under the restricted 55 donors eligible in both annotation schemes, BeamA and scPred TM_Fibro yielded 16,051 and 14,003 testable genes; **13,975** Ensembl loci were common. The log2FC Pearson correlation for common genes was **0.94990493**, **12,634/13,975** effect directions agreed, and **10,255** were BH-FDR significant in both models. All three leading candidates retained their direction and statistical significance: MYOC **+4.006/+3.525**, TGFBI **−1.868/−1.941** and MMP2 **−0.591/−0.459** (BeamA/scPred). Of 16 legacy preselected genes, only 15 were testable in both; all 15 were sign concordant and FDR significant, whereas **ERFE was filtered from the TM_Fibro model** and was not assigned a synthetic log2FC/FDR estimate.

Because the alternative annotations reuse donors, libraries and potentially some individual cells, their p-values cannot be regarded as independent replication tests. The high number of jointly significant genes likewise does not eliminate correlated technical biases.

### 3.5 Exploratory EMT direction stability across three experiments (Figure 3)

Among 194 jointly mapped Hallmark EMT genes, **141** effect directions agreed between the two 48-hour sources, **156** between 500 nM 48-hour and 500 nM seven-day sources, and **145** between 100 nM 48-hour and 500 nM seven-day sources. Across all three settings, **124/194 (63.9%)** had the same direction. The set was selected in follow-up analysis and the constituent genes are correlated; these results are descriptive and do not establish pathway-level independent replication or actual TM fibrosis.

## 4. Discussion

### 4.1 Central interpretation and contribution

The data support a narrow but reproducible claim: **glucocorticoid-associated MYOC induction and TGFBI/MMP2 reduction at the transcript level were directionally consistent across three experimental sources** and were robust to two tested cell-classification strategies within the large single-cell source. The principal contribution is not a novel DEX target or previously unknown EMT pathway; rather, it is a transparent cross-study and annotation-sensitivity audit that handles donor-level pseudoreplication, effects of selection thresholds, discrepancies in source expression estimators and incomplete cross-paper provenance.

### 4.2 Why context-discordant genes matter

ACTA2 increased in the two 48-hour sources but showed a negative average ratio in the seven-day third study. Several other candidate genes likewise exhibited different signed responses by source. These differences are not interpretable as direct evidence of temporal or dose-dependent reversal because the studies differ in cell strains, culture preparation, experimental conditions, treatment duration, library methods and statistical estimands. The negative evidence argues against promoting every strongly significant DEG in one study to a generalizable GC-response marker.

### 4.3 Main limitation: DEX and 10x library preparation are not disentangled

Even robust pairing of cell-derived expression profiles by donor is not the same as fully randomized multiplexed experimental preparation. The original 18 10x libraries were assigned by condition, and treatment and processing/libraries are aligned; the high fraction of FDR significant genes raises concern that some apparent magnitude may be technical or related to exposure-driven cell selection. The two cell-classification reruns change annotation but retain the same underlying libraries and donors. Accordingly, the observed numerical stability of MYOC/TGFBI/MMP2 is **not sufficient proof of their causal response in an independent, batch-balanced same-donor experiment**.

### 4.4 Pathway annotation and functional endpoint boundaries

The EMT label refers to a curated gene set and should not automatically be interpreted as a histological epithelial-to-mesenchymal transition in already mesenchymal-like TM cells. A decrease in MMP2 RNA is compatible with hypotheses about reduced matrix remodeling but says nothing definitive about extracellular metalloproteinase abundance, enzymatic activity or outflow facility. TGFBI and MYOC also require orthogonal protein evidence and, for physiological mechanism, targeted functional assays. Earlier publications have already established substantial DEX-related matrix, oxidative-stress and senescence biology, so transcript associations alone cannot be marketed as an entirely new mechanism.

### 4.5 Reproducibility and publication limitations

The source GSE313736 GEO summary reports 63 increases and 17 decreases beyond a stated high-effect threshold, but its deposited Series materials do not expose the complete original authors' list of gene identities. The present project's internal 83 versus 88 high-effect gene analysis therefore does not verify the original 80. In addition, the external journal article has an October 2026 erratum; its substantive correction must be checked against the frozen supplemental gene-level table before final publication. The third independent public study currently available as a preprint used normalized CPM values and three primary TM strains, so the present cross-paper comparison is a descriptive sign check, not a third independent gene-wise statistical test. No linked clinical steroid-induced ∆IOP outcome was available for prediction.

### 4.6 Outlook

The logical wet-lab next experiment would use independent human TM cultures from the **same donors** split into randomized vehicle, 100 nM and 500 nM DEX arms, measured with a common 48-hour time point and independently randomized extraction and assay batches. A 1-hour/6-hour/48-hour factorial could be added after an appropriate pilot sample-size calculation. Pre-specified RT-qPCR of MYOC/TGFBI/MMP2, secreted MYOC measurement, MMP2 protein/zymography and ECM functional phenotyping could test whether stable transcript-associated candidates are relevant to matrix mechanics and aqueous humor outflow. Clinical prediction would require a separate phenotyped steroid-response cohort with actual ∆IOP observations.

## 5. Conclusion — permitted final wording

A donor-level reanalysis of public TM single-cell RNA-seq identified MYOC, TGFBI and MMP2 as transcript-level candidates with concordant effect directions in two other public TM studies. The observed signed association persisted across the tested confidence and cell-annotation sensitivity analyses, while ACTA2 and other genes illustrated the limits of contextual generalization. These results constitute **computational prioritization and robustness evidence**, not causal, proteomic or patient-level glaucoma-risk validation.

## 6. Sources (verify correction before final submission)

1. GEO accession [GSE313736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736), public deposition Sep 2026; specific source-author peer-reviewed journal article for its 80 high-effect genes not verified.
2. [Mehrotra S et al. IOVS 2026;67(2):13](https://pubmed.ncbi.nlm.nih.gov/41636427/), DOI 10.1167/iovs.67.2.13. **Erratum** [IOVS 2026;67(12):34](https://pubmed.ncbi.nlm.nih.gov/42895266/), DOI 10.1167/iovs.67.12.34. Original preprint PMID 41509484 is not the final journal version.
3. [Singh N et al. 2026 bioRxiv preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015409/), DOI 10.64898/2026.03.12.711121, source GEO [GSE301517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517).
4. [Short and long-term effect of dexamethasone on the transcriptome profile of primary human trabecular meshwork cells in vitro (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117214/).
5. [Dexamethasone-induced senescence / TGFβ–Smad3–NOX4 experimental work (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10272683/).
6. [Verified in-repo results, statistical sources and final 3 figure artifact](THESIS_FIGURES_AND_MANUSCRIPT_RESULTS_VERIFIED_20261031.md), frozen [Actions #41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011).

## 7. Explicit manuscript validation gates

| Gate | Status and condition |
|---|---|
| Original 18×10x and 57 donor source provenance | **GO**: integrity/aggregation/model CI artifacts verified |
| Classifier score 0.50/0.75 and fixed shared 55 scPred | **GO**: actual count-based refits, 3 core same signed/FDR |
| Three-source transcript directions | **GO descriptive**, **HOLD** as final definitive published-source evidence until Mehrotra erratum table audit |
| Prove manuscript source's author-reported 63↑/17↓ gene-ID overlap | **HOLD**: exact 80-member list not verified |
| Third-study independent DE p-values | **HOLD**: only normalized CPM and 3-strain direction data analyzed |
| Biological independent protein/activity or functional ECM effects | **HOLD**: not measured |
| Clinical steroid-response prediction | **STOP**: no clinical steroid challenge ∆IOP linked to these TM cell data |
| Claim mechanism beyond existing 2022–2026 literature | **HOLD**: no novel orthogonal functional or causal evidence |

**Do not change these states merely because all programming jobs return SUCCESS.**
