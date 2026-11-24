# Master’s thesis scientific GO / HOLD decision — Ocular TM DEX (2026-11-02)

**Decision: GO for computational reproducibility/methodological audit master's thesis writing; CONDITIONAL HOLD for formal submission approval and journal-ready original mechanistic article.** This conclusion relies on real public-data computations in GitHub Actions. The thesis is not yet a finished composed graduation document. Academic program eligibility (registration, advisor, coursework and examinations) is a distinct decision to verify with the student's current graduate-school department.

## Exact proposed framing

**Korean:** 덱사메타손 처리 인간 섬유주 전사체의 연구 간 재현성과 세포 분류 및 실험 설계 편향에 대한 계산적 평가

**English:** Cross-study Reproducibility and Robustness of Dexamethasone-Associated Transcriptomic Signatures in Human Trabecular Meshwork

**Central research question:** How reliable are the *directions of transcript associations* in DEX-exposed cultured human trabecular meshwork when the biological unit, cell annotation threshold, cell classification, independent published strain context and treatment-linked library design are explicitly audited?

**Primary endpoint:** methodological robustness and experiment-level identifiability of transcript-associated signatures (NOT drug-specific causal DEX effects, protein mechanism or clinical steroid-induced intraocular-pressure response).

## Evidence fit, true added value, and defensible statements

| Dimension | Verified original result | Scholarly strength | Decision |
|---|---|---|---|
| Primary discovery/source reanalysis | GSE313736: 57 paired donor pseudobulk, 16,050 tested Ensembl genes, 13,717 BH-FDR<0.05 and 83 high-effect 67↑/16↓ internally defined | Reproducible but 85.46% apparent significance is unusually broad and may reflect arm-specific 10x library effects | GO association only, HOLD causal |
| Annotation robustness | identical 57 donors BeamA score 0.75/0.50, 16,050 common genes; identical 55 donors BeamA/scPred, 13,975 common genes; original low-effect <0.25 n=5,349 shows 78.258% sign agreement | Strong original-methodological diagnostic. Not external biological replication | GO |
| Independent published-context comparison | Published Mehrotra 100 nM 48h TM and GSE301517 500 nM 7d TM: MYOC↑, TGFBI↓, MMP2↓ sign concordance; ACTA2 context disagreement | Directional corroboration only; source CPM+1 ratios are different estimands, no 7d gene-wise independent FDR | GO directional, HOLD statistical replication |
| GSE301517 strain sensitivity | 3 strains TM62/TM70/TM72, the 3 candidate signs persist in each; 194 Hallmark EMT: 124/194 all study aggregate signs, only 90/194 also all 3 7d strain signs | Distinguishes mean-signal and strain heterogeneity; correlated genes/post-hoc gene set | GO descriptive |
| Original study library design | 9 source pools × 2 original 10x libraries; 19 columns for original-library FE + DEX, rank 18; DEX adds zero independent rank | Strong formal **identifiability limitation**, not a new drug mechanism or an estimate of how much drug effect is batch | GO audit, HOLD drug-only causal coefficient |
| Pool composition weights | +2.406339 pp pooled BeamA fraction = +0.868638 within source pools +1.537701 pool weight changes, compared with paired donor median −0.421 pp | Valid arithmetic explaining difference between weighted estimands, not causal bias fraction | GO |
| Figures, claims and manuscript references | Four actual final PNG visually checked; original S1 overlap detected and corrected; GitHub Actions original 4 archive source QA 2/2 and manuscript source gate 2/2 | Reproducible assets/technical prose checks; not a composed thesis proof or committee review | GO technical, HOLD submission |
| Original GSE313736 authors' 63↑/17↓ IDs | Genuine 80 gene IDs unavailable in verified source material | Exact overlap with our internally defined 83 impossible to claim | NON-BLOCKING HOLD; DO NOT CONTACT authors |

## Minimum viable thesis structure

1. **Introduction** — steroid-associated glaucoma and TM biology; previously known responses and need to audit computational reliability rather than claim new genes.
2. **Data and study design** — source cohort/technical libraries; explicit **74 reported source donors / 59 original author-QC-related persons / 57 own paired profiles / 55 fixed annotation pairs** labeled as non-equivalent *definitions* (DO NOT assume 59→57 excludes exactly two identifiable donors); include 9 pool × 18 library experimental-design diagram.
3. **Methods** — donor/cell QC, demultiplexing, true integer UMI pseudobulk and donor-blocked edgeR; BH-FDR/effect thresholds; fixed-55 annotation and same-57 score sensitivity; cross-study Ensembl ID matching; 7d three-strain leave-one-strain-out; design-matrix rank; symmetric pool composition formula.
4. **Results** — (i) primary transcript association and high apparent FDR count; (ii) within-source annotation/abundance/size sensitivity; (iii) external signed candidate and 124→90 EMT/strain sensitivity; (iv) rank18 identifiability and cell-fraction weighting discrepancy. Separate inferential from purely descriptive results.
5. **Discussion** — four explicit limitations: original DEX/library confounding; post-treatment label selection; experimental context/non-comparable log2 scales and source strain independence unknown; absent phenotype/protein and author gene80 IDs. Explicitly distinguish methodology contribution from established MYOC/ECM biology.
6. **Conclusion** — transparent, reproducible priorities and *bounds* of interpretation. MYOC↑, TGFBI↓, MMP2↓ are **transcriptional candidates**, not new genetically validated glaucoma drug targets.

## Thesis scientific acceptability gates

**A. Mandatory to progress to full writing (already PASS):**
- Genuine aggregate biological outputs and unique-source provenances; checked real CI.
- Fixed 55 donor annotation sensitivity and source-level gene-effect robustness.
- Independent studies labeled correctly by concentration, time, strain and effect estimand.
- Technical alias rank issue is openly disclosed and cannot be hidden as a routine edgeR batch covariate.
- Figures 1–3 and corrected S1 are actually visually inspected; outputs backed up in verified project Drive folder.

**B. Mandatory before final university submission (pending):**
- Advisor confirms this project is acceptable as a Digital Health master's research contribution, not merely a published-literature review.
- Institution-specific actual thesis examination eligibility and regulations: current registration/earned credits/qualifying exams/assigned adviser, distinct from numerical CI status.
- Assemble *complete* formatted Word/PDF thesis (not just Markdown) containing final text, original 3 main figures, corrected supplemental figure and appropriate separate tables with all source captions and bibliography.
- Manual human verification of every cited scientific claim, final styled bibliography, image/page readability, numbering, discussion completeness and originality/ethics/attribution.
- Use file-backed SHA/CI and actual Drive copies; any applied committee changes require rebuilt, re-reviewed final files.

**C. Not mandatory for this limited master's thesis (optional for stronger paper):**
- Source author's unpublished 80-gene list and researcher outreach, both not feasible.
- Independent new wet-lab proteomics/MMP2 activity or IOP challenge data, unless changing claims to mechanism/clinical prediction.
- Full new raw 10x re-download just to re-audit already certified aggregate numbers.
- New 7d per-gene inferential p/FDR unless manuscript asserts statistical replication (it must not).
- Main PR merge before writing; PR #253 is currently open and may be merged only with ordinary review and user approval.

## Scientific NO-GO triggers

- Thesis title, abstract or claims assert discovery of a novel causal DEX target, a patient glaucoma-risk predictor, proven ECM fibrotic mechanism or drug-dosage/time reversal without appropriate experiments.
- Treat 13,717 BH-significant genes as independent causal drug effects despite 18 source libraries condition-specific preparation.
- Treat 124/194 *aggregate-sign* EMT as equivalent to 90/194 *strict 3-strain-sign* or perform unsupported pathway p-value interpretation.
- Count 7 technical repeats as 7 donors, or claim the three TM strains were authenticated independent patients without evidence.
- Reconstruct the original authors' 80 gene identities from our internally defined 83.
- Equate passing software CI with faculty acceptance, manuscript peer review, or graduation eligibility.

## Final project decision

**Research as a bounded computational master’s thesis: GO.** The verified methodological question and multiple defensible results should be sufficient foundation to *write and submit for advisor review*. Whether the faculty accepts the novelty level as a completed thesis cannot be decided from code alone; formal adviser approval is the key external gate.

**Journal-ready biological/clinical research: HOLD.** A methodological reliability/bias audit might instead be a better fit than claiming novel transcriptional biology, and a new independently balanced culture experiment would be essential for mechanistic causal conclusions.

**One remaining artifact task:** build and proof the final thesis DOCX/PDF with archived figures and manuscript references, checking the user's exact SKKU program template and eligibility before claiming graduation readiness.

**Reference source links (live bibliographic status externally confirmed for context):**
- GEO GSE313736: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736
- Mehrotra 2026 peer-reviewed IOVS: https://pubmed.ncbi.nlm.nih.gov/41636427/
- Source publication/erratum and GSE301517 are separately labeled in the manuscript bibliography.
- Existing proof inventory and Drive backup: [2026-11-02 evidence index](THESIS_PUBLICATION_READINESS_EVIDENCE_INDEX_20261102.md).
