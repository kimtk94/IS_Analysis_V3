# HCC Stage29 — original Figure2b effect-vector and biopsy-ID gate protocol

Updated: 2026-10-31. **Code and source provenance work in progress. Scientific gates V03 / V07 remain PARTIAL.** Original Figure2b 5,481 ID source verified separately (#41185931131), original per-gene rho comparison is NOT PASSED until its latest CI completes.

## Hypothesis / estimand separation

1. `Ng 2022 Figure2b` **published** gene-wise Spearman rho effect vector: 5,481 actual Ensembl IDs, original source XLSX SHA256 `16376438aa16a09a13db417a4ae1f4cc16b87ec2c47630aeae651004501be420` (actual `Fig2b=sheet8.xml`, independently verified median rho 0.287059). This is a **published number**, not a newly computed effect.
2. `Suter HCCfull` **computed** gene-wise rho: one-to-one reviewed-UniProt/current mapping 5,343 RNA/protein pairs from 50 biopsy columns/48 patients. Exactly 5,253 Ensembl IDs overlap the published effect list. All 5,253 must be compared by their original IDs, not rank or row position.
3. Separate three panels: **50 biopsies naive**, **48 patient-level mean of repeat biopsies**, **46 singleton-patient sensitivity**. For every panel report `n`, same-gene paired median published/current rho, rho-vector Pearson/Spearman correlation, sign agreement, median signed and absolute delta, delta 10/90 percentiles and fraction |delta|>0.2.
4. The **primary paper-effect comparison** is 48-patient rather than biopsy-naive. Published Ng 2022 processing, normal control, assay batch and original annotation were not recreated; **even high rho-vector correlation is not independent external replication**. Only novel scientific significance if a prespecified independent cohort confirms a priori pathway-level discrepancies.
5. Source clinical patient index is anchored to **original `HCCfull$T` biopsy labels** from official supplement, not an assumed position in other assay matrices. Source `HCCfull$P` can be in a different original column order; require exact equal 50-name **sets**, then explicitly reindex P by T biopsy names before any RNA–protein patient-aware correlation. Fail if names are missing/duplicated.
6. `HCCfullz` may have the same name set but a shared column permutation relative to `HCCfull`; Stage28/28B proved this is sufficient to explain the continuous feature scores, **but cannot yet establish which object's label orientation is the true biological biopsy identity**. Quantify 50×50 cross-object matched patient-group consistency under the official Ng patient index; do not persist a guessed new patient label.
7. `HCCdata.rds` in the same Suter GitHub repository is used by `HCCclutsering.R` with 188 RNA, 116 protein and 158 phosphoprotein input nodes, per source code. It can be used as another source-object consistency check only, **not** independent patient or clinical validation. Only if all 50 sample names and enough common IDs exist should a gene-wise comparison be computed.
8. Separate Stage29 RCA: full RNA 20,356 -> fullz RNA 19,768 excludes exactly 588 samples-constant features. CI #41187760168 independently proves excluded 588/588 gene rows have SD=0, retained 0/19,768 constant (min retained SD=0.008226619). This is source-object transform membership provenance, not separate HCC biological signal.

## Required controls, risk and current execution

- Source immutable commit `cbg-ethz/HCC@438f52271d0a86a6701baa7a96e80c516ba575c2`. Source HCCfull SHA `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`, HCCfullz SHA `9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007`.
- Original Gene namespace integrity: 5,481 published, 5,343 strict current, 5,253 intersection, 228 published-only, 90 current-only. Mixed annotation versions **not assumed equal**.
- Private transient files in GitHub runners: original published 5,481 gene/rho rows; reviewed strict 5,343 gene/protein pairs; official 50-biopsy patient index; **not** uploaded to CI artifacts/logs or committed. Export only aggregate effect statistics, source hashes, and PASS/HOLD gate values.
- Previous #41187890718 failures were implementation errors (R JSON literal typo; `identical` order was imposed where same 50-name sets suffice; publisher original correlation numeric header checking too restrictive). Subsequent corrections are in PR #251; failure is **not** evidence for biology or irreproducibility.
- Last checked corrected multi-job CI: [#41188805830](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41188805830), queued at document time. Updated independent HCCdata author-source audit [#41188805879](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41188805879) queued; independent already succeeded RNA 588-zero-SD check [#41187760168](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41187760168).
- The original paper rho parser now guards source via official publisher SHA, sheet physical relationship, 5,481 ENSG digest and original rho median; publisher's abbreviated Excel correlation header is **not** a robust unique semantic proof on its own. Rho source interpretation stays verified by actual numerical benchmark and published source panel.

## GO / REVIEW / HOLD, no auto-upgrade

**GO_PUBLISHED_GENE_ID_AND_PROCESSING_OBJECT:** official workbook membership and 588 zero-SD omission + same sample permutation + zscore source-object relationship checked.

**REVIEW_PAPER_EFFECT:** wait for actual per-gene concordance CI; if unexpectedly low (<0.5 rho-vector Spearman or <70% sign concordance), examine source gene mapping, assay batch, source normalization and mismatched patient labels before any mechanism claim. These numeric heuristics are **triage thresholds only**, not inferential p-value cutoffs.

**HOLD_BIOPSY_ORIENTATION:** even if full/fullz source-object profiles perfectly align under inferred permutation, original biological labels require an independent original clinical/sample provenance anchor and correct array/metadata direction. Never overwrite source object names based only on maximal profile similarity.

**HOLD_EXTERNAL_REPLICATION_AND_CLINICAL:** actual same-person Gao HBV HCC RNA+proteome sample-matched data must be obtained and independently examined. GSE149614 scRNA is cellular context only, not protein concordance replication.

**Next conditional Stage30:** independent cross-cohort raw/matched matrix manifest and prespecified pathway hypothesis *after* source 5,253 gene-wise comparison and biopsy orientation check. STOP novelty claims if same Ng findings merely repeated without new independent evidence.
