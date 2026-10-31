# HCC Stage24 — Prespecified processed RNA–Protein correlation audit

**Status:** engineering pilot, not original Ng 2022 reproduction or new biological discovery.
**Created:** 2026-10-31.
**Implementation:** research feasibility PR #251, GitHub Actions #41173073247.
**Evidence/decision governance:** PR #252; keep Draft until reviewed.

## Research question and denominators

Verify whether the complete Ng/Suter processed multimodal matrices support technically interpretable **within-gene / across-patient** RNA–Protein rank correlations, after conservative molecular ID alignment. This is *not* an RNA DE analysis, not etiology inference, and not an incremental clinical prediction model.

- Immutable public Suter source `cbg-ethz/HCC` commit `438f52271d0a86a6701baa7a96e80c516ba575c2`; `HCCfull.rds` SHA-256 `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`.
- Full feature axes: RNA Ensembl 20,356 × 50 tumor biopsies; Proteome 5,585 × 50.
- Stage23 current reviewed UniProt snapshot: 5,390 mapped proteome features; 5,437 distinct gene *union* including ambiguity; **5,343 strict one-feature ↔ one-gene pairs**, excluding 43 multi-gene proteins and duplicate singly-mapped proteins. Exact paper 5,490/5,481 identities unverified.
- Ng original source Supplementary Data 1 SHA-256 `edf1d2afd3baaedd56edee8536c16f10e4a62cb6cdac929863a062671ee73d25`; 50 exact biopsy labels map to 48 distinct patients. There are two extra tumor biopsies. Sample labels, patient identifiers and gene–protein row mappings are runner-temporary and never uploaded.

## Prespecified analysis

1. Retrieve immutable `HCCfull.rds` and verify hash, numeric feature axes and exact RNA/Protein sample order. This Suter RDS is already normalized and potentially imputed/batch corrected; it is not Ng raw RSEM expected counts.
2. Build **current** UniProt reviewed Ensembl GeneId crosswalk; require strict 1:1 gene–protein feature mapping, deduplicate source gene targets and exclude multi-gene assignments. Retain separate batch digest and annotation timestamp.
3. For each mapped gene, calculate Spearman correlation between its RNA and protein tumor values **across** the 48 independent source patients. Primary descriptive version averages the two biopsies of patients with duplicate samples *on the processed value scale* (an explicit sensitivity assumption).
4. Repeat with only singleton-biopsy patients (expected 46 when two patients each have two biopsies) and with all 50 biopsies as a **nonindependent sensitivity diagnostic**, never as N=50 independent patients.
5. Export only anonymized cohort-wide aggregates: number of finite/undefined per-gene correlations, median, q10/q90, fraction positive and fraction |rho|>=0.5 for each scenario, plus median change when duplicate patients are excluded. Raw patient, biopsy, gene and protein identifiers, individual expression and per-gene `rho` or p-values are **not** published.
6. Require Stage23 provenance validation, `>=5,000` finite mapped gene pairs, 50/48 denominators, consistent source hashes and no inferred clinical conclusions before the technical workflow can PASS. Median change >0.05 is a **sensitivity review flag** to interpret, not a statistically validated cutoff.

## Original paper benchmark and limits

Ng et al. (2022), Nature Communications 13:2436 (doi:10.1038/s41467-022-29960-8), used 5,490 comparable entries for dysregulation quadrant analysis and **5,481** for mRNA/protein Spearman tests; its reported **median mRNA–protein Spearman rho=0.287** was computed on original tumor-vs-normal relative/normalized measurement definitions. These numbers are **not** a target threshold or a proof of identical gene membership for this processed Suter reanalysis.

Suter et al. (2022), PLoS Computational Biology (doi:10.1371/journal.pcbi.1009767), explicitly report additional preprocessing/imputation/batch correction, sample exclusion and two individuals contributing multiple tumors. We cannot claim original Ng's count-based edgeR/TMM quasi-likelihood DE or proteome limma DE by analyzing `HCCfull.rds`; normal reference, preprocessing lineage and historical GRCh37 annotation must be audited separately.

## Interpretation gate

- **Workflow PASS:** technically reproducible descriptive correlation distribution. **Never:** biochemical mechanism, etiologic interaction, independent replication, causal treatment target, prospective clinical biomarker or original Ng 5,481-gene reanalysis.
- **V03 remains PARTIAL** until full batch/imputation/provenance review. **V07 remains NOT_STARTED** until real correlation artifact is verified, then at most **PARTIAL** (descriptive only).
- No per-feature inferential p-value, multiplicity correction, CI or clinical event analysis in Stage24.
- Next stage: inspect the successful Stage24 artifact, compare patient-mean and singleton-only distributions, verify Suter S5 Appendix and original Ng normal-reference transformations, then decide whether to perform FDR-controlled RNA/protein concordance with independent replication. No automatic manuscript biological PASS.

### Reproducible sources

- Ng full paper: https://www.nature.com/articles/s41467-022-29960-8
- Suter paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009767
- Source code PR: https://github.com/kimtk94/IS_Analysis_V3/pull/251
- Stage24 Actions: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173073247
