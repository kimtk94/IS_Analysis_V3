# HCC Multi-omics Feasibility — Evidence Gate (2026-10-30)

**Status: HOLD for definitive mechanism/recurrence study design; GO for a metadata-only pilot.** This is a research feasibility audit, not a claim that all omics data have been downloaded or joined.

## Existing related work
- Radiogenomics is **a separate project**: Google Drive `HCC_Radiogenomics_Public` and [Codex PR #228](https://github.com/kimtk94/Codex/pull/228) (still open at audit).
- This branch adds **RNA/protein + etiology + clinical + single-cell/spatial** dataset qualification to the research repo, avoiding changes to the existing IS/CKD pipelines.
- The 2026-10-22 HCC Research Proposal template reportedly prohibits AI-written submission text. This document is **internal technical evidence**, not a submission.

## Key findings and public-source cautions
1. **Ng et al 2022**: 114 patients / 122 tumors; RNA across 122 tumor samples, but proteome/phosphoproteome on **51 tumors** (paper-level intersection). Its RNA sequencing is **EGA controlled access** (EGAS00001005074); public PRIDE datasets PXD025705/PXD025836 are proteome/phosphoproteome, respectively. Supplementary/source data may enable a restricted pilot without raw RNA, but actual sample-level accessible matrices and mapping are **unverified**. https://www.nature.com/articles/s41467-022-29960-8
2. **Gao et al Cell 2019**: 159 HBV-HCC patients paired tumor/adjacent tissue in paper; public processed matrix access and RNA–protein linking not yet verified. https://pubmed.ncbi.nlm.nih.gov/31585088/
3. **TCGA-LIHC**: potential RNA/clinical, plus RPPA which is a *targeted antibody-based assay*, not whole MS proteome. Need to inspect etiology and clinical completeness. https://portal.gdc.cancer.gov/projects/TCGA-LIHC
4. **GSE14520**: 488 *microarray samples* across two platforms, not 488 unique patients; RFS/time/event columns must be verified from the actual supplement. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520
5. **GSE149614**: 10 patients / 21 GEO samples; processed scRNA matrix exists, raw data EGA controlled. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614
6. **GSE242889**: 5 patients (3 MVI+, 2 MVI−), 10 tumor/adjacent samples. MVI is a proxy biological phenotype, **not observed recurrence**. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242889
7. **GSE319709**: 18 scRNA sample records from PBMC/tumor/peritumor, HBV/MASH, CD45+ only. Cannot infer malignant hepatocyte origin. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319709
8. **GSE238264**: 7 processed spatial samples **after neoadjuvant cabozantinib/nivolumab**; not a treatment-naive external recurrence cohort. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238264
9. **PXD050969**: iProX proteome/lipidome/metabolome atlas; public raw matrix access and matching to RNA are unverified. https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD050969

## Criteria (pilot engineering thresholds, NOT universally accepted clinical criteria)
- At least one observed, documented **patient/specimen-level RNA↔protein join with N≥30** (paper abstract is insufficient).
- At least two verified etiology categories **each N≥10 within those joined data**. Need inspect overlap and confounders (stage, treatment, batch, platform).
- For optional recurrence prediction: one eligible event-defined dataset with ≥30 events, **independent** external cohort with ≥20 events. These are *minimum exploration triggers*, not statistical power guarantees.
- RFS, DFI, PFI, OS, MVI, and treated response are different endpoints. TCGA-CDR recommendations: https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/

No `verified_*` flag may be changed without a non-sensitive supporting audit artifact with counts and join rules. Never commit raw restricted files or identifiable patient keys. `HOLD` = evidence missing, not `NO GO`.

## CI usage
```bash
python -m unittest discover -s tests -p 'test_hcc_feasibility.py' -v
python scripts/hcc/feasibility_audit.py --out-dir results/hcc/feasibility
# Optional metadata URL reachability only (NOT a file download):
python scripts/hcc/feasibility_audit.py --probe-links --out-dir results/hcc/feasibility
```
GitHub Actions: `hcc-multiomics-feasibility.yml` on PR / branch main change / manual dispatch; writes `audit.json` and `audit.md` as artifacts. Success of CI means the guardrails and generator ran; **not** biological/clinical data availability.

## Next concrete tests (in order)
1. Review Ng 2022 paper Source Data and PRIDE file listing; establish whether 51 tumor RNA & protein values can be paired **without** EGA controlled data. If no, request authorized data and consider alternative public pairs.
2. Download **metadata only** from Gao/TCGA/GSE14520; calculate distinct patient counts, eligible HBV/HCV/MASH vs alcohol categories, missingness, outcome definitions.
3. Check outcome follow-up/event counts for eligible discovery/validation cohorts, separated by training/validation and primary therapy.
4. Only after stage 1–3, design low-cost RNA/protein pathway concordance pilot; scRNA/spatial are contextual validation, not artificially patient-matched modalities.

Repository registry: `research/hcc/feasibility/datasets.json`.
