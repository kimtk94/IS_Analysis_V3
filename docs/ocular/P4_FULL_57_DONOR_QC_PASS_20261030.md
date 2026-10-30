# Ocular P4 full raw-count pseudobulk verification — 2026-10-30

## Actual execution (not just a planned pipeline)
- [GitHub Actions successful full-cohort run #41069563630](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41069563630)
- [Aggregate-only evidence artifact ID 12984842522](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41069563630/artifacts/12984842522)
- Status **GO_57DONOR_FULL_RAW_PSEUDOBULK_QC**: all 18 original GSE313736 10x matrices physically downloaded and completely streamed, full gzip CRC / NNZ and feature order verified; raw UMI counts summed by donor and arm in runner memory only, not committed.
- Primary eligibility pre-registered: `SNG`, Chen2025 fine `BeamA`, score >=0.75, >=30 cells per condition/donor. **57 paired biological donors**, with 57 preserved after all actual barcode joins.
- Metadata eligible cell rows: control 40,115 / DEX 46,454. Actual source-matrix-matched cells: **control 40,066 / DEX 46,011**, 86,077 total.
- Pseudobulk UMI sum: control **1,018,027,543**, DEX **909,649,764**; combined 1,927,677,307 integer UMIs.
- Median UMI per donor: control 17,667,495 / DEX 15,397,698. Median detected genes per donor: control 22,999 / DEX 23,767.
- 36,601 aligned feature IDs; 33,222 genes detected in at least one condition; 26,220 genes expressed in at least 10 donors in **each** arm (descriptive coverage only, not the formal edgeR `filterByExpr` set).
- All sparse entries from the 18 matrices verified: **1,751,655,185** total source nonzeros. No individual donor-level expression vectors or donor identifiers uploaded.
- **Only QC exception:** GSM9373291, Pool G DEX, selected BeamA barcode recall **92.743%**; all other 17 library recalls >=95%. Full P4 initial attempt used an over-strict >=95% gate and stopped at Pool G; independent preflight [#41069377720](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41069377720) showed 57/57 paired donors preserved, all 18 recall >=90%, and prompted a documented pre-analysis threshold update. Corrected whole-cohort run passed; keep Pool G leave-out sensitivity mandatory.

## Explicit research interpretation and limits
- **GO** primary donor-blocked BeamA raw-count pseudobulk differential-expression *testing*. Model to follow `docs/ocular/P5_PAIRED_BEAMA_STATISTICAL_PLAN_20261030.md`.
- **NOT DONE:** edgeR / DESeq2 count-based negative binomial modeling, FDR, gene hits, effect heterogeneity estimates, 9-pool sensitivity, independent treated/control replication, clinical genotype-based steroid IOP association.
- **HOLD** declaring novel pathogenic mechanism, treatment-response predictor, or patient adverse effect until above inferential and external validation gates pass.
- Treatment and 10x sequencing libraries are different for each arm; donor blocking alone cannot rule out technical batch confounding. Assess replication across pools and independent 100 nM/48h bulk study (GSE313736 used nominal 500 nM/48h).
- PRJNA1243143 external metadata check [#41068968781](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41068968781) found 46 paired-end ENA runs, but sample-level treatment/time labels unavailable. Do **not** interpret the generic DEX word in experiment titles as proof all 46 are treated. Independent validation HOLD.
- The original GEO description reports 74 donors / 316,799 singlets vs final public 59 / 169,398 SNG metadata. This study attrition remains an audit item; 57 are **selected BeamA donor pairs**, not the full original cohort.
- Clinical glucocorticoid-induced ΔIOP GWAS from Lama et al. 2026 (PMID 42388255) already implicates UBE2E3/SSFA2; the clinical phenotype is not present in these cultured cells and controlled-access dbGaP is separate.
- Existing Mehrotra et al. 2026 100nM DEX RNAseq (PMID 41636427) already includes paired DESeq2 and ECM/adhesion response; simply rediscovering FKBP5 or FAM107A is **not novelty**.

## Next explicit acceptance gate
1. Process all 57 donor-paired raw integer pseudobulks in **one trusted environment** with edgeR TMM normalization, `filterByExpr`, paired donor design `~ donor + treatment`, robust QL dispersion, BH-FDR 0.05, QC plots and transcript provenance.
2. Verify rank/confounding of donor plus pool and 10x library effects; repeat leave-one-pool-out and donor leave-one-out on prioritized genes.
3. Replicate effect direction in PRJNA1243143 **only after** sample↔condition/time↔biological strain metadata are verified.
4. Verify novel findings beyond GSE313736's original DEX/sc-eQTL analysis and 2026 TM time-course papers before advancing to MR/coloc/clinical claims.
5. Keep raw individual-level UMI matrices, donor identities and patient data off public Git and public Actions artifacts. Published gene-level aggregate statistics require a separate data-sharing assessment.

## CI resource governance
The expensive full-cohort stream workflow was changed to manual `workflow_dispatch` after the successful trial; normal PR commits should run only lightweight checks. Source is `scripts/ocular/ocular_p4_fullcohort.py`, pilot `ocular_p4_pseudobulk_pilot.py`, preflight `ocular_p4_prescan.py`.
