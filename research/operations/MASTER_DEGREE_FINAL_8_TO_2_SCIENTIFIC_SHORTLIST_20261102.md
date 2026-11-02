# MASTER DEGREE — 8 research projects → one primary + one conditional second (2026-11-02)

**Scope:** 2026-11-02 real original-study/output and recent public literature evidence review. This is NOT an arbitrary percent of scientific completion. Codes/metadata CI PASS do not certify new scientific discovery or graduate faculty approval. The recommendation below does not itself alter automated daily tasks; no project deleted, no PR merged.

## Recommended resource allocation

- **ACTIVE PRIMARY = GSE313736 DEX/TM glaucoma (computational analysis reliability/identifiability thesis).** Already completed source-count donor-paired edgeR, alternative 55-donor annotations, 3-study direction evaluation, 7d 3-strain leave-one-out, rank proof and actual figure QA. Remaining academic work is data-unit diagrams, selected candidate and FDR threshold reporting, manuscript package, advisor approval. **Not DEX molecular target discovery or clinical prediction.** Main manuscript at PR #253.
- **ACTIVE CONDITIONAL SECOND = Ischemic stroke QTL/GWAS, only if research question can be materially distinguished from prior proteome-wide MR/colocalization publications.** Verified exact official EAS AIS original GWAS 6,725,216 rows (source SHA/MD5, every row basic numeric QC), allele direction and p concordance 32/32 pre-selected non-palindromic SNP anchors; actual known QTL-GWAS matching includes FGF5 and CALHM2 candidates. The official five EAS stroke subtype study/ancestry accession identities all passed bounded CI. **No completed modern independent replication, new full data-approved MR or independent colocalization is inferred from the provenance audits.**
- **RESERVE ONLY, NOT ACTIVE THIRD = combined HCC tumor proteotranscriptomics + MRI radiomics.** Some real signals are genuine: Ng/Suter 48 unique patients, 5,343 exact one-to-one matched RNA/protein feature pairs and median per-feature Spearman 0.284086; MRI 14 registered PV patients, 168 source-reconstructed 2-reader×2-grid×3-bin extractions from official masks. But independent Gao PDC 330 quantitative TMT channels' original aliquot mapping and OEP000321 RNA-to-patient join remain unproven, and HCCDB stage58 fetched JS route chunks rather than RNA matrix (rna_matrix_downloaded=false). HCC-Radiomics MRI–RNA clinical matched 14 patient endpoint testing remains **NOT EXECUTED**. Without those joins, no three-modal HCC paper claim; avoid further expensive source-web scraping as if that were RNA data.
- **PAUSE / ARCHIVE:** stand-alone ophthalmology P0 (merge useful Mehrotra source metadata into glaucoma), muscle strength GWAS, OSCC/oral cancer, skin/cosmetic MR. Archive source links/code so these can be resumed, but do not represent source-only checks as biological discovery. No need to delete their repositories.

## Evidence comparison — 8 tracks

| Track | Verified latest real result | Scientific missing gate | Scientific maturity and decision |
|---|---|---|---|
| 1. Glaucoma TM/DEX | GSE313736 57 donor paired edgeR 16,050 genes/13,717 FDR hits; 55 donor BeamA/scPred, 78.258% sign agreement in small effects (n=5,349); MYOC↑, TGFBI↓, MMP2↓ in 3 source directions; strict EMT 90/194; 18 original arm-specific 10x library rank alias; main figures plus corrected S1 visually inspected | Causal drug vs 10x prep NOT identifiable, no clinical IOP; thesis program supervisor acceptance and composed reference proof | **HIGH completion as bounded METHODS/reproducibility master's; PRIMARY** |
| 2. IS plasma pQTL/eQTL MR+coloc | Original official GIGASTROKE EAS AIS 260,689,208-byte gzip 6,725,216 variants fully downloaded/checksummed and numerically screened in GitHub CI; signed original-to-archive crosswalk 32/32; 5 official EAS study accessions/traits/ancestries; archived 4-gene eQTL–GWAS crosswalk and earlier FGF5 evidence | Confirm original pQTL/QTL upstream models, effective N/ancestry-specific LD, eliminate BBJ-discovery cohort overlap, rerun fully validated MR/coloc, replicate in independent cohort; **existing 2025–2026 proteome-wide MR prior art overlaps** | **High future value, moderate current inferential readiness; CONDITIONAL SECOND** |
| 3. HCC matched RNA–protein | Suter original 48 patient means with 5,343 protein/gene pairs, Spearman median ~0.284; PDC 171 cases/331 aliquots and 330 nonreference TMT channels; independent feature-parity audit done | Same-person independent Gao RNA-to-PDC aliquot/source matrix not verified; original HCC paper covers many conclusions; no final gene-wise external biological reproduction | **Technical effort high, scientific independent join HOLD; RESERVE** |
| 4. HCC radiomics + RNA | Official registered portal venous MRI 14 patients with 2 expert ROIs, 2 voxel grids×3 bins×2 readers=168 extractions; GLCM contrast parameter changes large (e.g. 3mm binWidth25 median abs 86.59%) | Image-to-tissue timepoint and patient clinical eligibility + MRI–RNA correlation not executed; independent radiomics ROI validation missing | **GO imaging methods only, HOLD multiomics; merge with HCC, not standalone SECOND** |
| 5. Ophthalmology multiomics P0 | Official Mehrotra PRJNA1243143 46 ENA runs, 23 library pairs / 92 FASTQs, source metadata/erratum | Run-timepoint/biosample independence not resolved, raw counts external re-fit not performed | **SOURCE ONLY; integrate into glaucoma, pause as separate project** |
| 6. Muscle GWAS | GTEx v10 Muscle_Skeletal rsID→variantID defensive code, expression-QTL route knowledge | Missing exact 20 initial loci and actual GTEx variant source with ancestry-matched LD; no reproducible completed outcome effect | **PAUSE** |
| 7. Oral cancer OSCC | GSE103322 publicly documented 18 patients/5,902 single cells, viable reference | Donor-level expression/clinical metadata and independent patient-matched proteomics unverified, no tested claim | **PAUSE** |
| 8. Skin beauty MR | Bristol perceived-age GWAS public 423,992-person study with readme allele/effect schema and ~1GiB GWAS source | Full GWAS download/allele QC, applicable phenotype, instrument strength, independent MR/coloc unperformed | **PAUSE** |

### The critical IS novelty collision (updated 2026-11-02)

The 2025 Pang et al. Stroke report already uses EAS and EUR plasma proteomic MR with shared-signal colocalization for 2,923 proteins, naming **FGF5, F11, MMP12, FURIN, TMPRSS5, PROCR** among others (https://pubmed.ncbi.nlm.nih.gov/40304040/).

A 2025 BMC Neurology report already uses UKB-PPP plasma pQTL and deCODE validation, stroke subtypes and colocalization; FURIN/F11/DDHD2/VSIR were highlighted (https://pubmed.ncbi.nlm.nih.gov/40596960/).

**Particularly relevant newly published 2026-10-22 Song et al., Translational Stroke Research** performs cis-pQTL MR for AIS/LAA/CES/SVD, deCODE replication, Bayesian colocalization and a 65,479-cell human brain vascular single-nucleus atlas; includes MMP12, SWAP70, SCARA5 and SVD SPATA20/KDELC2 results: https://pubmed.ncbi.nlm.nih.gov/42911303/.

The 2026 Communications Medicine report uses >3,500 cis-pQTL protein profiles across platforms, MR and colocalization for ischemic stroke and subtypes: https://www.nature.com/articles/s43856-026-01734-z.

**Conclusion:** IS thesis MUST NOT be described as novel simply because it is cross-ancestry, stratified, uses deCODE, colocalization or brain scRNA. FGF5 etc overlap prior literature. The next IS scientific gate is a **claim-by-claim prior-art comparison** and a specific *new*, empirically falsifiable niche. Candidate niche is a rigorously revalidated East Asian GWAS/QTL source-lineage, allele/LD and cohort-overlap non-independence audit—NOT another simple druggable-target candidate list. Assess whether that is itself sufficiently novel for advisor's digital health genetics thesis; if not, **stop IS as an active second**, rather than proceed merely because of sunk effort.

## Conditional retention stop/go criteria

**Glaucoma PRIMARY:** advisor accepts a transcriptomic methodological robustness thesis with explicit treatment-library identifiability limitation. The final scientific claim cannot exceed signed transcript association, and the original author DEG80 ID list is non-blocking without investigator outreach. Deliverable is completed thesis draft in university approved format.

**IS SECOND:** retain only when all four conditions can be met from original/verified sources:
1. **Novelty** — a single primary result clearly not already demonstrated in Pang 2025, Zhao 2025, Song Oct 2026 and 2026 Communications Medicine; no generic "pQTL MR and coloc" novelty claim.
2. **Source** — original EAS file GWAS source and upstream cis-pQTL/eQTL are allele/assembly/LD/N/cohort-overlap validated. Whole-file QC and 32 SNP matches are excellent progress but insufficient alone.
3. **Inference** — independent validated MR/coloc or an explicitly declared *methods/provenance* study with falsifiable positive/negative controls and reproducible real-data results.
4. **University fit** — supervisor confirms this is a distinct, independently valuable second line of research, or decide to prioritize only one thesis.

If these gates fail, move IS to archive. **Do not automatically replace it with HCC** unless official independent same-patient RNA/proteomics or valid MRI–RNA matched clinical cohort is actually proved, rather than more HCCDB JS/API reconnaissance.

## Working model and automation

- Limit ACTIVE research work to **1 full thesis + at most 1 conditional high-risk project**, preserving all frozen results and code for the other 6.
- Stop treating repository count, commit count, CI green statuses, Stage numbers, pseudo-scientific percentage completion or web source discovery as biological completion measures.
- There is an existing daily eight-track task; **do not silently change/disable it without the user authorizing the revised coverage**. Once user confirms a shortlist, change the scheduled task from 8-track one-job-each to 2-track execution plus one compact weekly source gate for paused projects if needed.
- GitHub Actions is preferred for new real-data evaluation; new output must include source IDs, run IDs, falsifiable gate and independent-science status, not simply metadata PASS.

## Verified evidence anchors

- [Eight-track actual bounded source CI #41363141116](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41363141116): nine jobs SUCCESS, eight explicitly separate **scientific_gate=HOLD** from source probe PASS.
- [Eight-track prior code/metadata readiness #41215363085](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41215363085), all scientific gates HOLD at that earlier date; may be superseded by actual later GSE313736/IS original-data runs.
- [Glaucoma final claim/source/figure CI #41432752451](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41432752451) and [manuscript source check #41436138351](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41436138351), with archived corrected Figure.
- [IS original EAS source full-file QC #41436048672](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41436048672) and [32 non-palindromic signed GWAS anchors #41437991270](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41437991270), official 5 East Asian stroke subtype metadata [#41439023415](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41439023415). **No new MR/coloc from source audit.**
- [HCC Stage58 actual chunk CI #41439146697](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41439146697): 3/3 SUCCESS source JS hash, explicitly **rna_matrix_downloaded=false**, official patient RNA–PDC linkage=false.
- [MRI registered PV 14-patient radiomics factorial #41208617950](https://github.com/kimtk94/Codex/actions/runs/41208617950): 9/9 SUCCESS actual image parameter extraction; MRI–RNA association remains unexecuted.

**Decision owner:** user and academic advisor. This document recommends limiting work, not changing any scheduled automation or merging a research PR without permission.
