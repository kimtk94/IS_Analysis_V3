# Ocular P6a/P6b/P6c — external factor gates, actual 3-way Actions results (2026-10-31)

## Execution and immutable evidence
- [**GitHub Actions #41173047662**](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173047662) **3/3 matrix jobs SUCCESS**.
- P6a [artifact #13031786504](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173047662/artifacts/13031786504): `P6a_9_library_pairs_public_qc.tsv`, `P6a_contamination_batch_gate.json`.
- P6b [artifact #13032041651](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173047662/artifacts/13032041651): `P6b_9_hallmark_ranked_gene_set_enrichment.tsv`, `P6b_12_candidate_hallmark_membership.tsv`, `P6b_ranked_hallmark_summary.json`.
- P6c [artifact #13032016694](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173047662/artifacts/13032016694): `P6c_3_study_ena_metadata_verification.json`.
- Reusable code: `scripts/ocular/ocular_p6_external_factors.py`, on-demand matrix `.github/workflows/ocular-p6-threeway-gate.yml` (one-time PR trigger disabled following successful launch).
- All three artifact ZIPs directly downloaded and contents reviewed. No private donor identifiers or individual read-level outputs uploaded.

## P6a — culture and technical contamination: input QC PASS; contamination/causality HOLD

Reconstructed, previously checksum-verified `source_qc.json` from P5q source GSE313736:
- **18/18** original public 10x libraries checked, exactly **9 pools A–I ×2 treatment libraries**.
- One and only one vehicle library and one DEX library *per pool*: treatment assigned to a distinct GEM 10x library. Pool/donor pairing does **not** eliminate preparation-batch confounding. Source counts and barcodes do not imply independent randomization.
- DEX/vehicle selected barcode count ratios across nine pools **0.67459–1.44201**; these differences describe composition/measurement, not causal DEX effect.
- No mycoplasma PCR status, bacterial culture sterility, LPS/endotoxin assay, or condition-balanced independent library-prep randomization were demonstrated by the **aggregate QC artifact**. Absence from this artifact **cannot prove** the original lab did not run such tests. Do not state contamination exists or was absent.
- Human-focused 10x 3' poly(A) scRNA-seq counts are **not a 16S assay**. Its negative bacterial signal is not a sterile culture certificate.

**Decision:** GO sequencing/count-library provenance; **HOLD** contamination-negative and DEX-vs-batch attribution. Next evidence: original lab SOP/QC records and experimental design, or independently designed matched randomization experiment.

## P6b — full-universe published Hallmark gene-set test: exploratory GO

P5b full edgeR aggregate **16,050 tested Ensembl gene rows**, duplicate gene symbol collapse yields **16,045 unique symbols** (5 collisions resolved deterministically). The official public MSigDB **2025.1.Hs Hallmark GMT** was downloaded at runtime and SHA256 checked: `f22066af72e215ccb7b89d88e492c07e1eef17534c2ca7b0f9902cfecbbdd8e9`.
- 9 *preselected for this stage* Hallmark families; score signed `-log10(PValue)` capped at 50, by primary DEX–Vehicle effect sign.
- Weighted running-sum enrichment score (ES) against 999 fixed-seed **random-gene-set-label assignments of the same size**. Two-sided empirical p minimum resolution **0.001**; Benjamini-Hochberg correction **across only these 9 examined sets**. This is **not a sample-label or donor-wise permutation test**, does not preserve gene–gene dependence and cannot establish randomized treatment causality.

| Published gene set | ES | Empirical p | BH-FDR across 9 |
|---|---:|---:|---:|
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | **-0.44713** | **0.001** | **0.009** |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | **-0.38868** | **0.002** | **0.009** |
| HALLMARK_HYPOXIA | -0.28091 | 0.049 | 0.147 |
| HALLMARK_TGF_BETA_SIGNALING | +0.31103 | 0.240 | 0.360 |
| HALLMARK_INFLAMMATORY_RESPONSE | -0.27406 | 0.163 | 0.3456 |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | -0.21448 | 0.579 | 0.579 |

Other sets APOPTOSIS, IL6/JAK/STAT3 and reactive oxygen species are also non-significant. **Do not equate this result with LPS, bacterial taxa, a verified DEX mechanism or direct protein-level ECM change.** In a rank test using gene-level statistics, the negative ES indicates rank clustering toward downregulated DEX–Vehicle side *for this particular dataset and selection*.

Exact overlap of 12 previously prioritized opposite-sign candidates with the selected 9 MSigDB sets (membership only, no protein mechanism):
- **CHAC1**: UNFOLDED_PROTEIN_RESPONSE.
- **EGR3**: TNFA_SIGNALING_VIA_NFKB and APOPTOSIS.
- **LPAR1**: INFLAMMATORY_RESPONSE.
- **SDC1**: EPITHELIAL_MESENCHYMAL_TRANSITION.
- Other 8 genes: not members of this *selected nine* Hallmark sets (not evidence they lack other pathway relationships).

**Decision:** GO numerical reproducibility and exploratory gene-set prioritization; **HOLD** causal interpretation, clinical relevance, 100 nM vs 500 nM attribution, and validation in a genuinely independent cohort.

## P6c — actual ENA raw-run metadata visibility: GO accession availability, HOLD clinical alignment

ENA official `read_run` report queried live for each exact accession; actual response content SHA256 saved. Data are **metadata only**, no FASTQ downloaded.

| Study BioProject | Platform per ENA | Runs / unique BioSamples | Runs with `fastq_ftp` links | Data suitability |
|---|---|---:|---:|---|
| [PRJNA1145861](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1145861) | AMPLICON (conjunctiva) | **139/139** | **139/139** | conditional GO for independent glaucoma/DED ocular microbiome reanalysis |
| [PRJEB74248](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB74248) | WGS (conjunctival shotgun) | **62/62** | **31/62** | P6e FOLLOW-UP: **62/62 have `submitted_ftp` paths**; binary download, file type and patient joins remain HOLD |
| [PRJNA1416399](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1416399) | AMPLICON (oral wash) | **129/129** | **129/129** | conditional GO for separate oral microbiome glaucoma project |

Corrections/controls:
- The previously mentioned 'PRJNA1180380' **is not the BioProject accession** for tear proteome/ocular shotgun paper; **PRJEB74248** is the valid accession, while `1180380` is NCBI's internal numeric BioProject ID.
- A lack of `fastq_ftp` is **not** proof the data are absent from SRA, ENA `submitted_ftp`, CRAM or controlled-access sources. A second accession-level audit is required.
- 62 ENA runs are **not** the number of independent glaucoma patients; no per-patient case/control or repeated-eye mapping has been validated.
- Shotgun data being public does **not** prove the paper's tear proteome is independently deposited or sample-paired at public access level.
- Zero verified matching primary GSE313736 DEX-treated TM donors, zero validated donor-level clinical steroid IOP response, no verified per-sample antibiotics/glaucoma eye-drop adjustment in these metadata results.
- For low-biomass ocular microbiome, experimental extraction and sampling negative blanks must be independently confirmed; the paper alone is insufficient for raw-data contamination correction.

## P6d/P6e verified follow-up (supersedes earlier metadata-only HOLD on missing `fastq_ftp`)
- [P6d cross-study pathway gate](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173246350): **12,251 mapped Ensembl IDs**, 9 common-background Hallmark tests, **5/9 signed pathway scores agree**; EMT signed **-/-** but UPR **-/+**. These are descriptive pathway effect rankings, not replicated pathway p-values.
- [P6e ENA alternative archive route gate](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173368834): The 31 runs without FASTQ links have `submitted_ftp`; **62/62 WGS runs have a submitted archive link**. This resolves metadata route discovery, **not** binary file access or clinical sample annotation.
- Full artifact-reviewed results in `P6DE_CROSS_STUDY_AND_SHOTGUN_ARCHIVE_VERIFIED_20261031.md`.

## Decision matrix
- **GO P6a source library QC and contamination-disclosure; GO P6b exploratory pathway ranking; GO P6c three metadata accessions.**
- **HOLD** sterile culture assertion, DEX mechanistic specificity, low-biomass blank subtraction, drug/age/disease covariates, matched tear proteomics, and published 100nM vs 500nM pathway replication.
- **STOP** microbial causal role in steroid-induced glaucoma and patient steroid response prediction from cross-sectional microbiomes and cultured TM without linked clinical outcomes.

## Next proposed computational gate
A **P6d independent-study pathway direction comparison** can test whether official 100nM TM 48h bulk gene-expression pathway signs corroborate the P5b 500nM BeamA Hallmark signal. It is only cross-cohort descriptive comparison, not matched dosing or a causal replication.
