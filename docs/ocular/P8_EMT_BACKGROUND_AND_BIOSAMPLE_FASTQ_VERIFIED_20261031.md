# Ocular P8 — EMT background specificity, BioSample XML, bounded shotgun FASTQ probes (2026-10-31)

## Immutable execution evidence

- **[Actions #41181407360](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41181407360): 3/3 matrix jobs SUCCESS**, with all three mode-specific artifacts actually downloaded and inspected, not merely assumed from green workflow status.
- P8a [EMT artifact #13035460980](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41181407360/artifacts/13035460980).
- P8b [BioSample XML artifact #13035176752](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41181407360/artifacts/13035176752).
- P8c [bounded FASTQ artifact #13034852571](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41181407360/artifacts/13034852571).
- Source code `scripts/ocular/ocular_p8_specificity_biosample_fastq.py`, reusable on-demand three-way workflow `.github/workflows/ocular-p8-specificity-biosample.yml`. Lightweight ocular PR CI continues automatic. Only public gene-level aggregate and **counts of metadata field types** are uploaded. No raw donor counts, patient identities, accession-to-patient mappings or FASTQ sequence snippets are stored.

## P8a — EMT direction specificity versus non-EMT background

Fixed source P7a [run #41179413138](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179413138), 12,250 unique matched Ensembl gene symbols from original 500nM 48h BeamA edgeR vs published 100nM 48h bulk TM.

| Fixed subset | Same sign | Denominator | Fraction | Wilson 95% interval *under independent-gene approximation* |
|---|---:|---:|---:|---|
| MSigDB Hallmark EMT | **141** | **194** | **72.68%** | 66.02–78.46% |
| All non-EMT background | **6,810** | **12,056** | **56.49%** | 55.60–57.37% |
| Difference | | | **+16.19 percentage points** | Not an independent-gene inferential CI |

An exploratory 2x2 gene-count Fisher exact test yields **odds ratio 2.0494, one-sided nominal p=2.75e-6**. This **must not be reported as a valid independent pathway-replication p-value or used as a causality claim**: genes within EMT and the background are co-expressed, and cross-study cohorts/dose/platforms are different. The Wilson CIs also adopt independent-gene assumptions and are diagnostic only.

**Both-cohort absolute log2FC thresholds**, applied to each gene in the original 194 EMT set (thresholds chosen *post hoc*):

| Minimum abs(log2FC) in both studies | Genes passing | Matching directions | Fraction | Both-study FDR<0.05 | Both-study significant & same sign |
|---:|---:|---:|---:|---:|---:|
| 0 | 194 | 141 | 72.7% | 49 | 46 |
| 0.25 | 75 | 68 | **90.7%** | 47 | 45 |
| 0.50 | 34 | 31 | **91.2%** | 24 | 23 |
| 0.75 | 23 | 21 | **91.3%** | 17 | 17 |
| 1.00 | 12 | 10 | 83.3% | 8 | 8 |

These post-hoc direction fractions are **conditional on effect-size thresholds in both cohorts** and do not create a held-out discovery. The 46 two-cohort FDR-positive signed-concordant genes in the output are exploratory molecular prioritization evidence; they **are not** 46 validated causally regulated proteins or 46 clinical biomarkers.

Examples of strong same-direction effects (500nM selected BeamA / 100nM bulk TM log2FC):
- NTM −2.518/−2.366; VCAM1 −1.690/−1.399; MXRA5 −1.375/−2.141.
- IGFBP3 −1.851/−1.277; TNFRSF11B −2.168/−1.269; CXCL12 −1.368/−1.228.
- WNT5A +0.996/+1.347; LUM −2.507/−0.990; **TGFBI −1.869/−0.843**.
- Divergent within the 49 both-study FDR-significant EMT genes: **SDC1, ITGA2, FAS**.

**GO:** direction-consistent ECM pathway and gene-level candidate prioritization. **HOLD:** a properly randomized 100/500nM donor-matched factorial study, functional ECM deposition/contractility/aqueous-outflow validation, exclusion of treatment-linked library preparation and marker protein confirmation. **STOP:** clinical steroid response prediction.

## P8b — live ENA BioSample XML audit [62/62 accessible; clinical group mapping HOLD]

Official ENA Browser API `/ena/browser/api/xml/{sample_accession}` returned **62/62** publicly available XML sample records for `PRJEB74248`, with no fetch errors. We read `SAMPLE_ATTRIBUTE/TAG/VALUE` internally and uploaded **counts only**, never raw attribute values or sample identifiers. Explicit named disease/diagnosis/treatment/health-state/control attribute screen:

- Explicitly assignable glaucoma cases: **0/62**.
- Explicitly assignable healthy controls: **0/62**.
- Ambiguous class: **0/62**.
- **UNKNOWN: 62/62**.
- Explicit negative-blank-attribute count: **0/62**.

This is stronger evidence of an *available-metadata labeling gap* than P7b title/description screening because the source sample attributes were checked directly. It **does not prove the authors lacked clinical records or blanks**, nor does it exclude unpublished sample group crosswalks or clinical data in other supplements. The paper describes 16 glaucoma and 16 healthy subjects, but its 62 runs should not be treated as 62 independent participants.

**HOLD** any differential taxa analysis, eye-drop/steroid stratification and microbiome→tear-protein linkage until case/control, donor-eye nesting, medication and low-biomass negative-control sample mappings are verified externally.

## P8c — bounded on-wire FASTQ accessibility [2/2 GZIP header PASS]

Representative `PRJEB74248` `submitted_ftp` files selected deterministically from both categories: a run with direct ENA `fastq_ftp` available, and a run with **submitted-only** FASTQ metadata. HTTPS requests used `Range: bytes=0-65535`, read **at most 65,536 bytes each** and closed connections; the original byte contents were never saved/uploaded.

- **2/2 HTTP 206 Partial Content** (server honored requested range).
- **2/2 GZIP magic bytes verified**, first 65,536 bytes accessible.
- **0** whole-file downloads; **0** full gzip CRC or MD5 checksums, **0** paired FASTQ parsing, **0** host read depletion.
- Source raw submission manifest still lists 62 runs, 155 FASTQ.gz files, and 155 MD5 records totaling approximately **528.9 GB**, so full-scale shotgun analysis would be a separate resource-intensive project.

**GO** network range readability/format magic for representative files, **HOLD** complete FASTQ integrity, R1/R2 pairs, host DNA contamination and group-aware differential microbiome analysis. No microbial role in steroid-induced glaucoma has been tested by these endpoints.

## Research conclusions and gates

P8 adds controlled **numerical** evidence that high-magnitude EMT genes often have consistent signed effects in two non-identical TM experiments, together with a concrete infrastructure audit showing that shotgun source reads can be accessed. It **does not bridge the cohorts**: the microbial study is ocular surface glaucoma patients, the DEX–TM study is cultured trabecular cells and has no patient steroid ΔIOP phenotype. Continue the drug-response mechanistic axis using orthogonal functional validation and keep the microbial cohort as a *separate, presently clinical-metadata-HOLD* project. Do not impute case/control labels by sample order, run parity, or paper figure numbering.
