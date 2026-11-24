# Mehrotra IOVS 2026 independent SRA — treatment timepoint and raw count readiness gate (2026-11-01)

## Frozen verified results

**Actual source audit:** [GitHub Actions #41217462028](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217462028) **SUCCESS**, [artifact #13049025050](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217462028/artifacts/13049025050). Artifact downloaded and ZIP CRC, exact underlying JSON and 23 distinct paired BioSample audit records inspected.

- **Original primary article:** [Mehrotra et al. *IOVS* 2026 67(2):13, PMC12882097](https://pmc.ncbi.nlm.nih.gov/articles/PMC12882097/), [official full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML), SHA256 `87a80692f387318037da1d6b7aa45b7a7e546e888abac92a37796e1b274f64c2`. The latest erratum [PMC13712379](https://pmc.ncbi.nlm.nih.gov/articles/PMC13712379/) adds only a Data Availability line with [SRA PRJNA1243143](https://www.ncbi.nlm.nih.gov/sra?term=PRJNA1243143); no numeric source changes stated.
- **Official ENA data access**: 46 unique RNA-seq paired-end RUNs, 92 gzipped FASTQ files, complete 92/92 MD5 metadata, 89,677,930,411 compressed bytes (~89.68 GB). The original archive has **23 nominal matched Vehicle/DEX sequencing-library pairs** (16 TM + 7 SCE), with **9 TM + 4 SCE distinct strain-label families**. These are *not* 23 independent donors.
- The original paper reports the following **experimental design sample counts**, NOT public-run-to-timepoint IDs: 1h = 3 TM +3 SCE strains, 6h = 4 TM +4 SCE strains, 2d/48h = 10 TM +5 SCE primary strains screened (9 TM +4 SCE pass source paper QC). Explicitly note the publication's screened vs QC-passed counts. These numbers cannot be summed and assumed to match the SRA sample totals because the source record has 23 candidate pairs, not 27 expected stage-specific pairs and repeated strains; original sample QC/source map may be incomplete.
- **Official XML BioSample examination**: 23/23 representative BioSample XML fetches PASS, **0/23** yielded a machine-identifiable explicit 1h/6h/48h treatment-time attribute. All 23 candidate arms had treatment confirmed before this audit (17 by exact aliases, 6 by 12 original BioSample fields). Thus DEX/Vehicle **GO**, RNA sample timepoint **HOLD**.
- **Official paper supplementary inventory**: the actual Europe PMC PMC12882097 XML declares **four** main supplements, `iovs-67-2-13_s001.docx`, `_s002.docx`, `_s003.xlsx` and `_s004.docx`. Paper XML does **not directly link to a source-level per-sample raw integer-count matrix**. It reports Broad RNA-Seq 2.0, alignment **STAR 2.4.1a / GRCh37 / GENCODE 19**, gene quantification **featureCounts** and a paired-donor DESeq2 analysis, but the direct paper gene-effects workbook is **not proven to be sample-level integer count data**.

## Supplementary file download retries — separate transport from absence

The first [supplement retrieval audit run #41217825435](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217825435), [artifact #13048932753](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217825435/artifacts/13048932753), inspected the two official PMC individual-file URL forms for all four files: **all four returned HTML/non-Office bytes or HTTP 404, so none was actually parsed**. It cannot establish lack of count tables or publication metadata. Source filenames exist in the primary XML. To avoid this false absence inference, the source verifier was corrected to first request the official Europe PMC `/PMC12882097/supplementaryFiles` ZIP endpoint, then fall back to individual routes. [Retry CI #41218079101](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218079101) **SUCCESS**; all four source files were actually retrieved and examined via official Europe PMC ZIP. The remaining HOLD is *source-specific run/timepoint evidence and sample-level raw counts*, not supplementary-file availability.

## Subsequent official supplement access and candidate pairing — ACTUAL VERIFIED UPDATE

The first PMC individual-file URL errors were **transport/access issues**, not missing supplements. The alternate official Europe PMC ZIP [run #41218079101](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218079101) **SUCCESS**; [artifact #13049046005](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218079101/artifacts/13049046005) downloaded and independently inspected. The ZIP **23,188,342 bytes**, SHA256 `47adf77cbb2472003f222c04e4a626b742a434dd3aa453d009ba3dc8a1966fd2`, contains and passes CRC for all four official source files (3 DOCX, one XLSX). The official `s003.xlsx` is **16,686,303 bytes**, SHA256 `3c30e93461812d05f42a1978b8e446235709f2283e11f620c212beb26478698d`.

The 63-sheet source workbook contains author-published differential-expression and pathway tables for each timepoint; `SuppTable82_TM_MarkerExpression` and `SuppTable83_SC_MarkerExpression` additionally enumerate the **individual strain/timepoint-specific CONTROL TPM columns**. No clearly labeled complete raw integer gene-by-46-sample count matrix was detected among those named sheets. This is a narrow inventory result, **not proof that source featureCounts files exist nowhere else**.

The per-strain source marker columns give:
- **TM 48h/2Days (9):** TM123, TM134, TM135, TM137, TM140, TM151, TM155, TM92, TM96.
- **TM 6h (4):** TM135, TM151, TM92, TM96.
- **TM 1h (3):** TM123, TM151, TM96.
- **SCE 48h/2Days (4):** SC67, SC86, SC89, SC91.
- **SCE 6h (4):** SC67, SC77, SC89, SC91.
- **SCE 1h (3):** SC86, SC89, SC91.

Using only the frozen publisher-derived marker-list headers and exact ENA run aliases, the inexpensive [temporal candidate crosswalk Actions #41218779997](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218779997) **SUCCESS**, [artifact #13049575131](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218779997/artifacts/13049575131) downloaded, three output files and ZIP CRC verified. The **23 candidate paired RNA libraries** distribute in a pattern fully consistent with the supplement's strain-by-timepoint combinations: **TM 48h 9 pairs, TM 6h 4, TM 1h 3, SCE 48h 4, SCE 1h 3**. No SCE 6h pairs are identifiable in the **current 46-RUN accession**, despite four SCE 6h strains in the published expression supplement. This apparent archive-coverage difference should be documented; do not assume those four are unavailable in every repository, or that the archive is complete.

**Critical scientific distinction:** the SRA BioSample entries themselves contain **0/23 directly verified treatment-time attributes**. Consequently, these 23 assignments are **source-consistent *candidate* mappings using naming conventions plus strain/timepoint set intersection**, not an author-confirmed, sample-specific run/time map. The 48h TM candidate cohort comprises **9 paired library sets = 18 RUNs = 36 FASTQ files**, totaling **32,707,309,685 compressed bytes (~32.71GB)**, as opposed to 46 RUNs / 89.68GB in the broader source project. See artifact `MEHROTRA_9_TM_48H_CANDIDATE_RUN_PAIRS_HOLD.tsv` for exact accession pairs; do not claim independent raw-count replication until timepoint and QC sample assignment are positively corroborated.

A first lightweight candidate crosswalk [run #41218660217](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218660217) failed an internal zero-Count comparison (`SC_6Hrs` 0 is omitted from a sparse Counter), not source-metadata inconsistency. It was corrected; the *successfully validated* second run #41218779997 is the authoritative artifact.

**Updated scope recommendation:** author-published 48h/6h/1h DEG tables are accessible, so they can be compared without re-aligning FASTQs. To promote a 48h TM *raw-count* validation, obtain authors' sample sheet or an explicit run ID ↔ timepoint mapping; then check source-level featureCounts integer tables before deciding on any 32.7GB pilot. A large-scale 89.7GB download remains STOP without that verification.

## Scientific GO / HOLD / STOP

| Checkpoint | Status | Meaning |
|---|---|---|
| Primary 2026 Erratum content | **GO** | Adds source SRA project only, no stated DEG or source-table correction |
| 46 SRA RNA-seq RUNs and 92 FASTQ locations/MD5 | **GO metadata** | Raw study exists; no FASTQ downloaded |
| 23 nominal treatment-control library pairs | **GO** | 17 exact alias pairs plus six pairs documented by official BioSamples |
| Assign each run to 1h/6h/48h and TM/SCE QC-passed cohorts | **REVIEW candidate / HOLD direct provenance** | Supplement-supported strain/timepoint and alias-style match 23/23; no per-RUN timepoint attribute in ENA |
| Determine whether already deposited *sample-level integer counts* exist | **HOLD** | All four supplemental files parsed; 63 author gene-DE/pathway sheets and control TPM examples, no full per-sample integer count sheet identified |
| New independent 100nM 48h TM edgeR/DESeq2 | **HOLD** | Need exact timepoint, strain identity, inclusion/QC and counts before treating any file as an independent human donor |
| Download and reprocess all 89.68GB on spec | **STOP for now** | Unjustified expense and high risk of incorrect sample inclusion/time contrast |
| Use frozen Mehrotra publication values in descriptive three-study Figure | **GO** | Corrected erratum did not claim a numerical update; keep source-file provenance |

## Next narrowly scoped work

First positively identify an original *sample name → timepoint (1/6/48h) → TM/SCE → strain → treatment → QC pass* crosswalk in supplementary files, deposited source metadata or author-supplied sample sheet. If a direct 48h TM *raw integer-count* matrix exists, use it with paired strain DESeq2 instead of re-aligning 90GB FASTQ. If only FASTQ is available, audit one matched 48h TM pair and pilot QC before scaling a reproducible GitHub-hosted processing plan. If a source timepoint map cannot be obtained, retain the published table signed comparison and do not claim an independent raw-data replication.

**Study contribution caution:** even a new 48h TM re-fit cannot break the original GSE313736 treatment-specific library confounding, and the original source's 63↑/17↓ named gene identity list remains independently unverified.
