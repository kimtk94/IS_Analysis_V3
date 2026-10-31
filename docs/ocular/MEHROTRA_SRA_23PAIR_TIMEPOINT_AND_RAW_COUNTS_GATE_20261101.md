# Mehrotra IOVS 2026 independent SRA — treatment timepoint and raw count readiness gate (2026-11-01)

## Frozen verified results

**Actual source audit:** [GitHub Actions #41217462028](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217462028) **SUCCESS**, [artifact #13049025050](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217462028/artifacts/13049025050). Artifact downloaded and ZIP CRC, exact underlying JSON and 23 distinct paired BioSample audit records inspected.

- **Original primary article:** [Mehrotra et al. *IOVS* 2026 67(2):13, PMC12882097](https://pmc.ncbi.nlm.nih.gov/articles/PMC12882097/), [official full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML), SHA256 `87a80692f387318037da1d6b7aa45b7a7e546e888abac92a37796e1b274f64c2`. The latest erratum [PMC13712379](https://pmc.ncbi.nlm.nih.gov/articles/PMC13712379/) adds only a Data Availability line with [SRA PRJNA1243143](https://www.ncbi.nlm.nih.gov/sra?term=PRJNA1243143); no numeric source changes stated.
- **Official ENA data access**: 46 unique RNA-seq paired-end RUNs, 92 gzipped FASTQ files, complete 92/92 MD5 metadata, 89,677,930,411 compressed bytes (~89.68 GB). The original archive has **23 nominal matched Vehicle/DEX sequencing-library pairs** (16 TM + 7 SCE), with **9 TM + 4 SCE distinct strain-label families**. These are *not* 23 independent donors.
- The original paper reports the following **experimental design sample counts**, NOT public-run-to-timepoint IDs: 1h = 3 TM +3 SCE strains, 6h = 4 TM +4 SCE strains, 2d/48h = 10 TM +5 SCE primary strains screened (9 TM +4 SCE pass source paper QC). Explicitly note the publication's screened vs QC-passed counts. These numbers cannot be summed and assumed to match the SRA sample totals because the source record has 23 candidate pairs, not 27 expected stage-specific pairs and repeated strains; original sample QC/source map may be incomplete.
- **Official XML BioSample examination**: 23/23 representative BioSample XML fetches PASS, **0/23** yielded a machine-identifiable explicit 1h/6h/48h treatment-time attribute. All 23 candidate arms had treatment confirmed before this audit (17 by exact aliases, 6 by 12 original BioSample fields). Thus DEX/Vehicle **GO**, RNA sample timepoint **HOLD**.
- **Official paper supplementary inventory**: the actual Europe PMC PMC12882097 XML declares **four** main supplements, `iovs-67-2-13_s001.docx`, `_s002.docx`, `_s003.xlsx` and `_s004.docx`. Paper XML does **not directly link to a source-level per-sample raw integer-count matrix**. It reports Broad RNA-Seq 2.0, alignment **STAR 2.4.1a / GRCh37 / GENCODE 19**, gene quantification **featureCounts** and a paired-donor DESeq2 analysis, but the direct paper gene-effects workbook is **not proven to be sample-level integer count data**.

## Supplementary file download retries — separate transport from absence

The first [supplement retrieval audit run #41217825435](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217825435), [artifact #13048932753](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41217825435/artifacts/13048932753), inspected the two official PMC individual-file URL forms for all four files: **all four returned HTML/non-Office bytes or HTTP 404, so none was actually parsed**. It cannot establish lack of count tables or publication metadata. Source filenames exist in the primary XML. To avoid this false absence inference, the source verifier was corrected to first request the official Europe PMC `/PMC12882097/supplementaryFiles` ZIP endpoint, then fall back to individual routes. [Retry CI #41218079101](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218079101) was launched; scientific gate remains HOLD until its actual artifact confirms source-level file metadata or document content.

## Scientific GO / HOLD / STOP

| Checkpoint | Status | Meaning |
|---|---|---|
| Primary 2026 Erratum content | **GO** | Adds source SRA project only, no stated DEG or source-table correction |
| 46 SRA RNA-seq RUNs and 92 FASTQ locations/MD5 | **GO metadata** | Raw study exists; no FASTQ downloaded |
| 23 nominal treatment-control library pairs | **GO** | 17 exact alias pairs plus six pairs documented by official BioSamples |
| Assign each run to 1h/6h/48h and TM/SCE QC-passed cohorts | **HOLD** | No explicit timepoint in any 23 representative BioSamples |
| Determine whether already deposited *sample-level integer counts* exist | **HOLD** | Publication contains DE results; no raw-count matrix directly linked in paper XML |
| New independent 100nM 48h TM edgeR/DESeq2 | **HOLD** | Need exact timepoint, strain identity, inclusion/QC and counts before treating any file as an independent human donor |
| Download and reprocess all 89.68GB on spec | **STOP for now** | Unjustified expense and high risk of incorrect sample inclusion/time contrast |
| Use frozen Mehrotra publication values in descriptive three-study Figure | **GO** | Corrected erratum did not claim a numerical update; keep source-file provenance |

## Next narrowly scoped work

First positively identify an original *sample name → timepoint (1/6/48h) → TM/SCE → strain → treatment → QC pass* crosswalk in supplementary files, deposited source metadata or author-supplied sample sheet. If a direct 48h TM *raw integer-count* matrix exists, use it with paired strain DESeq2 instead of re-aligning 90GB FASTQ. If only FASTQ is available, audit one matched 48h TM pair and pilot QC before scaling a reproducible GitHub-hosted processing plan. If a source timepoint map cannot be obtained, retain the published table signed comparison and do not claim an independent raw-data replication.

**Study contribution caution:** even a new 48h TM re-fit cannot break the original GSE313736 treatment-specific library confounding, and the original source's 63↑/17↓ named gene identity list remains independently unverified.
