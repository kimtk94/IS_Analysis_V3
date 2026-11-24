# IOVS 2026 source-specific contradiction and SRA 23-library triple-source audit

**Verified 2026-11-01**, [GitHub Actions run #41219957609](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41219957609) **SUCCESS**, [artifact #13049965603](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41219957609/artifacts/13049965603) actual ZIP CRC PASS and three summary/TSV tables directly checked.

## Exact immutable source records

1. Mehrotra et al., IOVS 2026 original [PMC12882097 official fullTextXML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML), SHA256 `87a80692f387318037da1d6b7aa45b7a7e546e888abac92a37796e1b274f64c2`, **Table 1: Demographic Data of Donors for Eyes From Which TM and SCE Cells Were Isolated**. Exact cells were parsed and cross-checked programmatically, not hallucinated.
2. Original [PMC12882097 official Europe PMC supplementaryFiles ZIP](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/supplementaryFiles), XLSX `iovs-67-2-13_s003.xlsx` SHA256 `3c30e93461812d05f42a1978b8e446235709f2283e11f620c212beb26478698d`, source **Supplementary Table 82 TM marker expression** and **Table 83 SCE marker expression**; per-strain CONTROL TPM column headers include 1h, 6h and 2Days.
3. SRA original source [PRJNA1243143](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1243143), 46 ENA RUNs, 23 source alias/control DEX candidate library pairs. Frozen previous [run #41218779997](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218779997), [candidate artifact #13049575131](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41218779997/artifacts/13049575131). **0/23** public BioSample representative source records explicitly identifies actual exposure hour, so the alias-based timepoint is **not a sample-specific experiment log**.

## Triple-source findings

| Cell lineage | Time | Original paper Table 1 strain count | Supplement S003 named strains | Source agreement |
|---|---|---:|---|---|
| TM | 48h/2Days | 9 | 9 | **PASS** |
| TM | 6h | 4 | 4 | **MISMATCH** |
| TM | 1h | 3 | 3 | **PASS** |
| SCE | 48h/2Days | 4 | 4 | **PASS** |
| SCE | 6h | 4 | 4 | **PASS**; none of these four has a matched 6h library pair in the 46-run SRA archive |
| SCE | 1h | 3 | 3 | **PASS** |

**Exact conflicting names in TM 6h:**
- Table 1 **only**: `TM155` (the published demographics row lists `TM155 | 2 d, 6 h`).
- Supplementary Table 82 **only**: `TM151` (the published gene expression marker-column list includes `TM151_CTRL_6Hrs`, whereas Table 1 lists `TM151 | 2 d, 1 h`).
- ENA alias-based 6h candidate `TM151.4` contains `SRR32897704;SRR32897705`, but neither RUN provides an independently verified timepoint. This particular candidate is **HOLD_PAPER_TABLE1_VS_SUPPLEMENT_DISCORDANT**, not author-confirmed timepoint/strain.
- There are **22/23 source-consistent candidate SRA library pairs** and **1/23 inter-document discordant candidate**. This is not a 22/23 direct time-attribute validation: **0/23 direct per-RUN duration fields**.

**Critical positive finding for thesis**: all nine **2-day TM** strains are unchanged across paper Table 1 and supplement S003 (TM123, TM134, TM135, TM137, TM140, TM151, TM155, TM92, TM96). The corresponding SRA candidate subset remains **9 paired library sets / 18 RUNs / 36 paired-end FASTQ files / 32,707,309,685 compressed bytes**. This is stronger source *consistency*, but **not** author-confirmed run-specific treatment duration or QC acceptance. These nine strains occur as source names in the article, and the 9+9 SRA libraries show source design consistency. Use no stronger label than **REVIEW_SAME_STRAIN_THREE_SOURCE_CANDIDATE_NOT_RUN_TIME_PROVEN**.

## Interpretation and paper-facing caveats

- The manuscript may state that the *published timepoint-specific strain sets are almost consistent*, with one demonstrable 6h inconsistency. **Do not silently 'correct' TM155 to TM151**; the correct row or sample metadata must be obtained from authors. The October 2026 erratum only added Data Availability; no correction of Table 1/Supplementary Table 82 is reported.
- The nominal 2-day 9-strain subset is supported by two parts of **the same research article**, not two independent sample-cohort studies. It therefore remains candidate run/time mapping.
- SCE6h strata and available SRA accession coverage differ; do not assume the four published 6h SCE samples must be in this archive.
- The complete source raw integer featureCounts matrix is not identified by the existing four official supplements. No FASTQ or donor-level counts were downloaded, and no independent raw replication results were generated in this gate.
- Original GSE313736 author 63↑/17↓ named gene list remains a **separate HOLD**; do not confuse the two source papers.

## Next precise gate

Request or find one **authoritative original run-to-time sample sheet**, with fields `RUN, BioSample, TM/SCE, strain, 1h/6h/48h, DEX/vehicle, QC inclusion, sample replicate/batch`. The single TM151-vs-TM155 6h discrepancy makes this particularly important. If genuine complete integer gene-counts are available, prioritize reproducing the nine **48h TM** pairs without aligning ~33 GB FASTQ; only afterwards consider a one-pair pilot. **Do not infer raw-read independent replication from source-consistent accession labels alone**.
