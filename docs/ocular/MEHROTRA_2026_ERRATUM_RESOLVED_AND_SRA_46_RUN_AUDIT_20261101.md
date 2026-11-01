# Mehrotra 2026 Erratum Verified and Newly Disclosed Original RNA-seq Evidence

**Audit date:** 2026-11-01 (KST). **Final official metadata audit:** [GitHub Actions #41214240628](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41214240628), **SUCCESS**, [artifact #13047557508](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41214240628/artifacts/13047557508). **Latest 12-BioSample mapping:** [run #41214740781](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41214740781) **SUCCESS**, [artifact #13047522395](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41214740781/artifacts/13047522395), ZIP and individual JSON/TSV opened and checked.

## 1. Primary erratum full text — RESOLVED

The original October 1, 2026 IOVS published erratum, [PubMed PMID 42895266](https://pubmed.ncbi.nlm.nih.gov/42895266/), is freely archived as [PMC13712379](https://pmc.ncbi.nlm.nih.gov/articles/PMC13712379/) / [Europe PMC official full text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13712379/fullTextXML).

**Verbatim substantive correction category**: *addition of the Data Availability paragraph to the Acknowledgments*. It discloses that the raw RNA-seq was deposited in NCBI SRA under **PRJNA1243143**. The entire short erratum body was programmatically parsed, and the public correction does **NOT** state a modification of DEGs, statistical testing, figures, tables, source cohort counts, or supplementary gene-level effects. The earlier flags marking this erratum's **content** unknown are now resolved.

- FullTextXML SHA256 **`271d8a0c20775d1e58afdd34f06ce5b9d20bcf93a93a40c7109e88f74fdf2594`**.
- Official notice body: three XML paragraph elements, published **October 1, 2026**, DOI **10.1167/iovs.67.12.34**.
- Automated status: `PASS_ERRATUM_DEPOSITION_STATEMENT_VERIFIED`; `notice_only_data_availability_addition=true`; `explicit_numeric_deq_correction_in_this_notice=false`.
- **Correct analytical interpretation:** No reason to recalculate Figure 1/3 based **solely on this published erratum**, because the notice does not correct those values. This is **not** proof that every byte of the published online original supplement remains identical to the previously frozen workbook; direct supplement version/hashes remain distinct provenance checks. Retain the original source workbook hash from earlier P6–P7 work.
- Prior source-audit Actions [#41213817707](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41213817707) correctly retrieved the notice, but was conservative and left correction type unclassified until body inspection. The current Actions #41214240628 resolves the issue.

## 2. Original raw RNA-seq PRJNA1243143: real accession, verified metadata

The new study accession is linked from the **official published correction**, not inferred from a similar-looking GEO study. The official [ENA read_run manifest](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&format=tsv) was downloaded and validated, without accessing actual FASTQ files.

| Audit | Verified result |
|---|---:|
| Project | **PRJNA1243143** |
| Read runs | **46** unique |
| BioSamples | **46** unique accessions |
| Library type | **46 RNA-seq**, **46 paired-end** |
| FASTQ files | **92** (`fastq_ftp` metadata for all runs) |
| File MD5 metadata | **92/92** (hex MD5 values present) |
| Compressed FASTQ total | **89,677,930,411 bytes = 89.68 GB decimal** |
| Submitted FASTQ data downloaded in this audit | **ZERO** |
| Aliased two-run sequencing groups | **23** |
| TM annotated groups | **16** |
| Schlemm canal endothelial annotated groups | **7** |
| Treatment pairs explicitly labeled in sample aliases | **17 / 23** |
| Additional `-1/-2` groups resolved from official BioSample treatment XML | **6 / 6** |
| Vehicle/DEX paired groups after official BioSample XML | **23 / 23 treatment pairs mapped** |
| Timepoint per run (1h, 6h, 48h) | **NOT verified** |
| Unique biological donors/strains per 46-run file universe | **NOT independently verified** |

Official ENA response SHA256 **`8108331d9d9d836486591354503a4767d16f07920dbf6b3ed24487cd06197a9e`**. Stored validated files:
- `MEHROTRA_ERRATUM_ADDENDUM_VERIFIED.json`
- `MEHROTRA_PRJNA1243143_METADATA_READINESS.json`
- `MEHROTRA_PRJNA1243143_ENA_RUN_MANIFEST.tsv`
- `MEHROTRA_PRJNA1243143_PAIRING_QC.tsv`

**Important hierarchy:** There are **23 matched sequencing-library candidate pairs, NOT 23 independent donors**. The same TM strain name recurs in distinct experimental groups. For six pairs originally named only `-1/-2` (TM96, TM151, TM123, SC91, SC89, SC86), **12 official BioSample XML records establish six Vehicle and six 100 nM DEX samples**. Therefore, adding the 17 explicitly labeled pairs, all **23/23 library-pair treatment assignments** have original-source evidence. This does **not** establish exposure duration, independent donor identities, or the 48h QC-eligible subset. The publication describes nine QC-passing TM and four SCE strains.

## 2b. Nominal source cell-strain group consistency

The actual ENA alias stems group into exactly **nine distinct TM prefixes** (`TM92, TM96, TM123, TM134, TM135, TM137, TM140, TM151, TM155`) and **four SCE prefixes** (`SC67, SC86, SC89, SC91`). These *label-level* distinct group counts match the nine QC-passing TM and four SCE strains described in the published study. Several strain labels occur in multiple library pairs (e.g., TM96 and TM151), reinforcing that 23 experiment pairs are **not** 23 independent biological donors. This nominal match alone **does not** prove all deposited experiments passed publication-level QC, link any one run to the 48-hour subset, or establish independent donor genomic identity.

The latest [12-BioSample official XML check #41214740781](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41214740781) resolved treatment attributes for the previously ambiguous `TM96-1/-2`, `TM151-1/-2`, `TM123-1/-2`, `SC91-1/-2`, `SC89-1/-2` and `SC86-1/-2` groups, but none of those 12 BioSample treatment records provided a verified 1h/6h/48h exposure-duration attribute. **Actual source-paper experiment-level metadata is still needed for timepoint mapping.**

## 3. Research ramifications — important priority change

**Previously** the IOVS 100 nM 48h study had been used only from a frozen paper supplementary workbook, and its raw sample-level replicates had not been accessible through the known workflow. **Now** there is an authenticated path to actual source RNA-seq FASTQ. This may enable more convincing independent validation of MYOC/TGFBI/MMP2, batch/donor effects, cell-state compatibility and EMT pathway robustness.

However, direct reprocessing of **~89.7 GB compressed FASTQ** with reference alignment can consume large I/O/storage/CPU and must **not** begin indiscriminately before resolving the sample design.

**Gates for a GO to independent raw-data reanalysis:**

1. **P0 GO:** Official erratum resolved, 46 RNA-seq runs and 92 files with checksums metadata verified. No significant data-value correction reported by notice.
2. **P1 PARTIAL GO / HOLD TIMEPOINT:** confirm a 46-run `BioSample → tissue → strain → run → 1h/6h/48h` crosswalk from original sample metadata or publication. All 23 candidate paired treatment arms are now identified through 17 alias labels and six BioSample XML pairs, but the actual exposure duration, source culture-strain identity and QC-approved 48h subset remain unknown. Do not infer timepoint from `-1/-2`.
3. **P2 HOLD:** isolate the **100nM 48h TM subset** (paper's nine QC-passing TM strains, subject to raw metadata). Determine whether raw gene-count tables are already deposited in an official secondary analysis resource; prefer published, auditable processed raw counts over unneeded 90GB FASTQ realignment.
4. **P3 CONDITIONAL GO:** only after P1/P2 pass, choose low-cost primary validation with independent TM strains and proper donor pairing, negative/positive control panels, and the manuscript's fixed candidates MYOC/TGFBI/MMP2/ACTA2. Use original same-level integer gene counts and a paired contrast, not normalized CPM fitted as count data.
5. **STOP:** treating these 23 groups as 23 independent donors, assuming all are 48h TM, or including SCE in TM evidence without stratification.

## 4. Other blocker unchanged — GSE313736 author 80 gene identities

Official Europe PMC publication searches were executed with `GSE313736`, an exact deposited title, and 74-donor-related terms; the current query returned no matching identified journal article or author-specific 80-gene table. **HOLD** exact overlap with our 83 gene hits. Search failures cannot establish nonpublication, and the correct original 80-gene set must not be synthesized from our own results.

## 5. Thesis wording correction

Prior thesis versions said *Mehrotra correction contents still unknown*: **outdated as of this audit**. Replace with: *“An October 2026 erratum added the raw RNA-sequencing availability statement (SRA PRJNA1243143) to the article; it did not report changes to the published gene-level results. We used the publisher's previously frozen supplementary gene summary and record the original raw SRA manifest as an independent-verification opportunity.”*

The manuscript **Methods and Results numeric conclusions are unchanged** by the erratum notice. The biological-study source claims remain separate: actual TM/SCE sample QC and new raw-data replication are not yet performed.
