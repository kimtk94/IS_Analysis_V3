# Ocular DEX/TM — minimal original-author evidence request packet

**Prepared:** 2026-11-01. **Status:** documentation prepared only. **No messages sent or external permissions requested.**

## A. GSE313736 original investigators: actual author 80 DEG IDs and 59/57 QC

**Original author provenance:** [2024 Australian Vision Research report, page 20](https://australianvisionresearch.org/wp-content/uploads/2025/08/2024_AVR-Annual-Report.pdf) attributes the single-cell mapping work to PI Alex W. Hewitt; [2024 University of Tasmania showcase abstract, page 34](https://mooc-media.utas.edu.au/files/7/5/171877818175/CoHM%20Research%20Showcase%202024%20-%20Presentation%20Abstracts%20Booklet.pdf) lists presenter Thomas Guinan and collaborating authors. Official [GSE313736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736) identifies its original submission contact at the University of Tasmania (the official GEO record provides an address; avoid relying on guessed contact details).

**Short specific requests:**
- Actual `63 upregulated + 17 downregulated` original Ensembl gene IDs and symbols, `log2FC`, raw `p`, BH `FDR` and comparison/annotation; *no raw donor genotype or personal identifiers needed*.
- The exact **decreased-gene cutoff** underlying the 17: AVR 2024 report incorrectly/ambiguously prints `>3` for the decreases and public GEO prints `<3`, so the actual intended log2FC sign and exact inequality are not authenticated.
- Original single-cell/genotype QC donor count/table distinguishing the 2024 **59 QC-passing samples** from the original 74 donors and our paired **57 BeamA-selected donors**. Request anonymized aggregate QC reason groups, not identifying per-person sensitive records.
- Preprint/manuscript/DOI, differential-expression model, cell-state definition, donor as replicate and source 10x library batch covariates, and whether a table or public GitHub source code already exists.

**What will qualify as GO:** author-supplied/published gene table with stable locus identifiers and reproducible contrast/threshold, and evidence of its provenance. A single narrative number 63+17 or a list inferred from our own 83 is **NOT sufficient**.

## B. Mehrotra et al. 2026: 1h/6h/2d run/timepoint/QC crosswalk and internal inconsistency

**Official peer-reviewed source:** [IOVS 2026 PMC12882097](https://pmc.ncbi.nlm.nih.gov/articles/PMC12882097/), with erratum [PMC13712379](https://pmc.ncbi.nlm.nih.gov/articles/PMC13712379/) revealing raw SRA [PRJNA1243143](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1243143). The journal paper lists corresponding authors and their contact information in the published article; use **the publisher-provided correspondence address**, not a guessed personal email.

**Short specific requests:**
- Original 46 `SRA RUN → BioSample → TM/SCE → strain → vehicle or 100 nM DEX → 1h / 6h / 2d` mapping, plus culture batch and **QC-included** samples (9 passing 2-day TM strains versus source screening samples).
- Clarify whether **TM151** or **TM155** was actually studied in the **six-hour TM** arm. Both [Jan 2026 preprint Table 1](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12776393/fullTextXML) and Feb 2026 [journal Table 1](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML) state **TM155 at six hours**; official [Supplementary Table 82](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/supplementaryFiles) instead has TM151 six-hour marker columns. A sample sheet is stronger than a narrative correction.
- Ask whether a published/downloadable **complete per-sample integer featureCounts matrix** (preferable to 89.7GB raw FASTQ reprocessing) exists; if yes, supply source URI/checksum/reference version and document QC scope.
- Clarify why published SCE six-hour source groups appear in the expression supplement yet no six-hour SCE sample pair is identifiable in the 46-run PRJNA1243143 subset; avoid asserting a data omission without original deposit clarification.

**GO for independent count-level 48h TM:** source-verified exact 9-sample paired strain design with 18 run IDs, treatment/time and QC inclusion; then raw counts. **Before GO:** retain the candidate 32.71GB source subset as REVIEW, no speculative FASTQ download.

## Contact/outreach log (to be filled after actual outreach)

| Study | Evidence request drafted | Sent? | Official reply or artifact | New GO/HOLD decision |
|---|---|---|---|---|
| GSE313736 original author group | YES | **NO** | None | HOLD author 80 IDs / exact original QC |
| Mehrotra IOVS 2026 original authors | YES | **NO** | None | HOLD authoritative 1h/6h/2d sample crosswalk |

No individual human-genotype/medical record data should be requested. Do not claim contact has occurred until an actual user-authorized communication and reply are verified.


## C. Ready-to-send inquiry to GSE313736 source investigators (UNSENT)

**Official route (verified 2026-11-02):** The [GSE313736-linked original GEO GSM record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9373290) lists University of Tasmania / Menzies Institute contact **Alex Hewitt**, institutional email **alex.hewitt@utas.edu.au** (also lists a personal email; prefer institutional). Address is sourced from the live official NCBI GEO record, not inferred from other publication authors. **Outreach remains UNSENT.** Do not automatically send from CI.

**Suggested subject:** Request for original 80 DEX-responsive TM genes and methods (GSE313736)

> Dear Professor Hewitt and research team,
>
> I am conducting an independent, methods-focused reanalysis of the public GSE313736 trabecular meshwork single-cell dataset. I located your group's 2024 Australian Vision Research report (63 DEX-upregulated and 17 downregulated genes), the July 2024 University of Tasmania showcase abstract (59 samples passing genotype and single-cell QC), and the 16 February 2024 ANZGS conference program featuring Thomas Guinan's presentation.
>
> Would you be able to share, or point me to a published supplementary source containing, the original **80-gene DEG results** (Ensembl gene IDs, symbols, log2 fold changes, nominal P values and adjusted P values)? The report and GEO abstract differ in the stated downregulated log2 fold-change cutoff, so a clarification of the actual threshold, DE contrast, cell-type subset, statistical model and unit of replication would be extremely helpful.
>
> If available, I would also appreciate the public analysis code or package versions and **aggregate, de-identified QC criteria/counts** explaining the original 59 passing samples. I am **not requesting individual genotypes or donor-identifying clinical records**.
>
> I would cite the original work and any relevant preprint, manuscript or published DOI appropriately. Thank you for considering this reproducibility request.
>
> Kind regards,
> [Researcher name / university affiliation]

**Evidence receipt checklist (not assumed received):** (1) author-owned 80-row gene table with stable IDs; (2) directional cutoff and effect scale; (3) original sample/cell-state eligibility and contrast/model; (4) file or publication provenance and version. **Do not** promote overlap from HOLD to GO on the basis of the reply text alone unless the actual gene-level table is attached or linked and checked. **Sent:** NO. **Reply:** NONE.


## D. Submission-ready minimal request — verified contact, not yet sent (2026-11-02)

- To: `alex.hewitt@utas.edu.au` (verified from [NCBI GEO GSM9373290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9373290)).
- Subject: `GSE313736 — request for original 63 up / 17 down DEG table and analysis methods`.
- Requested artifact: exact author-owned 80-gene IDs, gene symbols, signed log2FC, adjusted p-values, contrast and 59-sample QC criteria; optionally R analysis script and paper/preprint DOI.
- Correction to draft language: **do not claim** ANZGS conference talk in February 2024 disclosed the DEG 80 — its official program only provides a scheduled title and time. The 2024 AVR report provides the aggregate 63/17 numbers.
- Source search status: targeted public indexing and official GEO Series table listing checked. No authenticated original 80-row table retrieved. Parallel baseline P5 57 donor and Gate B matched 55 donor remain separate reanalysis, not an author-specific reproduction.
- Approval/sent/replied timestamps: **none**. Do not write to an author automatically from GitHub Actions or include genotype/clinical identifiers.
