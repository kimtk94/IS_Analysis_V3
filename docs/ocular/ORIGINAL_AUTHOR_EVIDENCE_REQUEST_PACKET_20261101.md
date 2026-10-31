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
