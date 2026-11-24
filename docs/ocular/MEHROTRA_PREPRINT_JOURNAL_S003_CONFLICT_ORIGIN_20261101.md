# Mehrotra preprint vs journal Table 1 — verified provenance of TM151/TM155 6-hour discrepancy

**2026-11-01 final official source-level audit**: [GitHub Actions run #41221426210](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41221426210) **SUCCESS**; [artifact #13050700071](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41221426210/artifacts/13050700071) downloaded and JSON inspected (ZIP CRC PASS).

## Three official source versions distinguished

| Source | Record | Relevant exact source-level statement |
|---|---|---|
| Original 2026 **bioRxiv preprint** | [Europe PMC PMC12776393 fullTextXML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12776393/fullTextXML), SHA256 `1c8a4aea9d564868d11262fdd647ec92eb0010c95fe51b8de4e4c0de9ed6b0ba` | Demographic Table 1: **TM151 2 days + 1 hour; TM155 2 days + 6 hours** |
| 2026 Feb **IOVS peer-reviewed article** | [PMC12882097 fullTextXML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/fullTextXML), SHA256 `87a80692f387318037da1d6b7aa45b7a7e546e888abac92a37796e1b274f64c2` | Demographic Table 1: **TM151 2 d + 1 h; TM155 2 d + 6 h** |
| Official 2026 IOVS **Supplement Table 82** | `iovs-67-2-13_s003.xlsx`, SHA256 `3c30e93461812d05f42a1978b8e446235709f2283e11f620c212beb26478698d`, from [official source ZIP](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/supplementaryFiles) | Control TPM columns for 6h **TM151**, *not TM155*. Thus supplement conflicts with BOTH preprint and peer-reviewed Table 1. |

**Corrected comparison**: the preprint spells `2 days, 1 hour` and `2 days, 6 hours`; the final article spells `2 d, 1 h` and `2 d, 6 h`. A naive string comparison falsely reported differences; normalized canonical time units demonstrate **NO biological duration differences** across preprint versus journal for either TM151 or TM155. The original first CI #41220863273 and intermediate #41221050316 were conservatively HOLD due source caption formatting and textual duration normalization; final CI #41221426210 resolves these issues, preserving failed/review history transparently.

## Evidence-based conclusion

The TM151-vs-TM155 **6-hour source disagreement predates the final journal article**: a preprint-era original Table 1 already specifies TM155 at 6h, as the journal does. The 2026 article's October erratum added only Data Availability (SRA PRJNA1243143) and did not address the table/supplement contradiction. There is **no basis to decide whether the demographic table or the supplied expression column label is correct**. Author sample sheet and original aliquot metadata are needed.

**48-hour TM is NOT affected in the named strain set**: the nine 2-day TM strains agree across preprint Table 1, journal Table 1, supplement S003 and ENA nominal strain aliases. The nine SRA pairs are still **candidate**, however, because each actual RUN lacks a direct 48-hour treatment attribute and raw source QC eligibility is not separately confirmed. Do not treat a multi-source match from one study as independent biological validation.

## Next narrow deliverable

1. Obtain original authors' `RUN / BioSample / culture strain / treatment / 1h,6h,48h / extraction batch / QC included` mapping; ask whether the 6h source strain in Table 82 should be TM151 or TM155.
2. Once this source crosswalk is available, decide whether sample-level original integer gene counts exist to validate selected **nine 2-day TM** pairs. Do **not** download 32.7GB/89.7GB FASTQ before that gate.
3. Original GSE313736 authors' **80-member** high-effect gene list is a different missing external artifact and remains **HOLD**.

**All evidence in this document is publication-source and metadata QC. No new RNA count-level model or functional glaucoma phenotype is claimed.**
