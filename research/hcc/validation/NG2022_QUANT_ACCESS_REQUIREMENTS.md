# HCC P0 — Ng 2022 quantitative RNA×protein access and verification contract

Updated: 2026-10-31
Decision owner: study governance; **access NOT approved and quantitative patient join NOT established**.

## Grounded dataset and source access

- Paper: Ng et al., *Nature Communications* (2022), DOI [10.1038/s41467-022-29960-8](https://doi.org/10.1038/s41467-022-29960-8).
- Paper cohort: 122 HCC tumor biopsies / 114 patients. RNA-seq performed on 122 tumors; proteome and phosphoproteome on a selected 51 tumors. **51 tumors are NOT 51 unique patients**. The public Supplementary Data 1 audit identified 49 distinct patients flagged for all three assays, but this does NOT establish matching gene-by-sample and protein-by-sample quantitative values.
- RNA-seq: [EGA EGAS00001005074](https://ega-archive.org/studies/EGAS00001005074) — source explicitly states restricted access; do not try automated bypass or equate landing page with expression matrix.
- Public proteome: [PRIDE PXD025705](https://www.ebi.ac.uk/pride/archive/projects/PXD025705); public phosphoproteome [PXD025836](https://www.ebi.ac.uk/pride/archive/projects/PXD025836). A public PRIDE *file list* is not proof of accessible processed quantitative matrix.
- Public supplementary data 1–14 / Source Data accompany the publication. A separate CI probe screens file layout and limited safe header families without storing row-level sample IDs.

## Quantitative data to request or verify

1. **Per-biopsy RNA normalized gene expression matrix**: rows stable Ensembl/Entrez gene identifiers, columns unique biopsy sample IDs, explicit count/TPM/RSEM measure, transformation and annotation version, exact 51 tumor biopsy intersection. Confirm whether a processed matrix can lawfully be shared independently of restricted raw sequencing.
2. **Per-biopsy SWATH protein-abundance matrix** for the same tumor biopsies: protein accession/gene mapping, normalized intensity/ratio definition, peptide aggregation, missingness, batch/instrument correction and quality-control exclusions.
3. **Approved, de-identified cross-modality lookup table**: unique patient key, unique biopsy key, RNA column, proteome column, phosphoproteome column, biopsy/tissue labels, duplicate patient/biopsy handling. Store mapping in an approved protected workspace; publish only aggregate overlap counts.
4. **Clinical variables**: source-defined nonexclusive ALD/HCV/HBV/NAFLD etiologies, stage, grade, treatment; define mixed-etiology rule *a priori*. Never assume the current 13 ALD-only and 11 HCV-only modality-flag counts are quantitative matched analysis cohorts.
5. **License/security/authorization**: institutional DUA/DAC approval, exact legal sharing scope, storage security, permitted analyses and permitted result disclosures. Do not upload restricted expression, patient-level clinical data or individual sample IDs into public GitHub Actions logs or artifacts.

## Testable release gates

| Gate | Evidence required | Current decision |
| --- | --- | --- |
| A. Public supplement enumeration | Download checksum + safe worksheet shape, audited filenames | IN_PROGRESS |
| B. RNA gene×biopsy matrix | Real numeric values, gene IDs, source checksum, quantitative scale, authorized access | BLOCKED |
| C. Protein×biopsy matrix | Real normalized quantitative values, protein IDs, documented batch | BLOCKED |
| D. RNA-protein same-biopsy join | Unique biological patient/biopsy keys in both matrices, one-to-one/replicate audit, actual joined n | BLOCKED |
| E. Etiology contrast | Verified exclusive and mixed subtype membership within matched n, covariate balance, minimum power | HOLD |
| F. Prognosis | Event/status codebook, time origin, censoring, sufficiently independent outcomes | HOLD |

## Interpretation and next action

Complete technical metadata discovery and public supplementary geometry screening through GitHub Actions, then request lawful processed RNA expression and normalized SWATH abundance/lookup from the corresponding investigator or through EGA DAC as applicable. The original paper requests contacting the corresponding author for restricted sequence access and says the author responds within four weeks. This is a **source-stated process**, not a commitment by any connected service or proof that access will be granted.

**Do not promote V02 to PASS** until both numeric modalities and same-biopsy linking are directly verified. Keep V04 PARTIAL and V12 PARTIAL until clinical definitions are documented. Avoid simple proteogenomic quadrant/pathway claims that the original paper already reported.
