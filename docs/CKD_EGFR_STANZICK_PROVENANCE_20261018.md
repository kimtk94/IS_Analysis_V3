# CKD eGFRcrea original source provenance check — 2026-10-18

Verified execution: [GitHub Actions 39122942790](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39122942790).

## Findings from existing analysis-ready provenance JSON

- source schema: `stanzick_egfr`; phenotype `eGFRcrea`; ancestry `EUR`
- original source rows reportedly parsed: **8,844,847**
- independent instrument rsIDs: **10,580**
- GWAS matched rsIDs in the selected instrument lookup: **7,245** (68.478%)
- duplicate matched rows: **0**
- provenance notes source-published allele orientations not yet harmonized at the recorded stage.

These are **pipeline provenance counters**, not a fresh independent parse of all 8.8 million original GWAS rows. No full raw EUR eGFRcrea file was identified in the read-only directory inventory so far. The four outlier rsIDs were not present in the limited analysis-ready eGFRcrea subset.

## Prior confirmed UKB-PPP pQTL results

All four target records were recovered directly from original per-protein UKB-PPP archives in [original pQTL audit](./CKD_TOP3_FOUR_OUTLIERS_ORIGINAL_PQTL_VERIFIED_20261018.md). In those archives, `GENPOS` and the coordinate portion of `ID` diverge at nearly all scanned loci. The evidence shows dual coordinate representations; **genome-build attribution remains to be confirmed with original study documentation** rather than inferred from these values alone.

## Scientific verification boundary

- Original UKB-PPP pQTL allele and frequency on four targets: **SOURCE_VERIFIED**
- eGFRcrea original study schema: **PROVENANCE_IDENTIFIED**
- eGFRcrea original four target records: **PENDING**
- Source-provided GWAS allele orientation vs pQTL orientation: **PENDING**
- 1KG EUR external LD vs original study LD equivalence: **NOT ESTABLISHED**
- UMOD SuSiE nonconvergence cause: **UNRESOLVED**

No original research data or live settings were changed.

Next direct technical work: locate original Stanzick eGFRcrea GWAS download source from ingestion script or download manifest and target only the four variants, then independently verify eGFRcrea allele/effect estimates. Do not substitute eGFRcys.
