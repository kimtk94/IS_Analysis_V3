# CKD four MAF outliers — direct UKB-PPP raw archive audit (2026-10-18)

**Verified run:** [GitHub Actions 39122708873](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39122708873), SUCCESS. This runs entirely read-only against the original local UKB-PPP full-locus TAR files and scans both source `GENPOS` and coordinate embedded in `ID`.

| Gene | rsID | Source GENPOS | Source ID variant position | ALLELE0 / ALLELE1 | Raw A1FREQ | Source MAF | Raw INFO | Raw BETA | Raw SE | N |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| SDCCAG8 | rs35646972 | 242526841 | 242690143 | C / T | 0.0549205 | 0.0549205 | 0.804870 | -0.00896149 | 0.0179946 | 33777 |
| SDCCAG8 | rs6663222 | 242827333 | 242990635 | T / C | 0.0835498 | 0.0835498 | 0.691110 | -0.0330412 | 0.0159850 | 33777 |
| INHBC | rs11834807 | 58040191 | 58433974 | T / C | 0.695821 | 0.304179 | 0.975895 | +0.102113 | 0.00822668 | 33687 |
| UMOD | rs11860566 | 20197608 | 20208930 | A / T | 0.00141033 | 0.00141033 | 0.764296 | +0.216093 | 0.109696 | 34030 |

All four original source frequencies produce the previously observed harmonized pQTL MAFs. The large MAF differences from external 1KG EUR are *already present in original source data* and cannot be attributed simply to a Stage2B/Stage2C frequency-transformation bug.

**Critical format finding:** Source `ID` includes coordinate and alleles (e.g. `16:20208930:A:T:imp:v1`), distinct from `GENPOS` (e.g. `20197608`). Nearly all inspected rows have different values across these columns. This is consistent with differing build coordinates, but the precise GRCh37 versus GRCh38 field definitions **have not been independently proven by source documentation**. Variant matching was performed by ID-internal coordinate plus allele pair; it did not substitute external liftover. The previous direct GENPOS/rsID lookup returned zero because rsID does not appear as an ID and coordinates did not match.

**Interpretation:** UMOD `rs11860566` original pQTL A1FREQ=0.00141033 while 1KG EUR minor allele frequency=0.135189. The original pQTL INFO=0.764296 is not missing, but INFO alone does not resolve population-frequency discrepancy. Determine whether build mapping, variant identifiers, technical imputation, multi-allelic handling, ancestry or reference-panel processing explain divergence, before inferring a SuSiE cause. INHBC A1FREQ=0.695821 corresponds to MAF=0.304179, not 0.695821.

**Next independent evidence gate:** Verify the alternate GWAS consortium original eGFRcrea source (not summarized/harmonized and not eGFRcys), confirm builds with original study metadata, and inspect nearby variant signed Z/LD compatibility. No variant was dropped, reweighted or refitted; no server source files changed.
