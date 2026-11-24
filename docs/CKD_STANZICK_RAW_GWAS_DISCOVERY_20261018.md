# CKD Stanzick 2021 original eGFRcrea file — verified location (2026-10-18)

Actual successful read-only server inventory: [GitHub Actions 39123205308](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39123205308).

Source: `data/metadata/ckd_public_download_manifest.tsv`
- Manifest key: `ckdgen_stanzick2021_egfr_eur`
- Source: `https://ckdgen.imbi.uni-freiburg.de/files/Stanzick2021/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz`
- Declared genome build: `GRCh37`

Verified existing full original file on research server:
`/srv/is-analysis/data/ckd/rawdata/outcome/EUR/ckdgen_stanzick2021_egfr_eur/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz`, 290,933,654 bytes.

Distinguish from the restricted subset `data/ckd/analysis_ready/EUR/eGFRcrea.tsv.gz` (164,733 bytes), which contains none of the four target rsIDs. The presence of the raw GWAS is verified; original GWAS allele and beta for each target must be verified by an independent read from this file. Do **not** treat location confirmation as SNP-specific validation.

Raw GWAS script: `scripts/audit_ckd_stanzick_four_variants.py` (read-only, target rsIDs and GRCh37 chr:pos). GitHub Actions first scan [39123437003](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39123437003), optimized rerun [39123850576](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39123850576); execution/conclusion requires inspection before interpreting allele-level results.

Original UKB-PPP pQTL source-side records for four MAF outliers were previously verified. Study-matched LD, GWAS source records, and UMOD SuSiE explanation remain separate science gates.
