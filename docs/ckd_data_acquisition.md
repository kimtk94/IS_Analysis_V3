# CKD thesis data acquisition

Data acquisition precedes downstream MR/colocalization workflow development.

## Acquisition order

1. EUR outcomes: CKDGen Stanzick 2021 eGFRcrea and Wuttke 2019 binary CKD.
2. Renal support traits: BUN, eGFRcys, then UACR.
3. EAS outcome: Chen 2024 Taiwan + Japan eGFR meta-GWAS from Figshare article 24356587.
4. Plasma pQTL screening instruments: UKB-PPP cis/independent pQTL signals.
5. Candidate-only UKB-PPP full cis-region summary statistics for colocalization.
6. KoGES individual-level data only after CODA approval; never commit or publicly upload restricted data.

## Why UKB-PPP raw data are not bulk-downloaded first

The existing `ukb_ppp_download_manifest.tsv` contains per-protein archives that are commonly around 0.5 GB each. Bulk-downloading thousands of proteins would therefore create a TB-scale staging requirement. The thesis pipeline first uses compact cis-pQTL/independent-signal instruments to screen proteins by MR. Full regional pQTL data are acquired only for candidates that proceed to colocalization.

## Commands

List core public datasets:

```bash
python3 scripts/download_ckd_public_data.py \
  --dest "$WORK_ROOT/data/rawdata/ckd" --public-core --list
```

Download EUR eGFR and binary CKD:

```bash
python3 scripts/download_ckd_public_data.py \
  --dest "$WORK_ROOT/data/rawdata/ckd" \
  --data-id ckdgen_stanzick2021_egfr_eur \
  --data-id ckdgen_wuttke2019_ckd_eur
```

Resolve the EAS Figshare file inventory before selecting the genome-wide summary-statistics file:

```bash
python3 scripts/download_ckd_public_data.py \
  --dest "$WORK_ROOT/data/rawdata/ckd" \
  --data-id eas_chen2024_egfr_meta --figshare-list
```

Every downloaded file is gzip-tested when applicable and gets a sidecar `.download.json` containing source URL, byte size, SHA256, phenotype, ancestry, genome build and acquisition timestamp.
