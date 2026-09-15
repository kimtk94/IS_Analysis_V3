# CKD thesis data acquisition status

This project uses public summary statistics for discovery/replication and keeps
restricted KoGES individual-level data outside the public repository.

## Public datasets acquired and CI-verified

- UKB-PPP / Sun et al. 2023 Olink: Supplementary Table 16 independent pQTL signals.
  The canonical exposure set contains **10,750 cis signals**, **1,966 protein IDs**
  and **10,580 unique variant IDs** (GRCh37/hg19).
- EUR eGFRcrea: Stanzick et al. 2021 / CKDGen.
- EUR CKD: Wuttke et al. 2019 / CKDGen.
- EUR BUN: Wuttke et al. 2019 / CKDGen.
- EUR eGFRcys: Gorski et al. / CKDGen.
- EUR UACR: Teumer et al. 2019 / CKDGen.
- EAS eGFRcrea and BUN: Taiwan Biobank + BioBank Japan / Chen et al. 2024.

Every bulk download is gzip-tested where applicable and inventoried with SHA256.

## Analysis-ready design

The first screening stage does **not** download every UKB-PPP per-protein raw
archive. It uses the published independent cis-pQTL signals, matches those rsIDs
to the kidney outcomes, and preserves source alleles/effect estimates. Allele
orientation is deliberately harmonized only in the downstream MR stage.

Current self-contained analysis-ready bundle:

- exposure: UKB-PPP cis instruments
- EUR: eGFRcrea, CKD, BUN, eGFRcys, UACR
- EAS: eGFRcrea, BUN
- per-dataset summary JSON
- MATCH_SUMMARY.json
- SHA256SUMS.txt

Full per-protein pQTL regions should be fetched only for proteins surviving the
initial MR screen and then used for colocalization/fine-mapping.

## KoGES

KoGES is a restricted, approval-dependent Korean individual-level validation
layer. It is **not** downloaded by public CI and must never be committed to this
repository. Planned use is targeted validation of prioritized loci using
baseline eGFR, repeated eGFR/eGFR slope and incident CKD definitions.

## Server execution

From a checked-out repository:

```bash
chmod +x server/ckd_acquire_public_data.sh
CKD_DATA_ROOT=/srv/is-analysis/data/rawdata/ckd \
CKD_READY_ROOT=/srv/is-analysis/data/analysis_ready/ckd \
  server/ckd_acquire_public_data.sh
```

Override the paths if the production mount differs.
