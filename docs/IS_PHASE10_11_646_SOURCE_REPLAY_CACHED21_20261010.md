# IS Phase10/11 original 646 coloc assays — local cache21 and LD numerical audit

**Snapshot: 2026-10-10/11 KST. Source-only numeric reproducibility; not an independent biological replication.**

## 646-source input statuses

- Original source assays: **646** from four original BBJ GWAS loci; this does not represent the expanded 80 provisional loci or 2,225 candidate gene IDs.
- Direct SNP-level R coloc.abf 5.2.3 source replays **21 PASS / 0 FAIL / 625 MISSING_INPUT**.
- The 625 MISSING_INPUT files are absent from this server's local cache; the owner's original Google Drive folder inventory independently contains 646 TSV files.
- Maximum H0–H4 posterior error across 21 PASS inputs: **1.69e-15** (numeric tolerance 1e-8).
- An existing source posterior being computationally reproducible does **not** validate allele harmonization, ancestry reference LD, gene mediation or causal disease mechanisms.

## Source-replayed gene–tissue tests (all PASS)

| Gene | GWAS region | Historical GTEx tissue | SNPs | PP.H4 | Max |
|---|---|---|---:|---:|---:|
| CALHM2 | BBJ_IS_L002 | Brain Cerebellar Hemisphere | 1583 | 0.7977 | 1.7e-15 |
| FGF5 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.7787 | 4.7e-16 |
| NEURL1 | BBJ_IS_L002 | Brain Cerebellar Hemisphere | 1583 | 0.6776 | 1.2e-15 |
| C4orf22 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.6599 | 7.2e-16 |
| INA | BBJ_IS_L002 | Brain Anterior cingulate cortex BA24 | 1581 | 0.5364 | 8.3e-16 |
| FGF5 | BBJ_IS_L001 | Brain Cerebellum | 1694 | 0.4373 | 1.3e-15 |
| SH3PXD2A | BBJ_IS_L002 | Artery Tibial | 1583 | 0.3508 | 8.9e-16 |
| COL4A2 | BBJ_IS_L004 | Brain Putamen basal ganglia | 3015 | 0.3308 | 9.4e-16 |
| COL4A1 | BBJ_IS_L004 | Brain Amygdala | 3016 | 0.1171 | 3.3e-16 |
| PRDM8 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.0901 | 1.9e-16 |
| C4orf22 | BBJ_IS_L001 | Brain Caudate basal ganglia | 1694 | 0.0818 | 1.1e-16 |
| ANTXR2 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.0770 | 4.7e-16 |
| ANTXR2 | BBJ_IS_L001 | Brain Caudate basal ganglia | 1694 | 0.0692 | 2.5e-16 |
| PRDM8 | BBJ_IS_L001 | Artery Coronary | 1694 | 0.0501 | 3.6e-16 |
| ANTXR2 | BBJ_IS_L001 | Artery Coronary | 1694 | 0.0475 | 5.3e-16 |
| FGF5 | BBJ_IS_L001 | Artery Coronary | 1694 | 0.0441 | 9.7e-16 |
| PRDM8 | BBJ_IS_L001 | Artery Aorta | 1694 | 0.0325 | 3.3e-16 |
| PRDM8 | BBJ_IS_L001 | Artery Tibial | 1694 | 0.0256 | 4.4e-16 |
| ANTXR2 | BBJ_IS_L001 | Artery Tibial | 1694 | 0.0019 | 1.1e-16 |
| PRDM8 | BBJ_IS_L001 | Brain Cerebellum | 1694 | 0.0012 | 1.3e-16 |
| ANTXR2 | BBJ_IS_L001 | Artery Aorta | 1694 | 0.0008 | 1.1e-16 |

## 1000G EAS504 matrix and GRCh37 allele-ID numerical QC

- **21** local gene–tissue TSV source inputs, **37,771** source gene–SNP records, across **3** reference matrices (GWAS region labels L001, L002, L004). Repeated SNPs across gene–tissue inputs are counted repeatedly in this denominator.
- LD provenance: PLINK 2.0 genotype reference **1000 Genomes EAS with 504 reference samples**. A BBJ_IS_* filename denotes regional label, not BBJ participants or JPT-only genotype LD.
- All cached SNPs match the source reference by GRCh37 chromosome:position:REF:ALT ID. GRCh38 eQTL match keys are different and require a verified build conversion.
- The binary correlation matrices are finite and symmetric, unit diagonal, have correlations bounded by [-1,+1]. Deterministic submatrix eigenvalue QC uses 128 spread-out matrix indices; this is a **sampled PSD diagnostic only**, not a full-matrix eigenvalue guarantee.
- The local PLINK reference allele ID and input SNP ID agreement does not independently validate effect allele alignment of BBJ GWAS versus GTEx eQTL.
- This is **GWAS-side EAS reference**. It is not donor/cohort-matched GTEx molecular QTL LD; multi-signal joint coloc SuSiE remains blocked pending additional evidence.

## Provenance / status

- Original 646-row replay status SHA256: `0500d22529f9d60426bc0764cdea5c430fb23362a99178f0338a1065efca1a4c`
- Source 21-file LD audit SHA256: `d2aa482b99cd4204875383d7235ba2541b6ca2e030593c74e4b2e1cdd231e293`
- Full source SHA256 list, 21 gene–tissue audits and per-reference LD source hashes are stored outside Git at results/is/audits/is_646_cache_ld_numeric_20261010_v1/.
- Original replay matrix/status output is stored outside Git at results/is/audits/legacy_coloc_646_batch_cache21_20261010_v1/.

## Next experimental gate

Run the one-cell 646-original-source Colab against the complete Google Drive dataset. Review all 646 PASS/FAIL outcomes and preserve the full status file. The current 21/646 PASS should not be presented as the entire historical screen being independently rerun.

Following complete computational replay, validate variant/build/effect alleles and obtain appropriate QTL-LD and independent molecular evidence. Preserve the broader candidate discovery universe.
