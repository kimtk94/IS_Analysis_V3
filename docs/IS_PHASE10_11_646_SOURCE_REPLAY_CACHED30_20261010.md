# IS Phase10/11 — 30-source SNP-level coloc replication and LD numerical QC

**Status: 30/646 PASS, 0 FAIL, 616 source files uncached on server. These are 30 of the original four-locus 646 gene–tissue assays, not 30 independent loci or the expanded 80-region screen.**

## Direct SNP reproducibility

- Recomputed actual coloc.abf H0–H4 using pinned coloc 5.2.3, historical p1=p2=1e-4 and p12=1e-5, input beta/SE/MAF, BBJ original N and first-scalar GTEx QTL N.
- Maximum absolute H0–H4 posterior difference across all 30 PASS assays: 1.69e-15, below 1e-8 threshold.
- The original Google Drive folder contains all 646 source TSVs; 616 means not locally cached on server, not source-file loss.

| Gene | Locus | Original GTEx tissue | n SNPs | Direct PP.H4 |
|---|---|---|---:|---:|
| CALHM2 | BBJ_IS_L002 | Brain Cerebellar Hemisphere | 1583 | 0.7977 |
| FGF5 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.7787 |
| NEURL1 | BBJ_IS_L002 | Brain Cerebellar Hemisphere | 1583 | 0.6776 |
| C4orf22 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.6599 |
| C4orf22 | BBJ_IS_L001 | Brain Cerebellum | 1694 | 0.5597 |
| INA | BBJ_IS_L002 | Brain Anterior cingulate cortex BA24 | 1581 | 0.5364 |
| NEURL1 | BBJ_IS_L002 | Artery Aorta | 1583 | 0.4496 |
| FGF5 | BBJ_IS_L001 | Brain Cerebellum | 1694 | 0.4373 |
| NEURL1 | BBJ_IS_L002 | Brain Cortex | 1583 | 0.3653 |
| SH3PXD2A | BBJ_IS_L002 | Artery Tibial | 1583 | 0.3508 |
| COL4A2 | BBJ_IS_L004 | Brain Putamen basal ganglia | 3015 | 0.3308 |
| COL4A2 | BBJ_IS_L004 | Brain Amygdala | 3016 | 0.3063 |
| COL17A1 | BBJ_IS_L002 | Brain Cerebellum | 1567 | 0.2709 |
| SLK | BBJ_IS_L002 | Brain Spinal cord cervical c-1 | 1583 | 0.2296 |
| CALHM1 | BBJ_IS_L002 | Brain Cerebellar Hemisphere | 1583 | 0.2115 |
| RPL6 | BBJ_IS_L003 | Brain Cerebellar Hemisphere | 2016 | 0.2114 |
| SH3PXD2A | BBJ_IS_L002 | Brain Cerebellar Hemisphere | 1583 | 0.1986 |
| COL4A1 | BBJ_IS_L004 | Brain Amygdala | 3016 | 0.1171 |
| PRDM8 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.0901 |
| C4orf22 | BBJ_IS_L001 | Brain Caudate basal ganglia | 1694 | 0.0818 |
| ANTXR2 | BBJ_IS_L001 | Brain Cerebellar Hemisphere | 1694 | 0.0770 |
| ANTXR2 | BBJ_IS_L001 | Brain Caudate basal ganglia | 1694 | 0.0692 |
| PRDM8 | BBJ_IS_L001 | Artery Coronary | 1694 | 0.0501 |
| ANTXR2 | BBJ_IS_L001 | Artery Coronary | 1694 | 0.0475 |
| FGF5 | BBJ_IS_L001 | Artery Coronary | 1694 | 0.0441 |
| PRDM8 | BBJ_IS_L001 | Artery Aorta | 1694 | 0.0325 |
| PRDM8 | BBJ_IS_L001 | Artery Tibial | 1694 | 0.0256 |
| ANTXR2 | BBJ_IS_L001 | Artery Tibial | 1694 | 0.0019 |
| PRDM8 | BBJ_IS_L001 | Brain Cerebellum | 1694 | 0.0012 |
| ANTXR2 | BBJ_IS_L001 | Artery Aorta | 1694 | 0.0008 |

## 1000G EAS504 LD QC

- 30 original cached gene–tissue files contributed 53,979 source gene–SNP rows; SNPs are repeated across assays.
- L001, L002, L003 and L004 contain four 1000 Genomes EAS504 locus-labelled LD reference matrices; BBJ_IS file names designate loci, not BBJ subject-specific or Japan-only genotypes.
- All cached GRCh37 variant IDs align to the original local REF/ALT-labelled PLINK .pvar and .vars lists.
- The four entire correlation matrices have finite entries, unit diagonals, symmetric elements and pairwise correlations bounded [-1,+1].
- The deterministic 128-variant submatrix eigenvalue checks show no significant negative values, but are not a full-matrix PSD proof.
- This does not independently validate effect-allele direction, a cohort-matched GTEx QTL LD reference, or multi-signal causal mechanisms.

## Source provenance and next release gates

- Original 646-row status file SHA256: 2ef17ff27a092773abc8290f780ccc6953f1dfca1acf7d149000b6a13fa5a47c
- Current LD numeric audit JSON SHA256: 90109a51f051f6245148fc529f7271fd797f80e8816dd63c6d48f6ae759026c6
- Full source SHA256/harmonization counts/LD matrix hashes remain outside Git under results/is/audits/is_646_cache_ld_numeric_20261010_v2.
- Historical cache21 results preserved; current 30-row PASS / 616 MISSING_INPUT status preserved under results/is/audits/legacy_coloc_646_batch_cache30_20261010_v1.
- Next: execute complete Google Colab direct-Drive 646-input replay, review all errors, and validate GRCh37↔GRCh38 allele mapping plus appropriate GTEx QTL LD.
- No stroke causal gene established by numeric reproduction or healthy-reference expression.
