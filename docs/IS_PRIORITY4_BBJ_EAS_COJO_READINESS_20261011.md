# IS priority-four — BBJ IS × EAS504 GCTA-COJO input readiness

**Audited 2026-10-11 KST.** This audit builds **nonproduction, exploratory input files only**; no `gcta64 --cojo-cond` or `--cojo-slct` was run. The source data and previous coloc posteriors remain unchanged.

## Scientific and software gate (BLOCKED)

GCTA's primary [official COJO reference manual](https://yanglab.westlake.edu.cn/software/gcta/#COJO) states:
- LD reference genotype sample **more than 4,000 individuals is suggested**;
- **1000 Genomes and HapMap panels are explicitly not recommended** as COJO reference samples because their LD estimates have substantial sampling error;
- to calculate phenotypic variance, GCTA asks for **the GWAS summary data of all SNPs**, even when COJO tests are restricted to particular regions (use `--extract` to focus the analysis);
- for a binary case-control trait the effect `b` should be **log(odds ratio)** and `se` the corresponding standard error;
- `--diff-freq` can exclude GWAS/reference allele-frequency mismatches; a 0.2 absolute frequency difference is the documented default.

**Current server:** `gcta64` is not on PATH. Existing Phase 3 **EAS504 is far below the >4,000-person reference recommendation**. Thus a programmatically executable EAS504 COJO run must **not be presented as a validated independent-signal result**. The current audit does not launch the GCTA binary.

## What was actually verified

We scanned **all 13,435,531 records** in the frozen canonical BBJ IS GWAS `BBJ_IS_GRCh37.canonical.tsv.gz`; every potential SNP was aligned by exact GRCh37 `chr:pos:REF:ALT` to EAS504 genotypes previously fetched directly from the 1000 Genomes Phase 3 regional VCF. This starts with **regional BBJ GWAS data**, not only the narrower EBI QTL SNP intersection.

Filters and orientation: source `INFO >= 0.8`; GWAS `MAF >= 0.01`; valid numerical beta, SE, p, allele frequency and **source per-SNP N**; biallelic SNP; excluded REF/ALT palindromes (A/T and C/G) and multiallelic records; BBJ `effect_allele` explicitly aligned to genomic ALT and effect β flipped only when effect allele was REF; frequency `freq` oriented to the same ALT; Bim ID and alleles reconfirmed after PLINK `--make-bed`.

| Numeric QC | chr4 FGF5/C4orf22 | chr10 CALHM2/NEURL |
|---|---:|---:|
| Exact BBJ regional interval (GRCh37) | 80,087,872–82,356,697 | 104,113,502–106,353,104 |
| BBJ original GWAS rows in interval | 8,514 | 11,979 |
| Exact EAS reference SNP + GWAS INFO/MAF QC preliminary | 3,513 | 3,759 |
| SNPs surviving PLINK genotype extraction + EAS MAF>=0.01 | 3,374 | 3,448 |
| EAS-vs-BBJ absolute ALT frequency delta >0.20 additionally removed | 0 | 1 |
| **COJO-format allele and frequency QC-ready SNPs** | **3,374** | **3,447** |
| EAS sample count in PLINK BED/FAM | 504 | 504 |
| BBJ reported GWAS per-variant N | 174,686 | 174,686 |
| Absolute median EAS-vs-BBJ ALT frequency delta | 0.03447 | 0.03423 |
| Reference LD publication-readiness gate | **BLOCKED** | **BLOCKED** |

The provisional GCTA `.ma` columns are exactly `SNP A1 A2 freq b se p N`; genomic ALT is A1, genomic REF A2. Coordinates and effects are **not** derived from rsID-only matching. The frozen SNP outputs are not the full genome-wide BBJ summary required by the official GCTA variance-estimation procedure.

### Distinguish total N from potentially effective case/control N

The frozen BBJ source column `n` consistently contains **174,686** throughout these loci, and that raw source value is preserved unchanged in `.ma`. Because ischemic stroke is binary, **case/control design and the intended COJO sample-size convention must be independently reconciled** using study metadata and method-specific documentation. It is scientifically incorrect to silently replace the source `N` with a guessed effective N. This audit therefore marks **case/control N handling UNVERIFIED**, not a successful full scientific readiness state.

### Prepared research-server artifacts

- `scripts/is/prepare_is_priority4_bbj_eas_cojo_readiness.py` — full-source BBJ scan and exact EAS GRCh37 SNP join, preliminary `.ma`, `extract.snplist`, `lead.snplist`, manifest.
- `scripts/is/finalize_is_priority4_bbj_eas_cojo_reference_qc.py` — recomputes actual EAS504 PLINK `ALT_FREQS`, verifies SNP alleles against BIM/FAM, genotype sample size, per-variant frequencies; emits the final **format-QC** `.ma` and guarded manifest.
- Preliminary files: `/srv/is-analysis/results/is/audits/is_priority4_bbj_eas_cojo_readiness_20261011_v1/`
- Reference PLINK BED/BIM/FAM and allele-freq: `/srv/is-analysis/results/is/audits/is_priority4_bbj_eas_cojo_reference_20261011_v1/`
- Final nonproduction `.ma` with manifest: `/srv/is-analysis/results/is/audits/is_priority4_bbj_eas_cojo_final_20261011_v1/`

No downstream biological claim has changed. Original **2,225 positional IS candidate genes**, archived **646 coloc** runs, and separate **120 EBI–BBJ ABF** sensitivity scenarios are preserved.

## Next legitimate COJO step

Find or obtain a large, ancestry-matched **BBJ cohort / Japanese genotype LD panel** with enough participants (official GCTA guidance suggests >4,000), verify effective sample-size convention and complete BBJ summary statistics, then construct an appropriate COJO reference and rerun allele/frequency/source regression QC. Before that, the EUR503, EAS504 and JPT104 work provides **reference-specific LD diagnostic evidence** only, not valid independent conditional stroke signals. For QTL fine-mapping the separately missing GTEx molecular donor cohort LD is still an independent blocker.

## Identified larger Japanese LD reference (future controlled-access option)

**BioBank Japan official published inventory** (as of 2025-12-01): the 1000G-augmented BBJ imputation reference panel includes **7,472 Japanese BBJ WGS samples** under **NBDC JGAS000738 / JGAD000873**. The official published inventory also lists 15,646 BBJ WGS participants and 269,000 SNP-array individuals, with different levels of access. Dataset https://biobankjp.org/researchers/730 . A population-specific 7K BBJ reference-panel study is https://pubmed.ncbi.nlm.nih.gov/39702642/ .

The 7,472-person panel **exceeds GCTA's suggested >4,000 reference sample guideline**, but the raw variant-level genotypes are **controlled-access**: approval and permissible reuse must be checked. No claim of authorization or file possession is made. This option is explicitly **DEFERRED** to avoid blocking public-data-only thesis progress. It would still need overlap-with-GWAS assessment, variant/QC harmonization, effective case-control sample-size agreement, and appropriate genome-wide BBJ summary input before formal COJO.

Frequency panels such as jMorp 8.3KJPN / 54KJPN provide allele frequencies but do **not by themselves** supply the signed, participant-level LD correlation matrix required for this analysis. Source https://jmorp.megabank.tohoku.ac.jp/docs/guide-en/overview/access_to_datasets/ .

**Current runnable public-data task ends at reproducible COJO-format input QC**; no COJO command was executed. No controlled BBJ genotype access was attempted.
