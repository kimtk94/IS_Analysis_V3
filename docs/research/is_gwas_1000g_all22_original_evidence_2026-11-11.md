# Original GIGASTROKE EUR/EAS IS GWAS × ancestry-matched 1000G Phase3 LD reference — 22-autosome source verification

2026-11-11. **Internal reproducibility and source QC ledger, not course submission text.**

## Original source and 22 autosome execution evidence

- GIGASTROKE EUR IS `GCST90104540`, EAS IS `GCST90104545`; genome assembly GRCh37 in original GWAS filenames and published phenotype, same original files whose publisher-MD5 and gzip CRC were previously validated.
- Ancestry-specific LD variant reference from [S-LDSC 1000 Genomes Phase 3 Zenodo 7768714](https://zenodo.org/records/7768714), EUR `1000G_Phase3_plinkfiles.tgz`, EAS `1000G_Phase3_EAS_plinkfiles.tgz`.
- **Official compressed source publisher-MD5 matches in the full 22-chromosome run**. The 1000G BIM source build was not independently coordinate-proven across a genome reference and requires explicit verification before variant-level fine-mapping.
- First attempt [#43031126709](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43031126709): EUR **PASS**, EAS failed due to read timeout before source completion; *not* a biological or variant-matching disagreement.
- Resumable retry implementation [#43031787049](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43031787049): **EUR and EAS both SUCCESS**, both full original LD and GWAS source MD5s verified, every autosomal BIM and GWAS row staged and classified; previously validated chromosome 1+2 row and allele-bin counts reproduced exactly in both ancestries.

## Actual full-autosome source results

| QC category | EUR | EAS |
|---|---:|---:|
| Original GWAS rows (22 autosomes) | 7,482,032 | 6,725,216 |
| Original 1000G Phase3 BIM rows | 9,997,231 | 8,768,561 |
| **Coordinate in ancestry-specific BIM reference** | **98.7091%** | **97.6406%** |
| **Nonpalindromic unordered A/C/G/T allele-pair exact** | **83.8924%** | **82.7498%** |
| Palindromic A/T or C/G variants with effect direction ambiguous | 1,108,581 | 1,001,439 |
| GWAS coordinate absent from LD BIM | 96,584 | 158,676 |
| Publisher checksum GWAS and LD reference | PASS both | PASS both |
| Prior independently successful chr1+2 regression | PASS | PASS |
| Full association effect-direction harmonization | NOT_RUN | NOT_RUN |
| Actual LD correlation matrix construction | NOT_RUN | NOT_RUN |
| Suitable EAS complete cis-eQTL variants and effect statistics | HOLD | HOLD |
| GWAS–eQTL coloc / MAGMA / LDSC-SEG | NOT_RUN | NOT_RUN |

**The reference coverage difference is technical.** It does not imply differences in ischemic-stroke pathophysiology.

### Chromosomal technical flags

- **Chr5 outlier**: EUR 96.6929% and EAS 95.8492% coordinates in population-matched LD reference, notably below the usual ~98.6–99% EUR / 97.7–98.3% EAS. Check specific source variant classes and any GRCh37 reference mapping issues before analyzing chr5 loci. A genome-wide average masks regional/missingness problems.
- **EAS chr6**: 96.4614%, lower than many other EAS chromosomes. MHC/complex variation may be a factor; this explanation is untested and must not be asserted as cause without a localized match audit.
- **Chr16** has relatively low nonpalindromic match fraction even when coordinate overlap is high: EUR 81.2720%, EAS 80.2466%. This measure includes ambiguous palindromic exclusion; need bin breakdown before interpretation.

## Correct technical interpretation

1. BIM A1/A2 allele sets **do not align** the GWAS effect allele and beta sign. A nonpalindromic unordered allele-set match is a candidate for a subsequent harmonization step, not proof of orientational correctness.
2. Palindromic A/T or C/G remain **UNRESOLVED** without ancestry-appropriate effect-allele frequencies and strand/reference information; do not arbitrarily flip or delete without a predetermined rule.
3. SNPs not present in the 1000G LD panel are unevaluable with this specific LD reference. No biological conclusions from unequal ancestry-specific missingness.
4. GIGASTROKE original standardized GWAS lacks per-SNP N and the contribution of clinical vs ICD-based stroke ascertainment by ancestry is not independently fully audited; gene mapping by eQTL must control power/ascertainment before comparison.
5. 1000G reference PLINK BED/FAM can be used to compute actual local LD after reference build, population sample ancestry and allele direction have been audited. **No actual LD or fine-mapping has been computed here.**
6. MAEEA EAS whole-blood eQTL `Zenodo 21296030` previously had zero public attached files; Japanese JCTF QTL is significance/PIP selected, not a complete cis distribution. Do not use these to claim coloc completeness. Japanese NHA000051 n105 data and EUR BLUEPRINT full associations have only **remote-access/object checks** so far, not locus-wide QTL QC.

## Next run requirements

- Reconcile both original source-artifact files via `scripts/is/aggregate_gwas_1000g_all22.py` into a sample- and patient-free tabular comparison. The synthesis must fail if one ancestry's source MD5, all-chromosome status, chr1+2 regression, or row/bin totals disagree.
- Then run locus-based positional reference-build verification and **effect/eQTL allele orientation**. Explicitly record unresolvable palindromic variants and use effect-allele frequency only when a well-specified reference frame exists.
- Prioritize EAS **full cis-eQTL** and balanced EUR comparator; until evidence is available, locus-to-gene cross-population concordance is **NOT_EVALUABLE**, not discordant biology.
- MAGMA gene/property is optional/secondary; LDSC-SEG is conditional on ancestry-specific LD scores and heritability QC. Neither is independent replication.

**GitHub main PR:** [Draft #281](https://github.com/kimtk94/IS_Analysis_V3/pull/281), not merged.
