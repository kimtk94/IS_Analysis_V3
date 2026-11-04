# CKD MAF outlier 4-variant provenance audit (2026-10-18)

**Read-only verified run:** [GitHub Actions 39112255591](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39112255591), SUCCESS.

| Gene | rsID | GRCh37 position | pQTL effect / other allele | 1KG REF / ALT | pQTL MAF | EUR ref ALT freq | Ref MAF | Absolute delta |
|---|---|---:|---|---|---:|---:|---:|---:|
| SDCCAG8 | rs35646972 | 242690143 | T / C | C / T | 0.0549205 | 0.167992 | 0.167992 | 0.1130715 |
| SDCCAG8 | rs6663222 | 242990635 | C / T | T / C | 0.0835498 | 0.194831 | 0.194831 | 0.1112812 |
| INHBC | rs11834807 | 58433974 | C / T | T / C | 0.304179 | 0.590457 | 0.409543 | 0.105364 |
| UMOD | rs11860566 | 20208930 | T / A | A / T | 0.00141033 | 0.135189 | 0.135189 | 0.13377867 |

**All four:** Effect/other allele pairs agree with the *derived* 1000 Genomes GRCh37 REF/ALT pairs. These checks do not independently attest original UKB-PPP pQTL or kidney GWAS allele annotations. GRCh37 coordinate alone is not sufficient without chromosome context. Original source builds/alleles, minor-allele definition, reference population frequency, and rsID mapping remain to be revalidated.

**Priority:** UMOD rs11860566 is a much larger relative-frequency discrepancy; independently inspect genotype QC/imputation INFO and original pQTL summary frequency. This does not prove it causes SuSiE nonconvergence. The other three also require source-level auditing.

**No irreversible action:** No variant exclusions, refitting, result overwrites or LIVE changes were performed.

**Next:** Compare the four SNPs against original source pQTL/GWAS allele identifiers and cohort frequencies; then examine whether LD–GWAS Z-score diagnostics point to reference mismatch. Preserve source-matched candidate status as PENDING.
