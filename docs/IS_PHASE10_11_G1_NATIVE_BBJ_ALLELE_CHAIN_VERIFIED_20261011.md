# G1 original BBJ GWAS effect Allele2 and GRCh37→GRCh38 independent chain verification

**2026-10-11; direct read of original Japanese ischemic-stroke GWAS native archive; 8,309 of 8,309 distinct locally cached GWAS SNPs PASS. This does not constitute biological causality or full-646 eQTL source validation.**

## Evidence sources and interpretation

- Native BBJ original autosomal file: hum0197.v3.BBJ.IS.v1.zip (contains GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz). Independently streamed **13,435,541 native GWAS rows**.
- Original fields: v, CHR, POS, Allele1, Allele2, AF_Allele2, BETA; Allele2 is the source GWAS effect allele, with BETA and AF_Allele2.
- Original coloc-cache input: 30 gene–tissue files, 53,979 gene–SNP records, **8,309 distinct BBJ GWAS GRCh37 SNP IDs**. Duplicate gene–tissue appearances 45,670. All 8,309 unique IDs were found in the native ZIP.
- Independent hg19ToHg38 chain coordinates compared against GTEx GRCh38 match_key; each mapping had to be unique with chrom/position/strand and REF/ALT (reverse complement where needed).
- The BBJ native effect Allele2 is GRCh37 ALT for all 8,309 SNPs. Native BETA and AF_Allele2 match archived coloc GWAS beta/EAF within 1e-9 tolerance.
- All 8,309 chain mappings uniquely on plus strand; QTL GRCh38 REF/ALT strings match the mapped GRCh37 REF/ALT.
- GTEx Portal official FAQ confirms GTEx eQTL effect (NES/slope) is ALT relative to REF. However this does not independently attest the exact GTEx file transform, actual GRCh38 reference FASTA bases, or GTEx genotype LD.
- Liftover chain confirms coordinates and REF/ALT strings, not actual reference sequence bases; this proof does not establish gene mediation.

## Results by original GWAS region

| Original BBJ GWAS region | Distinct GWAS SNPs | Native beta/EAF proof | Unique chain/GTEx alleles |
|---|---:|---|---|
| BBJ_IS_L001 | 1,694 | 1,694/1,694 PASS | 1,694/1,694 PASS |
| BBJ_IS_L002 | 1,583 | 1,583/1,583 PASS | 1,583/1,583 PASS |
| BBJ_IS_L003 | 2,016 | 2,016/2,016 PASS | 2,016/2,016 PASS |
| BBJ_IS_L004 | 3,016 | 3,016/3,016 PASS | 3,016/3,016 PASS |

**Total: 8,309 PASS, 0 mismatch, 0 native variants absent, 0 ambiguous chain mappings.**

## Scientific validation boundaries

- **G0 original numeric coloc: COMPLETE 646/646**, unchanged.
- **G1 BBJ native effect allele and chain: COMPLETE for 8,309 distinct GWAS SNPs** represented by local 30 gene–tissue files; not evidence that all 646 GTEx source files are independently revalidated.
- G1 remaining: independently verify any additional SNP identities in 616 eQTL source tests, actual GRCh38 reference FASTA bases, official GTEx source transformation and any exceptional strand/allele patterns.
- G2 N sensitivity is complete only for 30 locally cached gene–tissue tests; full-646 follow-up Colab prepared but not executed. Appropriate GTEx-QTL LD and multi-signal inference pending.
- G3 disease-relevant regulatory validation pending. Four original BBJ loci are not the expanded 80 provisional genomic windows.

## Provenance

- GTEx official FAQ: https://www.gtexportal.org/home/faq (eQTL ALT effect allele).
- Original native ZIP size: 1,611,973,884 bytes; chain SHA256: 5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d.
- Verified 8,309-SNP TSV SHA256: 963d5632624c90209e92e7ae22351b2db6ce9199bea5c68cb1a037f2a1be6dfe.
- Full source audit files outside Git: results/is/audits/is_g1_native_bbj_chain_20261011_v1.
- Script: scripts/is/audit_is_g1_bbj_native_allele_chain.py; synthetic tests tests/test_is_g1_native_bbj_chain.py. Original inputs unmodified.
