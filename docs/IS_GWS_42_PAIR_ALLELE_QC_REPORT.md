# IS Broad Discovery — seven GWS regions, GWAS/ref overlap audit

Research snapshot: 2026-10-09 (KST).
Study comparison: BBJ ischemic stroke + 5 GIGASTROKE EAS sources (AS, AIS, CES, LAS, SVS).
References: 1000G EAS 504-sample, GRCh37 regional panels; G0004/G0008/G0021 new stable-ID PGEN, G0007/G0015/G0022/G0023 historical BBJ PGEN.

## Actual data-scan results

The complete canonical GWAS inputs were scanned, not sampled.

- Source rows scanned: BBJ 13,435,531; EAS GIGASTROKE AS 6,789,909, AIS 6,725,216, CES 5,169,600, LAS 5,993,005, SVS 6,173,637.
- Seven GWS groups × six source GWAS = **42** GWAS/group combinations.
- **145,949** GWAS/window variant *records* across source/group combinations (not unique variants).
- **107,697** records with matching reference alleles and computable orientation (73.79%).
- **1,465** palindromic SNP records unresolved (1.00%).
- **36,787** SNP records whose position is absent from the regional reference (25.21%).
- Explicit source REF/ALT conflict **0** records.
- Lead allele matched for all **42/42** GWAS/group pairs.
- Lead variant ALT EAF delta versus EAS 1KG AF was ≤0.10 for **42/42**, none >0.10, none missing.

Aggregate coverage by GWS group:

| Group | Chromosome | GWAS/window records | Reference/orientation matched | Palindromic review | Position absent |
|:--|:--|--:|--:|--:|--:|
| G0004 | chr2 | 11,817 | 10,922 | 91 | 804 |
| G0007 | chr4 | 13,430 | 12,080 | 110 | 1,240 |
| G0008 | chr4 | 17,927 | 16,397 | 358 | 1,172 |
| G0015 | chr10 | 29,404 | 13,541 | 173 | 15,690 |
| G0021 | chr12 | 13,016 | 11,710 | 62 | 1,244 |
| G0022 | chr12 | 37,589 | 21,661 | 371 | 15,557 |
| G0023 | chr13 | 22,766 | 21,386 | 300 | 1,080 |

**Priority flag:** L002/G0015 chr10 and L003/G0022 chr12 regional references cover much less of the *expanded* discovery window. The original four regional PGENs were designed around narrower BBJ loci; they are not a full expansion to 2+ and 4+ Mb merged regions. Re-extract those larger windows before comprehensive credible set claims.

## Reproduction

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
python3 scripts/is/audit_gws_gwas_reference_overlap.py
python3 scripts/is/build_gws_study_readiness.py
python3 scripts/is/audit_gws_lead_eaf.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Output root:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/`

- `IS_GWS_GWAS_REFERENCE_OVERLAP.tsv`
- `IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz`
- `IS_GWS_STUDY_INPUT_READINESS.tsv` and summary
- `IS_GWS_LEAD_EAF_QC.tsv` and summary
- `reference_frequency_qc/*.afreq`

No reference-aligned genotype LD covariance matrices or SuSiE credible sets are claimed. Matching lead alleles and EAF are necessary source-QC checks, not sufficient proof of valid colocalization.

## Remaining gates

1. Complete and MD5 verify EUR AIS GCST90104540 raw data. Its download is separately running and not considered validated until the expected source checksum matches.
2. Normalize EUR effect/other alleles without treating them as genomic REF/ALT. Audit build/allele and cohort overlap before comparison.
3. Expand 1000G EAS regional reference windows at G0015 and G0022. Also assess coverage of all seven groups by actual variant-position overlap rather than reference boundary only.
4. Compute matched study-specific signed LD covariance and signal-level fine-mapping; low LD sample size (504) and reference ancestry differences should be analyzed in sensitivity checks.
5. Screen eQTL, sQTL, pQTL and expression/cell-type evidence for all 869 candidate gene-region rows; NOT_TESTED is not a negative finding.
