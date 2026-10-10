# IS G0022/ALDH2 — Japanese ImmuNexUT molecular source audit
**Date:** 2026-10-10 KST
**Scientific classification:** NEW SUPPORTED_VARIANT_LEVEL_MARGINAL_EQTL; VALIDATED_COLOCALIZATION = 0; CAUSAL_MEDIATION = NOT ESTABLISHED.

## 1. Source release gates

| Dataset | Japanese/EAS sample scope | Source nominal QTL selection | Publication gate |
|---|---|---|---|
| NBDC NHA000172 JCTF v2 (2022) | 465 Japanese COVID patients | cis eQTL/sQTL records with p below 0.05 only | BLOCKED_FULL_QTL |
| NBDC NHA000193 JCTF v3 (2024) | Japanese COVID; 1019 eQTL, 1384 Olink pQTL | p below 0.05 OR PIP above 0.001 | BLOCKED_FULL_QTL |
| DDBJ E-GEAD-420 ImmuNexUT (2021) | 416 Japanese immune patients and healthy donors, 28 cell types | Complete nominal eQTL including non-significant pairs claimed officially | PROMISING_FULL_QTL_SOURCE, LOCAL_VARIANT_COVERAGE_UNVERIFIED |
| MAEEA Zenodo 21296030 (2026) | 2024 Chinese/Japanese donors | Complete SNP–gene resource described; previous server audit restricted; 2026-10-10 Zenodo API timed out | ACCESS_BLOCKED_LIVE_UNVERIFIED |

Official Japanese NBDC release descriptions establish that both JCTF v2 (749 MB) and v3 (1.7 GB) have restricted publication selection. Neither can be treated as an unbiased all-tested SNP denominator for new coloc.abf / coloc.susie from these released tables alone. No JCTF ZIP bulk download was performed in this phase.

## 2. ImmuNexUT full nominal archive: bounded original HTTP inspection

Official DDBJ E-GEAD-420.processed.zip:
- HEAD 200, content-length 43,507,433,191 bytes, Accept-Ranges bytes.
- First 4,096 bytes and final 262,144 bytes retrieved using HTTP 206 byte-range requests.
- ZIP first local entry and central end entry: nominal_eQTL_2.tar.gz. ZIP compression method is 8, ZIP64 end-record detected; official DDBJ file list corroborates the single tar.gz archive member.
- No 43.5GB archive was downloaded. A single doubly-compressed TAR cannot support ordinary region queries from the ZIP central directory: a targeted rs671-only query is not currently possible without a different indexed distribution or bulk extraction.

This technical access issue does NOT negate completeness of the nominal source metadata; it means the actual variant-specific contents still must be checked. The Japanese per-person genotype archive is not openly downloadable.

## 3. Direct 2026-10-10 ImmuNexUT ALDH2 eQTL browser verification

The official ALDH2 significant-QTL gene page was downloaded as HTML (source SHA256 in audit JSON) and parsed into 1,000 displayed records. The official UI shows only significant associations, with top 1,000 limitation for highly associated genes. Exact rsID plus GRCh38 coordinate and Neu cell-type verified all four GWAS AIS G0022 credible-set variants:

| AIS GRCh37 variant | Japanese eQTL rsID | GRCh38 position chr12 | cell | browser nominal p | browser beta (UNHARMONIZED) |
|---|---|---:|---|---:|---:|
| 12:112168009:G:A | rs11066015 | 111730205 | Neu | 4.88802e-9 | -0.207695 |
| 12:112241766:G:A | rs671 | 111803962 | Neu | 4.88802e-9 | -0.207695 |
| 12:112468206:C:T | rs11066132 | 112030402 | Neu | 4.88802e-9 | -0.207695 |
| 12:112736118:A:G | rs77768175 | 112298314 | Neu | 4.88802e-9 | -0.207695 |

Neu = neutrophils (per official site abbreviation table).

**Interpretation:** All four EAS credible-set SNPs have significant marginal Japanese ALDH2 neutrophil eQTL observations in a second molecular resource besides JCTF. Four identical p/beta reflect high LD and cannot separate causal SNPs. Source effect allele identity was not extracted; negative QTL beta cannot be aligned with AIS beta direction. The molecular cohorts are from different resources, but independent participant membership was not established. The ALDH2 eQTL could occur in other cell types; the web UI top SNP/cell summary is not a comprehensive cell-type null screen.

No valid GWAS–QTL coloc, causal gene or alcohol/BP/protein mediation was estimated or validated. ALDH2 catalytic loss, alcohol-related behaviors, QTL expression and Olink assay-binding artifact must remain distinct hypotheses.

## 4. Machine-readable files and reproducibility

Git worktree: /srv/is-analysis/worktrees/is-broad-discovery-20261009
- scripts/is/audit_is_eas_molecular_qtl_source_gates.py
- scripts/is/audit_is_immunexut_aldh2_browser.py
- tests/test_is_eas_qtl_access_and_immunexut.py

Read-only verified HTML:
- /srv/is-analysis/data/is/qtl/immunexut/G0022_v1/ALDH2_significant_eqtl_browser_20261010.html

Research outputs:
- results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/qtl_source_gates_20261010/IS_G0022_EAS_QTL_SOURCE_GATES.json
- results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/immunexut_aldh2_v1/G0022_IMMUNEXUT_ALDH2_FOUR_CS_NEU_EQTL.tsv
- results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/immunexut_aldh2_v1/G0022_IMMUNEXUT_ALDH2_NEU_SUMMARY.json

## 5. Next strict G0022 gate

1. Ask DDBJ/ImmuNexUT for a region-indexed BGZF or chromosome/cell-type extract of E-GEAD-420 nominal eQTL (including non-significant) before attempting a 43.5GB transfer.
2. Once data access is obtained, count the actual tested SNP–ALDH2 gene pairs per cell type and verify all four GWAS CS alleles, QTL ALT orientation, beta/SE/p, SNP-specific sample sizes, and imputation.
3. Match Japanese or GIGASTROKE EAS LD with cohort ancestry; EAS504/JPT104 remain low-rank exploratory references.
4. Study phenotypic and cohort overlaps. Distinguish coding enzyme/behavior effects from expression-mediated pathways and from Olink protein epitope effects.
5. Ask the MAEEA custodians for full nominal cis summary dataset and permissible access; NARD2 imputation service does not automatically yield raw signed haplotype LD.
6. Until all gates hold, retain 2,225-gene universe and all G0022 competing genes; no claim of ALDH2 mediation or drug-target directionality.

## Official sources
- JCTF v2: https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000172
- JCTF v3: https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000193
- ImmuNexUT source and 1000 display rule: https://www.immunexut.org/faqs
- ImmuNexUT ALDH2 actual marginal variants: https://www.immunexut.org/eqtlGenes?gene_symbol=ALDH2
- ImmuNexUT Neu label: https://www.immunexut.org/abbreviations
- NBDC E-GEAD-420 release: https://humandbs-production.ddbj.nig.ac.jp/en/research/hum0214/v7
- DDBJ file manifest: https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/E-GEAD-420.filelist.txt
- MAEEA: https://academic.oup.com/nsr/advance-article/doi/10.1093/nsr/nwag566/8779978

**No canonical IS GWAS, original prior QTL results, non-IS projects or active services were modified.**
