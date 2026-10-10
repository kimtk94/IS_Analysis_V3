# IS G0022 — regional EAS AIS GWAS × Japanese OASIS cis-eQTL overlap

**Audit:** 2026-10-10 KST (output also saved October 11 near midnight)
**Claim gate:** SOURCE_SNP_POSITION_REF_ALT_MATCHED; QTL_EFFECT_ALLELE_UNVERIFIED; WHOLE_CIS_UNIVERSE_UNVERIFIED; NO_VALID_COLOC.

## Original data and build transformation

GIGASTROKE EAS AIS source GCST90104545 for G0022 has **5,876 original GRCh37 variant rows**, **5,789 quality-controlled records** and **87 excluded by the existing source QC**. All 5,789 passed GRCh37→GRCh38 *coordinate* mapping using local official UCSC hg19ToHg38.over.chain.gz forward primary chromosome-12 chain ID **11** (score 12,341,034,620; 153 blocks). Each of the four known AIS 95% credible-set anchors mapped to its prior GRCh38 coordinates. This uses aligned blocks, not a guessed constant offset, and does not independently verify genome REF bases against an assembly FASTA for every SNP.

The previously archived Japan Omics OASIS Bokeh gene-level pages were joined by mapped GRCh38 position and **EXACT REF/ALT nucleotide strings**, not simple distance, correlation or surname. Multibase SNPs and any allele mismatches are not used as matched SNVs. Source GWAS ALT-normalized effect/se/p is retained but **OASIS beta effect-allele orientation is unknown**, so no signed harmonization or effects regression was performed.

## Results

| Gene | Cell type | OASIS plotted variants | Precisely matched EAS AIS SNV records | Browser "maf" values >0.5 |
|---|---|---:|---:|---:|
| ALDH2 | Mono-L1 | 2,254 | 1,736 | 1,082 |
| ALDH2 | B_Activated-L2 | 2,266 | 1,747 | 1,080 |
| BRAP | Mono-L1 | 2,140 | 1,632 | 1,048 |
| BRAP | B_Activated-L2 | 2,162 | 1,651 | 1,045 |
| RPH3A | Mono-L1 | 2,968 | 2,155 | 1,217 |
| RPH3A | B_Activated-L2 | No embedded Bokeh source | NOT ASSESSED | NOT ASSESSED |
| **Total** | Five populated source pages | **11,790** | **8,921 gene-cell-SNP records** | **5,472** |

**The 8,921 are not independent SNPs!** There are **2,866 distinct source GWAS variant IDs** in their union. The same variant can recur across several gene-cell results, and LD makes marginal association records strongly dependent. The 5,789 lifted GWAS variants are not all in each OASIS gene-level cis interval.

The three other top G0022 95%-CS SNPs rs671, rs11066132 and rs77768175 are **still unobserved on the five Bokeh pages**, even though within the respective displayed cis ranges; nearest displayed SNP gaps were 214, 211 and 35bp. This does not imply absent from original E-GEAD-1054 data or biologically negative.

## Critical allele-frequency mislabelling

The OASIS Bokeh JSON array is named **maf**, but **5,472/11,790 underlying values exceed 0.5**, so this field cannot possibly be a validated *minor allele frequency*. Examples include 0.808, 0.885 and 0.917. The official nominal matrix from the first E-GEAD-1054 internal file instead names the field **af**. This mismatch must be resolved with source documentation/validation before computing ALT QTL effect frequencies, using a p-value colocalization prior, or harmonizing beta signs.

The code labels this source field **oasis_browser_maf_column_AS_REPORTED_NOT_TRUE_MAF** and keeps OASIS effects UNHARMONIZED. Merely converting this column to min(af,1-af) is NOT sufficient to infer which allele the beta refers to.

## Scientifically valid conclusion and next barrier

There is now a large **molecular and stroke source-overlapping locus**, not yet evidence of shared causal variants:
- ALDH2 Monocyte-L1: 1,736 matched SNVs, rs11066015 eQTL p=0.000913.
- BRAP Monocyte-L1: 1,632 matched SNVs, rs11066015 p=0.303.
- RPH3A Monocyte-L1: 2,155 matched SNVs, rs11066015 p=7.62e-7.

The comparison is molecular-QTL marginal association annotation. Larger numbers and smaller p-values are not causal gene rankings, independent MR instruments or independent Japanese stroke replication.

Publication-grade colocalization requires *complete tested cis SNP–gene sets*, alleles/effect frequency semantics, QTL credible sets/conditional signals, per-cell N, matched signed LD and overlap/prior sensitivity. Current browser coverage is only the displayed Bokeh series; completeness is unverified. Distinct study cohorts can have different AF/QTL power. No valid coloc posterior or causal gene was computed.

Keep ALDH2, BRAP, RPH3A, ACAD10, NAA25, HECTD4 and the broader **2,225-gene IS candidate discovery set** intact.

## Descriptive SVG figure

**G0022_OASIS_EAS_AIS_SHARED_SNP_MONO_REGIONAL_FIGURE.svg** plots distinct GWAS AIS and OASIS monocyte -log10(p) panels for ALDH2, BRAP and RPH3A using only matched source SNVs. Dashed rs11066015 and rs671 mark the known AIS candidates; rs671 is missing from the OASIS plotted SNP list. Caption explicitly states **no causal or coloc claim**. Self-contained SVG validated with Python XML parser. No graphic library installation required.

## Reproduction, all input data read-only

New Git scripts:
- scripts/is/audit_g0022_oasis_regional_gwas_overlap.py
- scripts/is/render_g0022_oasis_regional_overlap_svg.py
- tests/test_is_g0022_oasis_regional_gwas_overlap.py
- tests/test_is_g0022_oasis_regional_figure.py

Unchanged original inputs:
- /srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_browser_gene_v1/

New output directory:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/oasis_gwas_fullregional_overlap_v1/

Outputs:
- G0022_OASIS_EAS_AIS_REGIONAL_ALLELE_MATCHED_SNVS.tsv (8,921 linked gene-cell rows)
- G0022_OASIS_EAS_AIS_SOURCE_COVERAGE.tsv (six source pages, SHA256 and coverage)
- G0022_OASIS_EAS_AIS_REGIONAL_OVERLAP_SUMMARY.json (full source and science gating)
- G0022_OASIS_EAS_AIS_SHARED_SNP_MONO_REGIONAL_FIGURE.svg

Official UCSC chain: https://hgdownload.gi.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
Original OASIS study: https://www.nature.com/articles/s41588-025-02266-3
Original E-GEAD-1054 source: https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-1000/E-GEAD-1054/

**Original GWAS, QTL HTML, prior canonical research records, existing LD and unrelated production services unchanged. No 39GB bulk download performed.**
