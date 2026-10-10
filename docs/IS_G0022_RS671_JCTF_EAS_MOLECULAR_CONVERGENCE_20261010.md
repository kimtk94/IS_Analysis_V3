# IS G0022: rs671-centered East Asian QTL convergence, coverage audit and mediation research plan
**Research date:** 2026-10-10 KST
**Research branch:** `research/is-broad-discovery-20261009`
**Verdict:** **NEW REAL East-Asian molecular-QTL variant-level evidence**; no validated AIS–QTL colocalization, no established stroke causal gene/mediated pathway.

## Research context

Current broad discovery remains **80 ancestry/study provisional intervals**, **2,225 unique positional Ensembl gene IDs**. G0022 is a chr12 >4Mb locus. GIGASTROKE East-Asian AIS `GCST90104545`, signed LD from 1000 Genomes EAS **504 samples**, whole-locus SuSiE-RSS, 4,649 SNPs, assumed study-level case-control effective N≈70,474. Baseline model converged and produced a single, purity-filtered exploratory 95% credible set of four variants.

This does not confer causality: n_ref=504 and 4,649 variant LD has rank ≤503; no fully verified per-SNP effective N / finite-LD correction / independent stroke replication.

## 1. Four GWAS credible-set variants converge on one highly correlated haplotype

| GRCh37 SNP | Known name/positional annotation | Exploratory AIS SuSiE PIP | 1KG EAS LD r² to rs671 |
|---|---|---:|---:|
| `12:112241766:G:A` | **ALDH2 rs671, coding missense** | ~0.37484 | 1.000000 |
| `12:112168009:G:A` | ACAD10 rs11066015 (positional) | ~0.27247 | 0.987528 |
| `12:112468206:C:T` | NAA25 positional overlap | ~0.29776 | 0.974237 |
| `12:112736118:A:G` | HECTD4 positional overlap | ~0.03778 | 0.869614 |

Biological literature confirms that rs671 is a coding change reducing ALDH2 enzyme catalytic function, particularly relevant to East Asian populations. **This is a separate and stronger functional category than cis-eQTL tagging**, but direct enzyme effect is not evidence for ALDH2–ischemic-stroke causality.

Using `z_j|k=(z_j-r_jk*z_k)/sqrt(1-r_jk²)` and reference-signed EAS LD after conditioning on rs671, **0 of 4,648 eligible remaining G0022 variants** retain GWAS significance at p<5e-8. Largest approximate conditional |Z|=4.3293 (p≈1.496e-5). This is *only a single-lead external-LD approximation*, not a formal cohort-level conditional GWAS or proof of one causal signal. The same 504-person reference can overfit high-LD residuals.

## 2. Why GTEx and six other public QTL datasets did not overlap

Original official GRCh38 eQTL Catalogue nominal source extracts from four GTEx tissues — Artery Aorta, Artery Coronary, Artery Tibial, Brain Cortex — contributed 110,739 association rows and 14,370 gene–SNP matching records. Direct source coordinate audit 4 GWS CS SNP×4 tissues: **0/16 source positions present**.

Extended to six *other* separate public eQTL Catalogue datasets using conservatively rate-limited 720kb chr12 region queries:

| Dataset | Study | Tissue / cell type | Catalogue study samples |
|---|---|---|---:|
| QTD000609 | OneK1K | CD14 monocytes | 959 |
| QTD000620 | OneK1K | NK cells | 981 |
| QTD000434 | ROSMAP | brain | 560 |
| QTD000051 | BrainSeq | brain | 479 |
| QTD000021 | BLUEPRINT | monocytes | 191 |
| QTD000110 | GEUVADIS | lymphoblastoid cell lines | 445 |

All 24 additional variant×dataset source-position checks failed to find the four GWAS CS variant positions. These studies yielded **16,350 allele-matched gene–SNP records** for the prespecified eight nearby genes, but not any GWAS CS variant itself.

Across **10 available QTL datasets, 40/40 exact source-position tests absent**. This reflects dataset-specific variant selection, frequency and source cohort composition; it does NOT prove absence of genetic regulation or a negative colocalization. OneK1K is a European-ancestry cohort as described by the source, so it does not solve rs671 ancestry ascertainment.

This is a strong reason to **stop repeatedly analyzing the same non-EAS QTL source subsets** and instead prioritize an East-Asian-specific QTL resource.

## 3. NEW Japanese molecular results: JCTF Japan Omics Browser

Used the official publicly browsable Japan Omics Browser `https://japan-omics.jp/` (Japan COVID-19 Task Force, published by Wang and colleagues in Nature Genetics 2024; Japan Omics Browser tool article in BMC Genomics 2025).

All **four GRCh37 GWAS CS variant identities** were translated to exact GRCh38 genomic positions with official UCSC reciprocal liftOver, and looked up in JOB by rsID or GRCh38 variant. The source pages are cached *as provided* (HTML), exact query identity checked, raw SHA256 recorded and all 72 eQTL/pQTL gene associations extracted into machine-readable TSV (without bulk source reconstruction).

For ALDH2, all four variants have **real Japanese blood eQTL and measured protein pQTL associations**:

| G0022 GWAS CS SNP | ALDH2 eQTL association p | Japanese eQTL SuSiE variant PIP | ALDH2 pQTL association p | Japanese pQTL SuSiE variant PIP |
|---|---:|---:|---:|---:|
| ACAD10 rs11066015 | 6.66e-12 | 0.1560 | 3.18e-48 | 0.0225 |
| **ALDH2 rs671** | **3.51e-12** | **0.3540** | **8.42e-49** | **0.0926** |
| NAA25-position SNP | 8.79e-12 | 0.1280 | 6.59e-49 | 0.1400 |
| HECTD4-position SNP | 3.42e-11 | 0.0196 | 5.37e-48 | 0.0190 |

The original JCTF 2024 cohort contained **1,019 eQTL RNA-seq and 1,384 pQTL Olink Explore participants**, both from Japanese COVID-19 patients; these are **study-level cohort numbers** and do not guarantee per-SNP effective sample sizes. Browser `effect size` alleles are not explicitly re-harmonized to AIS signed GWAS effect alleles, so we do not interpret directionality between traits.

**Critical interpretation:**
1. Variant-level co-association of rs671 and its LD partners with AIS, ALDH2 transcript and Olink-measured ALDH2 protein is supported as **three separate marginal genetic association observations**. It does not establish all are influenced by the *same causal variant* at the signal level.
2. eQTL association p-values (and pQTL p-values) cannot substitute for `coloc.susie` posterior or fine-mapped 95% CS overlap when the complete tested SNP universe and reference LD differ.
3. ALDH2 rs671 is a nonsynonymous protein mutation. A measured pQTL at a coding variant can reflect genuine abundance change, protein conformation or affinity/epitope of the assay; require orthogonal LC-MS quantification before using as a mediator.
4. Another strong cis-target in JOB is RPH3A eQTL at rs671 (p≈2.52e-20), but rs671 QTL SuSiE PIP=0 in the JOB table. This illustrates **strong association and true fine-mapped causal membership are not equivalent**. Do not select RPH3A as causal purely from marginal p-value.
5. Japan COVID Task Force patient selection may induce context-specific QTL effects; not the same sampling as healthy stroke controls.

## 4. Source download/access limitations explicitly established

**Japanese QTL NBDC public release**:
- NBDC **NHA000193 / hum0343.v3.qtl.v1**, unrestricted, ZIP 1,696,958,414 bytes (~1.7GB), with original `eqtl_sumstats_...` and `pqtl_sumstats_...` files. Official NBDC README inspected directly, specifying GRCh37 and GRCh38 variant IDs, per-row `slope` of ALT effect allele, nominal p, SuSiE and FINEMAP PIP.
- **CRITICAL MISSINGNESS:** Release README says rows are included only when **p<0.05 OR PIP>0.001**, and PIPs below 0.001 set to 0. Therefore the downloadable file is **not an exhaustive nominal cis-QTL summary SNP universe**. This ascertainment invalidates naïve full-summary `coloc.abf`/SuSiE without accounting for unreported non-significant alleles.
- The downloaded 2KB README appears to retain an older sample-count comment (465), whereas the 2024 source paper and NBDC v3 metadata specify 1,019 eQTL and 1,384 pQTL. Use the **cohort/assay-specific original paper metadata**, not the stale generic README count.
- The source ZIP is not downloaded by default because it is 1.7GB and still filtered; variant-level records required for current milestone are already available and have been archived from JOB.

**2026 East Asian meta-eQTL**:
- Wang et al., *National Science Review*, published 2026-09-02, reports MAEEA **N=2,024** Japanese/Chinese combined whole-blood eQTL.
- Zenodo DOI `10.5281/zenodo.21296030` publicly describes full cis SNP–gene association statistics; direct REST metadata verified `access_right=restricted`, `files=[]`. **Do not claim public downloadable file access.** Researcher permission/contact or a separately released accessible archive is needed.
- This is higher-priority for future ancestry-matched replication if complete source data becomes available.

## 5. Scientific claims matrix

| Claim | Status |
|---|---|
| Genome-wide East Asian AIS chr12 association | OBSERVED_GWS |
| One full-locus exploratory 4-SNP GWAS credible set | EXPLORATORY_FINE_MAPPING |
| rs671 experimentally known ALDH2 functional missense | EXTERNAL_FUNCTIONAL_EVIDENCE |
| rs671/other CS SNPs associated with Japanese blood ALDH2 mRNA and measured protein | **SUPPORTED_VARIANT_ASSOCIATION** |
| Same causal GWAS variant acts through ALDH2 expression/protein | **NOT_ESTABLISHED** |
| Direct enzymatic vs regulatory vs blood-pressure/alcohol mediation | **NOT_DISENTANGLED** |
| Independent East-Asian stroke cohort replication | NOT_ESTABLISHED |
| Valid multi-signal GWAS–eQTL–pQTL colocalization | **NOT_PERFORMED** |
| Established causal gene or approved therapeutic target | **NOT_ESTABLISHED** |

This work improves candidate interpretation **without narrowing the original genome-wide study to 9 genes**. Maintain the full 2,225 positional-gene universe and track HECTD4, NAA25, RPH3A, PTPN11 and other nearby genes as competing explanations.

## 6. Follow-up investigation priorities

1. **Source completeness:** Request comprehensive tested variant cis-eQTL/pQTL summary and allele QC from JCTF / MAEEA to enable actual signal-level `coloc.susie`, or use alternative EAS public dataset with complete non-significant SNP coverage. Source pages alone are insufficient.
2. **Protein orthogonality:** Compare Olink binding-dependent ALDH2 protein associations to LC-MS/MS measurements or orthogonal immunoassays; missense rs671 can alter assay affinity.
3. **Stroke mediation/confounding:** Contrast ALDH2 catalysis→acetaldehyde metabolism and alcohol consumption (with sex/ethnicity differences), blood pressure, lipid effects and direct vascular pathways. Using rs671 as an alcohol-consumption IV without exclusion-restriction sensitivity is unsafe because rs671 affects acetaldehyde metabolism independently of intake. Test separate conditional/multivariable MR plus pleiotropy/Egger where instruments allow.
4. **GWAS/LD:** Obtain larger EAS LD reference and study-specific `N_eff`; use proper finite-reference LD uncertainty model. Current one-lead conditional calculation is exploratory.
5. **Replication & tissue:** Independent EAS AIS cases if accessible, cerebral arterial/vessel wall / microglial/vascular single-cell tissue context, and enhancer/coding biochemical support.
6. **Existing broader evidence:** Keep all 80 ancestry/study intervals / 2,225 candidate gene IDs, avoid retroactively promoting rs671 result to all phenotypes (AS/CKD/etc.).

## Reproducible code and location

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
export OPENBLAS_NUM_THREADS=1
python3 scripts/is/prepare_g0022_non_gtex_qtl_sources.py
python3 scripts/is/download_g0022_non_gtex_qtl_regions.py   # only if source files missing; six EBI targeted regions
python3 scripts/is/audit_g0022_non_gtex_qtl_variant_coverage.py
python3 scripts/is/audit_g0022_rs671_conditional_ld.py
python3 scripts/is/audit_g0022_jctf_japan_omics_variants.py
python3 scripts/is/audit_g0022_jctf_integrated_evidence.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Data/source:
`/srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1/non_gtex/`
`/srv/is-analysis/data/is/qtl/japan_omics/G0022_v1/`

Main result:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/G0022_EAS_JCTF_INTEGRATED_VARIANT_EVIDENCE.tsv`
`G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json`
`jctf_japan_omics_variants_v1/G0022_JCTF_4SNP_JAPANESE_EQTL_PQTL_EVIDENCE.tsv`
`alternate_qtl_sources_v1/G0022_NON_GTEX_CS_SOURCE_COVERAGE.tsv`
`G0022_AIS_RS671_CONDITIONAL_DIAGNOSTIC_SUMMARY.json`

## External sources

- rs671 Japanese ALDH2 multi-trait locus: https://pheweb.jp/variant/12-112241766-G-A
- Japan Omics Browser rs671 gene association table: https://japan-omics.jp/variant?input_value=rs671
- Japan Omics Browser eQTL / pQTL paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11525184/
- JOB browser methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC12057183/
- NBDC Japanese JCTF QTL unrestricted 2024 archive: https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000193
- MAEEA n2024 Asian eQTL: https://academic.oup.com/nsr/advance-article/doi/10.1093/nsr/nwag566/8779978
- MAEEA repository metadata (restricted): https://zenodo.org/records/21296030
- eQTL Catalogue official sources: https://www.ebi.ac.uk/eqtl/Data_access/
- ALDH2 rs671 catalytic and brain biology (not stroke causality proof): https://www.nature.com/articles/s41467-024-46899-0
