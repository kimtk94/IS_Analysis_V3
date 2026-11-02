# Ocular P5j — exact 2026 IOVS human TM gene-level replication crosswalk (VERIFIED 2026-10-31)

## Execution and verified provenance
- **[P5j GitHub Actions #41158150851](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158150851) PASS**; aggregate-only [artifact #13025593642](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158150851/artifacts/13025593642). This independently fetched official public Europe PMC supplementary package and read the verified `iovs-67-2-13_s003.xlsx`; matched gene-level statistics to our prior P5b and P5d public aggregate-only artifacts. The preceding initial run #41157994029 also PASS but classified missing external adjusted p-values as non-significance; **use corrected P5j #41158150851 only**.
- Publication: Mehrotra et al., *Investigative Ophthalmology & Visual Science* 2026, PMID **41636427**, PMCID **PMC12882097**, Europe PMC authoritative supplementary package https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12882097/supplementaryFiles.
- Workbook sheets: `SuppTable3_TM2Days_HighConfDEGs` (857 data genes; 858 total rows including header); `SuppTable4_TM2Days_All_DEGs` (13,313 data genes; 13,314 total rows including header), `SuppTable6_SC2Days_HighConfDEGs` (2,086 data genes) and `SuppTable7_SC2Days_All_DEGs` (18,423 data genes). SC/SCE gene universe **must not be treated as the independent TM validation set**.
- Primary: GSE313736 57 paired donors, BeamA high-confidence `score>=.75`, 114 paired donor-condition pseudobulks, robust edgeR NB QL, 16,050 genes tested. Official external TM bulk used 100 nM DEX at 48h; GSE313736 nominal 500 nM DEX at 48h. Different culture conditions, cell-state mixtures and analysis models warrant limited effect-size comparability.
- Source code: `scripts/ocular/ocular_p5j_official_TM_overlap.py`; reproducible P5g shortlist from `scripts/ocular/ocular_p5g_candidate_screen.py`; bundled workflow `.github/workflows/ocular-p5j-published-overlap.yml`. Runs use **only aggregate gene-level published data**, not donor IDs or controlled genotypes.

## Full published human TM DEG set overlap
| Metric | Verified count / fraction |
|---|---:|
| Our primary genes tested | **16,050** |
| Official external TM all-gene table | **13,313** |
| Official published high-confidence TM DEGs | **857** |
| Ensembl gene IDs shared between the two all-gene sets | **12,251** |
| Primary BH-FDR<0.05 within shared IDs | **10,619** |
| Official TM BH-FDR<0.05 within shared IDs | **824** |
| Shared genes BH-FDR<0.05 in **both** analyses | **755** |
| Shared genes BH-FDR<0.05 **and** signed log2FC direction same | **659** |
| Signed replication fraction among 824 externally published high-confidence TM genes present in our primary tested set | **659/824 = 80.0%** |
| Signed concordance among 755 FDR-significant in both studies | **659/755 = 87.3%** |
| Opposite direction despite FDR significance in both | **96/755 = 12.7%** |
| Across all 12,251 shared genes: signed agreement | **6,951/12,251 = 56.74%** |
| All shared genes' log2FC Pearson r | **0.461345** |
| Primary high-effect abs(log2FC)>=1 among shared genes (both have available external FDR) | **772**, signed agreement **84.715%** |
| External TM genes with **missing adjusted p-values** in its all-gene table | **259** |

**Interpretation:** Strong high-confidence TM DEX molecular responses have substantial cross-study signed concordance, despite different dose and cell-state technologies. This is **published-result-level external concordance**, **not** an independent raw FASTQ reanalysis, a causal separation of DEX vs libraries, a matched donor proof, or a clinical pressure-response model. The global r=0.46 and only 56.7% same-sign among all shared genes prevent blanket transcriptome-wide replication claims. The 96 externally significant, direction-discordant genes warrant a separate biology/technical audit.

## P5g statistical robust top-86 candidates against official 857 TM DEGs
| External paper status for gene | Genes among our 86 |
|---|---:|
| **Already published in the 857 high-confidence TM DEG set** | **38** |
| Present in official TM all-gene table, external BH-FDR>=0.05 | **32** |
| Present in official TM all-gene table, **adjusted p unavailable** | **2** |
| **Not present in official TM all-gene tested universe** | **14** |
| **Total** | **86** |

- All **38/38** overlapping official high-confidence TM DEGs have **the same signed direction** in our high-stability shortlist. The 32 externally non-significant candidates have **28 same sign / 4 opposite sign**. The two with unavailable FDR have same signed effect but **cannot be classified as non-significant**.
- Prior manual literature audit had flagged 13 of 86 as known effects; 10 of those are part of the 38 official-high-confidence hits, and the other 3 (SAA1, SAA2, ZBTB16) had existing earlier-study support but are not in this official all-TM table. Hence at least **41/86** now have some specific prior evidence of TM/steroid expression response. The other **45** are **not proven novel**. Manual literature coverage and gene measurement scope are incomplete.
- Examples from 86: HIF3A, SLC16A12, NEBL, OTULINL absent from the external TM tested set (**not novel by absence**); CORIN present with a weaker same-direction effect but external FDR=0.1601 (not a strict statistical replication); CITED1 and MRO have external signed effects but **missing adjusted p-values**, so hold any significance conclusion. LINC01088, AFAP1L1 and CHI3L2 were already among official high-confidence TM DEX DEGs.
- These 86 are a **post-hoc exploratory prioritization of original P5 top100**, not an independently pre-registered genome-wide novel candidate list. Gene-level screen thresholds and manual literature flags are explicitly retained.

## GO/HOLD/STOP
- **GO** original raw-count pseudobulk and edgeR, 9/9 within-cohort leave-one-pool-out signed stability, 9 experimental-library-unit sensitivity and **official published independent study directional gene-level evidence** for 659 molecules. Note batch treatment confounding still applies.
- **HOLD** dose/time harmonization, original study paper's complete high-effect 63-up/17-down exact gene-level comparison, mechanistic novelty (including 45 insufficiently reviewed candidates), prospective external independent raw-data replication and causal independence from common library preparation.
- **STOP** clinical steroid-induced intraocular pressure prediction from this in-vitro GSE313736 experiment; no observed donor-level clinical delta-IOP phenotype linked.

## Next reproducible study gate
1. Compile the 96 genes **significant in both studies but opposite signed** into a conservative gene-level *derived* discrepancy audit using the P5j summary/official workbook; explicitly evaluate isoform annotations, treatment dose, cell state, reference assembly and per-pool effects.
2. Rank unreviewed candidates by lack of previous expression evidence **without** converting 'not tested' or 'non-significant in another study' into 'new'.
3. Investigate original GSE313736 63-up/17-down exact full supplement and sample attrition reasons (74→70→59), and perform orthogonal cell-state annotation sensitivity.
4. Clinical MR/coloc/IOP work must use independently accessed patient-level outcome/QTL evidence; these expression data alone cannot establish clinical mechanism.
