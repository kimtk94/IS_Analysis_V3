# IS — 30 East Asian GWAS regions and public nominal QTL expansion (2026-10-11)

**Branch:** `research/is-broad-discovery-20261009`
**Research status:** COMPLETED *input coverage* and source provenance; **coloc/fine-mapped shared causal variants NOT newly established**.
**Preserved canonical scope:** 2,225 distinct positional gene IDs, 30 provisional East-Asian cross-dataset genomic groups, 646 old BBJ–GTEx ABF assays, 8,956 old source-provenance SNPs.

## 1. Methods and provenance

The 30 original `IS_XDATA_G0001–G0030` groups were loaded from the frozen `CROSS_DATASET_INTERVAL_GROUPS.tsv`. These are **distance-clustered provisional GRCh37 groups**, **not independently LD-confirmed or independent causal loci**. Individual groups derive from stroke subtype sources (AS, AIS, CES, LAS, SVS, BBJ ischemic stroke) and some are suggestive rather than genome-wide significant.

The full local **processed/canonical** GWAS summary files were streamed once each, without new downloads:

| Source | Records scanned | Interpretation |
|---|---:|---|
| BBJ Japanese ischemic stroke | 13,435,531 | Original **normalized canonical** summary |
| GIGASTROKE EAS any stroke (AS) | 6,789,909 | Different phenotype from AIS |
| GIGASTROKE EAS ischemic stroke (AIS) | 6,725,216 | East Asian meta-analysis |
| GIGASTROKE EAS cardioembolic stroke (CES) | 5,169,600 | Subtype |
| GIGASTROKE EAS large artery stroke (LAS) | 5,993,005 | Subtype |
| GIGASTROKE EAS small vessel stroke (SVS) | 6,173,637 | Subtype |
| **Six files, total** | **44,286,898** | **30 × 6 = 180 region–source records** |

For each group, `[original GRCh37 interval start−500kb, end+500kb]` was used as a **contextual lookup window**, not a fine-mapped locus boundary. Within each window we retained number of records, minimum P, number of P<5e−8 records, and source-group original lead-variant presence. The 44 locally canonical BBJ records failing this strict per-row SNP/allele/metadata check were counted as **unresolved**, rather than rescued or silently used.

All outputs write to new immutable audit directories; SHA256 source manifests accompany output. **No archived 646-coloc posteriors, SNP sets, research canonical data, or runtime trading/production systems were changed.**

## 2. EAS AIS genetic signal tiers: 30 groups ≠ 30 validated causal loci

The source-derived region workflow queue has four tiers:

| Tier | Region count | Rule | Inference boundary |
|---|---:|---|---|
| **T1** | **4** | Both BBJ ischemic stroke **and** GIGASTROKE EAS AIS have a P<5e−8 variant within the same predeclared contextual interval | Does **not** imply same SNP, independent replication, or shared causal variant |
| T2 | 1 | GIGASTROKE EAS CES subtype P<5e−8 but not EAS AIS | CES-specific hypothesis, not broad AIS |
| T3 | 2 | EAS AS P<5e−8, but neither EAS AIS nor BBJ genome-wide significance in that window | Any-stroke association is not automatically AIS |
| T4 | 23 | No genome-wide signal meeting the above source-derived criteria | Retain as **EXPLORATORY**; not negative results |

The four T1 groups, strongest outcome-supported priority for an additional **QTL source-availability** audit:

| Group | Best BBJ P in ±500kb | Best EAS AIS P in ±500kb | EAS positional gene IDs | With archived GTEx ABF | Without archived ABF |
|---|---:|---:|---:|---:|---:|
| `IS_XDATA_G0022` | 1.9919e-18 | 8.426e-18 | 112 | 16 | **96** |
| `IS_XDATA_G0015` | 8.7620e-10 | 5.584e-11 | 71 | 17 | **54** |
| `IS_XDATA_G0023` | 3.8085e-10 | 1.007e-08 | 14 | 6 | **8** |
| `IS_XDATA_G0007` | 1.2119e-08 | 2.583e-08 | 10 | 4 | **6** |
| **T1 total** | — | — | **207** | **43** | **164** |

In contrast, `G0008` has EAS CES P≈6.963e−11 but BBJ and EAS AIS contextual minima around 0.0013: this is **not** an AIS-replicated signal. `G0004` and `G0021` are EAS AS-supported without a genome-wide EAS AIS window hit. `G0025`, despite its 92 candidate genes, is only suggestive for BBJ (P≈9.87e−8) with EAS AIS contextual P≈3.20e−4. **Gene counts are not causal prioritization scores.**

### Original lead SNP source-allele check

Exactly **28 of 30** cross-dataset group lead sites were reported in both BBJ and EAS AIS source GWAS summaries. Source allele pairs required **REF/ALT order reorientation in 14 of those 28** comparisons against the group-lead representation (11 BBJ-side, three AIS-side). After reorienting **only the recorded allele pair and beta sign, not complement-rescuing ambiguous strands**, the aligned effect direction agreed in **28/28**. The other two group leads (`G0009` and `G0016`) were not reported in the EAS AIS canonical file.

**Critical caveat:** BBJ and GIGASTROKE source cohort membership may overlap, so **28/28 sign concordance is neither independent replication nor evidence of a causal gene**; palindromic sites require original genomic reference validation before downstream allele harmonization.

## 3. New publicly obtainable nominal cis-eQTL source extraction

Only public, nominal (including nonsignificant association records) regional extracts from eQTL Catalogue were downloaded by indexed HTTPS `tabix`; **no controlled-access, paid CKB/TPMI or private genotypes**.

Three whitelisted datasets, GRCh38 coordinates:

- **QTD000131**, GTEx v8 Artery Aorta, catalogue N=387
- **QTD000171**, GTEx v8 Brain Cortex, catalogue N=205
- **QTD000434**, ROSMAP brain DLPFC, catalogue N=560

Official sources:
- `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000131/QTD000131.all.tsv.gz`
- `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000171/QTD000171.all.tsv.gz`
- `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz`

GRCh37 original group boundaries ±500kb were projected with the existing UCSC hg19→hg38 chain; both boundaries and internal group anchors required unique positive-strand monotonic mappings. **This is query-window alignment only**, not independent per-SNP GRCh38 reference/allele validation. The returned 19-column QTL variant tuples were schema-checked. Each raw file has a matching recorded SHA256 hash and was written to a new source directory. The three datasets' molecular-QTL ancestry and independence from BBJ/GIGASTROKE **remain unverified**.

| Provisional group | GWAS source category | Positional genes | Raw nominal QTL gene–SNP rows across 3 datasets | Gene IDs with one or more nominal-source rows (union) | Newly sourced nominal genes **without old ABF** |
|---|---|---:|---:|---:|---:|
| G0015 | T1 BBJ + EAS AIS | 71 | 645,751 | 46 | **30** |
| G0023 | T1 BBJ + EAS AIS | 14 | 352,866 | 9 | **3** |
| G0007 | T1 BBJ + EAS AIS | 10 | 58,691 | 4 | **1** |
| G0025 | T4 suggestive | 92 | 978,647 | 61 | **61** |
| G0018 | T4 suggestive | 53 | 641,157 | 29 | **29** |
| **All 5 pilot groups** | | **240** | **2,677,112** | **149** | **124** |

All **15 of 15** requested region–dataset public extracts succeeded; there were **no source download or schema errors**. Newly sourced 124 genes **do not have any new coloc, MR, functional verification or causal assignment** in this workflow. In all five regions, raw QTL availability does not prove that the SNPs overlap the relevant GWAS after allele + assembly harmonization.

### Consequence for the full 869 EAS positional genes

- **43:** previously had at least one of the 646 archived BBJ–GTEx ABF assays.
- **124:** were never in the 646 ABF assays, but now have one or more original nominal QTL rows in the targeted public datasets.
- **702:** neither archived ABF nor nominal-source rows from **these exact five-region/three-dataset pilots**. Some have other existing targeted QTL extracts; this label is **NOT biologically negative**, nor a claim that public QTL does not exist.
- The **original 826 EAS genes without archived ABF remain 826**; the 124 are newly source-ready, not newly coloc-tested.
- **1,356 other positional genes** remain in the non-EAS-region research queue; 2,225 unique IDs are unchanged.

For T1 four regions, **164 EAS genes** still lacked archived ABF. **34** of those acquired new targeted nominal QTL source rows in G0015, G0023 and G0007, while other T1 sources (`G0022`) already have earlier **partial local** QTL datasets whose completeness has **not** been asserted by this pilot. Do **not** subtract 34 from the old ABF-untested count; **no new ABF model was run**.

## 4. Figure and data products

- [30 provisional EAS regions: genome-wide outcome support vs archived QTL gaps](../figures/is/IS_EAS30_GWAS_QTL_COVERAGE_CONTEXT.svg)
- [869 EAS positional gene IDs: 43 archived, 124 nominal-source new, 702 not examined in this pilot](../figures/is/IS_EAS30_PUBLIC_QTL_GENE_READINESS.svg)
- `docs/data/is/IS_EAS30_OUTCOME_QTL_EXPANSION_PRIORITY.tsv`
- `docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv`
- `docs/data/is/IS_FIVE_EAS_REGIONS_NOMINAL_GENE_UNION.tsv`

Reproducible scripts (one-shot, non-overwriting output directories):
- `scripts/is/audit_is_30_eas_regions_six_gwas_sources.py`
- `scripts/is/audit_is_eas30_crossgwas_source_leads.py`
- `scripts/is/build_is_eas30_gwas_qtl_expansion_queue.py`
- `scripts/is/pilot_is_eas_two_region_public_nominal_qtl.py`
- `scripts/is/pilot_is_eas_t1_three_region_public_nominal_qtl.py`
- `scripts/is/merge_is_eas30_public_nominal_qtl_15_sources.py`
- `scripts/is/render_is_eas30_gwas_qtl_coverage.R`
- `scripts/is/render_is_eas30_2225_public_qtl_gene_readiness.R`

Server scientific output roots:
- `/srv/is-analysis/results/is/audits/is_eas30_six_outcome_source_coverage_20261011_v1`
- `/srv/is-analysis/results/is/audits/is_eas30_source_lead_alignment_20261011_v2`
- `/srv/is-analysis/results/is/audits/is_eas30_gwas_qtl_priority_20261011_v3`
- `/srv/is-analysis/results/is/audits/is_eas30_two_region_nominal_qtl_20261011_v1`
- `/srv/is-analysis/results/is/audits/is_eas30_T1_three_region_nominal_qtl_20261011_v1`
- `/srv/is-analysis/results/is/audits/is_eas30_15_public_qtl_2225_readiness_20261011_v1`
- QTL raw extracts: `/srv/is-analysis/data/is/qtl/eqtl_catalogue/eas30_public_two_region_pilot_20261011_v1` and `.../eas30_T1_three_region_pilot_20261011_v1`

Full public QTL row-level output files remain on the server, **not in Git**; only read-only small derived summaries, plots, documentation, scripts and tests are committed.

## 5. Next science gate before any new coloc result

1. **Per-SNP GRCh38 reference FASTA and GRCh38↔GRCh37 allele harmonization**, checking biallelic SNV/indels, swapped/ambiguous palindromic variants, duplicates, and actual source effect allele and `AC/AN` vs `MAF`. Do not automatically flip based on high population AF disagreement.
2. Match each nominal QTL gene–SNP test to BBJ and GIGASTROKE EAS AIS **at identical, reference-aware variants**, including nonsignificant tests. Count and report matched SNPs, GWAS information/MAF, QTL `SE`, N and allele-frequency discrepancies. Do not run coloc for genes without enough overlapping high-quality SNPs.
3. Distinguish full tested cis-QTL windows from published significant-only, window-truncated, or expression-undetected records; preserve QTL `sample_size` and actual sample ancestry/source, and verify GWAS/QTL cohort overlap and LD suitability.
4. If all gates pass, use predeclared baseline `p12=1e-5` and five-prior sensitivity in a **separate new analytical branch**, with matched-ancestry molecular-QTL LD required for multi-signal `coloc.susie`. Do not retrofit 646 archived input SNP sets.
5. Revalidate G0022's earlier partial public molecular-QTL datasets and the ALDH2 conditional effect gates independently. Keep ALDH2 mediation and causal conclusions BLOCKED.
6. Maintain all 2,225 positional candidates in manuscript/dashboard with distinct `ARCHIVED_ABF`, `NEW_NOMINAL_SOURCE_ONLY`, `NOT_IN_PILOT`, `NOT_ASSESSABLE` statuses; no candidate is promoted solely by nearest-gene positioning or healthy-control cell type expression.

**Computational outcome:** original fixed 646 ABF and prior sensitivity results unchanged. Results here constitute **new public data availability and source variant harmonization audit only**, not new causal inference.
