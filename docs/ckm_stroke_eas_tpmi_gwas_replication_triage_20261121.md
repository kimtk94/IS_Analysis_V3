# East Asian AIS GWAS replication source audit

Date: 2026-11-21. Technical source QC PASS. Scientific status HOLD. Draft PR #284 remains open and unmerged.

## Verified public data sources

### Taiwan Precision Medicine Initiative (TPMI) PheWeb, GRCh38

The official site https://pheweb.ibms.sinica.edu.tw/about states it provides GWAS summary statistics and PGS model resources from the TPMI Han Chinese program with GRCh38 variant coordinates. Official phenotype pages and case/control totals were independently retrieved in [GitHub Actions #44871010169](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871010169), source artifact [#14694240592](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871010169/artifacts/14694240592); 10/10 offline fail-closed tests passed.

| PheWeb phenotype | Observed source cases | Controls | Clinical interpretation |
| --- | ---: | ---: | --- |
| 433.21 Cerebral artery occlusion, with cerebral infarction | 9,249 | 304,660 | Highest-priority AIS-like infarction phenotype; not necessarily TOAST-imaging-confirmed |
| 433.2 Occlusion of cerebral arteries | 12,463 | 304,660 | Broader occlusion, can include no definite infarction |
| 433.3 Cerebral ischemia | 12,468 | 304,660 | Broad cerebral ischemia, not equivalent to acute ischemic stroke |
| China Kadoorie Biobank I63 | 14,302 | 67,954 | Ischemic-stroke code but participant overlap with earlier stroke consortia remains a problem |

Public pages: https://pheweb.ibms.sinica.edu.tw/pheno/433.21 , https://pheweb.ibms.sinica.edu.tw/pheno/433.2 , https://pheweb.ibms.sinica.edu.tw/pheno-filter/433.3 and https://pheweb.ckbiobank.org/pheno/i63 .

**Important**: no complete original outcome GWAS was downloaded; variant-level beta/SE, effect allele, imputation INFO, per-SNP effective N or full-file checksum were not checked. Merely finding a download link does not establish data access, complete variant schema, or license compliance. 2026 TPMI published original GWAS is not assumed identical to these public 2025 cohort PheWeb phenotypes.

### Distinguish original association statistics from polygenic score weights

PGS Catalog https://www.pgscatalog.org/score/PGS019250/ reports a score named TPMI_433.21_Lassosum2 and https://www.pgscatalog.org/score/PGS019234/ reports TPMI_433.2_SBayesR, both GRCh38. These are **polygenic score weights, not original SNP-level outcome GWAS beta/SE/P summary statistics**. The Catalog lists 30,664 association source cases for PGS019250, which differs from the live TPMI 433.21 phenotype page's 9,249 case count. Treat this as a source version/phenotype-cohort mapping issue to reconcile; do not conflate PGS evaluation, training, association and public PheWeb case-control cohorts.

### 2026 Taiwanese early-onset ischemic stroke study

Neurology 2026, DOI https://doi.org/10.1212/WNL.0000000000218126 , reports **5,546 cases and 143,017 controls** in a TPMI-based early-onset IS analysis and loci including rs541118668 near CYP4F3. This is **not established as the same summary-statistics dataset as PheWeb 433.21**. Its original full downloadable GWAS output, detailed cohort participant disjointness and SNP-wise N_eff remain unverified. The separate clinical TOAST-characterized patient group used for follow-up does not automatically redefine the public ICD-code-based discovery GWAS as TOAST-adjudicated.

### Other candidates and study origin risk

BBJ: OpenGWAS bbj-a-129 https://gwas.mrcieu.ac.uk/datasets/bbj-a-129/ lists ischemic stroke with 17,671 cases, 192,383 controls and 8,885,705 variants, but OpenGWAS VCF/API access can require a token and BioBank Japan participates in earlier major stroke consortia. CKB: official PheWeb I63 has 14,302 cases and 67,954 controls, but CKB is also represented in published EAS stroke meta-analysis; it cannot be called an independent GIGASTROKE replication without disjointness evidence.

GBMI 2023 https://doi.org/10.1016/j.xgen.2023.100345 has a 16-biobank meta-analysis of ischemic stroke (60,176 cases/1,310,725 controls) with BBJ and other overlapping data sources. The paper explicitly used a BBJ-excluded GBMI dataset when combining with MEGASTROKE due to overlap. A meta-analysis from GBMI or MEGASTROKE is not automatically independent from the original GIGASTROKE GWAS.

Published EUR GIGASTROKE original corrected supplement specifically includes FinnGen-named contributing cohorts (8,851 original AIS cases from two named source cohorts, or 14.25282% of reported EUR AIS cases); this is study-origin evidence and not a count of overlapping participants in FinnGen R13. See [FinnGen cohort audit](ckm_stroke_finngen_r9_r10_r13_overlap_audit_20261121.md).

## Source and scientific GO / HOLD matrix

| Gate | Status |
| --- | --- |
| Public TPMI 433.21/433.2/433.3 and CKB I63 phenotype metadata | PASS |
| PGS weights distinguished from original GWAS associations | PASS |
| Exact 2026 Taiwan early-onset GWAS equals 2025 TPMI PheWeb phenotype data | HOLD |
| Original whole-outcome summary statistics downloaded, source allele/beta/SE/P verified | NOT RUN |
| Original phenotype case definition independently validated against TOAST/MRI/CT | HOLD |
| Individual-level disjointness between TPMI and GIGASTROKE | HOLD |
| Original GIGASTROKE EUR effective N per variant | HOLD |
| New LDSC h2/intercept/rg, MR, colocalization | NOT RUN |

## Next practical scientific gates

1. Verify a usable official TPMI 433.21 association summary statistics download path, permissions/terms, original filename and header using a bounded source audit. The published site warns of download limits. Never repeatedly bulk-download or bypass restrictions.
2. Map public TPMI case ICD codes, controls and missingness for 433.21 to clinically defined AIS and ensure the 2026 Neurology GWAS is a different analysis if counts/build/source differ.
3. Confirm study disjointness from GIGASTROKE by comparing original contributing cohorts, calendar dates and overlap disclosures; do not claim that different biobank brand names mean disjoint participant IDs.
4. Prepare GRCh38-to-GRCh37 variant normalization, effect allele alignment, EAS-specific LD reference and candidate cis-pQTL MR/coloc gate separately. Do not use an arbitrary pooled case-control N for primary EUR LDSC.

Automation code: scripts/ckm_stroke/audit_eas_tpmi_ais_phenotype_sources.py; tests/test_ckm_eas_tpmi_ais_phenotype_sources.py; .github/workflows/ckm-eas-tpmi-ais-phenotype-source.yml. Draft PR #284 remains unmerged.


## Live official TPMI 433.21 original GWAS HEAD-only accessibility gate

A second source-backed [GitHub Actions #44871483498](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871483498) completed successfully with **14/14 negative and positive contract tests PASS**; [artifact #14694043597](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871483498/artifacts/14694043597). Only official original PheWeb phenotype HTML was downloaded (small pages), and **HTTP HEAD** rather than GET for the actual large GWAS result files.

| PheWeb 433.21 original association download | Advertised original source URL | HEAD | Reported compressed bytes | Result body downloaded? |
| --- | --- | --- | ---: | --- |
| Summary Statistics (PLINK2) | https://pheweb.ibms.sinica.edu.tw/download/gwas_plink2/433.21 | 200 | **58,302,393** | NO |
| Summary Statistics (SAIGE GWAS) | https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21 | 200 | **816,992,833** | NO |

Both HEAD requests remained on the official source hostname. HTTP 200 proves an advertised endpoint responded, **not** the file's internal format, full integrity, per-variant effect and sample-size columns, legal access authorization, or consistency with the 2026 Neurology early-onset IS discovery sample.

Important distinction: 2025 TPMI PheWeb **433.21 9,249 cases and 304,660 controls** versus 2026 published TPMI early-onset IS **5,546 and 143,017**. The PGS Catalog uses yet another association count for a score of the same named phecode; none is promoted as original paper-to-portal equivalence.

**Recommended next actual analysis gate:** obtain at most the smaller 58 MB PLINK2 source once, strictly respecting TPMI published download limits and terms, record response SHA256, verify the header and GRCh38 REF/ALT/BETA/SE/P/AF/missingness conventions, and compare the 433.21 ICD case definition. If the user opts to proceed with allowed downloads, carry out one reproducible CI or Drive-backed ingest. Do not run LDSC h2/rg or claim independent replication until overlap provenance and SNP-specific effective N are justified. The 817 MB SAIGE file is not required just to inspect access; avoid redundant downloads and do not bypass site throttling.


## Original TPMI 433.21 compressed PLINK2 association source downloaded ONCE — verified 2026-11-21

[**GitHub Actions #44871738962**](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871738962) **SUCCESS**, offline tests **9/9 PASS**. [Original header/full-compressed-byte SHA evidence artifact #14694277054](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44871738962/artifacts/14694277054). The official 58,302,393-byte original https://pheweb.ibms.sinica.edu.tw/download/gwas_plink2/433.21 was requested **one time only**, streamed to a temporary local file, full compressed bytes SHA-256 calculated, first 200 uncompressed association rows parsed, and source deleted before completing CI. Original source SHA-256 is **ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62**. Compressed type **gzip text**. The original source publisher does not provide an independent expected full-file checksum in the verified endpoint page; this hash pins the verified current retrieval but is not an external publisher-checksum comparison.

**Exact original header order:**

1. #CHROM
2. POS
3. ID
4. REF
5. ALT
6. PROVISIONAL_REF?
7. A1
8. OMITTED
9. A1_FREQ
10. TEST
11. OBS_CT
12. OR
13. LOG(OR)_SE
14. Z_STAT
15. P
16. ERRCODE
17. rsID

**Input interpretation:** The file provides an original **OR** and **LOG(OR)_SE** rather than a direct BETA column. A downstream, validated association preprocessing step could represent the effect as log(OR) for the original PLINK2 *A1* tested allele, but must first confirm TEST type, OR>0, source REF/ALT/A1 orientation, missing alleles, finite SE, P-vs-Wald consistency, possible variant duplication and association-model/error flags. None of these full-genome harmonization steps is marked PASS here. A1_FREQ is the source tested-allele frequency; it is not automatically MAF without computing min(f, 1-f).

**Per-SNP sample-size nuance:** OBS_CT **is present** in the source PLINK2 file, but is an observation count per variant/test, **not** a variant-specific case/control breakdown and **not by itself an approved effective N** for a case-control binary association. Across the first 200 variant records, all 200 P values parsed into [0,1] and there were 0 record-width/coordinate QC errors. The full decompressed genomewide association records were NOT processed or QC-approved in this preflight; first 200-row consistency does not prove whole-file content validity.

**Scientific status remains HOLD** for claiming independent EAS AIS replication: the 433.21 clinical ICD-code case definition, sample disjointness with GIGASTROKE, original 2026 Neurology stroke paper identity and whole-genome association QC remain pending. The original GIGASTROKE EUR SNP-wise meta-analysis effective N is still HOLD. No new LDSC h2/intercept/rg, MR, colocalization or causal inference is reported.

Workflow .github/workflows/ckm-tpmi-43321-original-header-once.yml, code scripts/ckm_stroke/inspect_tpmi_43321_original_sumstats_once.py, tests/test_ckm_tpmi_43321_original_once.py. This is a source ingest/QC milestone, **not research result replication**; the file was deliberately not redistributed or uploaded to repository artifacts to respect source use, download limits and space.


## 2026-11-21 UPDATED: full TPMI 433.21 original PLINK2 association QC and denominator discordance

The earlier 200-row original-file source check above has been **superseded for full row coverage** by [source-backed GitHub Actions #44872892525](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44872892525) **SUCCESS** with 11/11 offline source mutation and error tests PASS. Original-source aggregate artifact [#14694881855](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44872892525/artifacts/14694881855) includes no per-SNP effect values. It read the original compressed source from the same official TPMI PheWeb PLINK2 GWAS URL **once**, verified its complete **58,302,393 bytes** and pinned prior SHA256 **ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62**, decompressed and validated all original 1,739,048 association TSV records across 22 chromosomes, and deleted the original local gzip after the audit. The SHA is our repeat-verification anchor, not a separately published TPMI expected checksum. No SAIGE 817 MB file or variant effects were uploaded.

| Original source metric | Reproduced value |
| --- | ---: |
| All original chr1–22 variant/test rows | **1,739,048** |
| PASS numeric ADD + unflagged test (ERRCODE ".") + REF/ALT/A1 + AF + OR>0 + LOG(OR)_SE>0 + Z + valid P | **1,739,048** |
| Duplicate variant/test 128-bit digest | 0 (cryptographic digest collision caveat; not exact string identity proof) |
| Coordinate order violations | 0 |
| SNP P <= 5e-8 | **1** (aggregate count, not locus identity or replication evidence) |
| Original A1_FREQ MAF >= 0.05 | **1,235,405** |
| 0.01 <= MAF < 0.05 | **502,760** |
| 0.001 <= MAF < 0.01 | **883** |
| P=0 numerical underflow | 0 |
| P-versus-Z | All 1,739,048 source rows numerically comparable; none differed by >=0.1 in -log10(P) |
| Observed PLINK2 variant-wise OBS_CT | **204,099 exactly** in all 1,739,048 rows |

CRITICAL: Original PheWeb phenotype page reports 9,249 cases and 304,660 controls (**313,909 sum**), whereas the original **PLINK2 regression output OBS_CT is constant 204,099**. The gap is **109,810 samples, 34.98148% of published phenotype sample sum**. This is **not** proof of cohort data error, and it is not resolved by treating the page counts as the PLINK2 regression counts. The actual number of cases among those 204,099 regression observations is unknown. Never compute an "approved" binary N_eff from either 313,909 or 204,099 without source-specific analysis sample case/control metadata. A valid source N_eff case-control analysis requires Ncase and Ncontrol for the actual analysis, not just total OBS_CT.

Potential explanation (not proven): TPMI's 2025 Nature cohort methods report 205,779 independent versus 62,304 related participants, and PLINK2 Firth association on an independent subset: https://www.nature.com/articles/s41586-025-09680-x . The 204,099 original PLINK2 regression observation count is close to but 1,680 below the paper's 205,779 independent individuals. One cannot assume equality or infer the actual case/control breakdown, and the PheWeb page phenotype N may refer to different cohort inclusion or SAIGE results. This is a reproducible **sample denominator provenance HOLD**.

Independent fail-closed re-evaluation of the entire original aggregate JSON (without re-downloading the 58MB source) was performed by [GitHub Actions #44873247155](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873247155) **SUCCESS**, 9/9 contract tests PASS, source-derived [artifact #14695131212](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873247155/artifacts/14695131212). It reproduces chromosome/MAF and effect QC totals, single genome-wide significant P, and explicitly preserves the 34.98148% sample mismatch and scientific HOLD.

**Methodological next gate:** obtain original PLINK2 actual-analysis Ncase/Ncontrol for 433.21 from official authors/metadata (or reproducible analysis cohort source), clinically validate phecode to AIS, independently verify GIGASTROKE participant disjointness, and then implement ancestry-specific EAS LD plus GRCh38/37 and original A1/logOR harmonization. A single aggregate genomewide-significant hit is not sufficient to claim population-level independent replication or validate a causal MR instrument. No new LDSC h2, intercept, rg, MR or coloc executed. The previously verified GIGASTROKE EUR SNP-specific N_eff remains unresolved independently of this TPMI audit. Draft PR #284 remains unmerged.

Reproducible code and workflows:
- scripts/ckm_stroke/audit_tpmi_43321_full_gwas_qc.py and tests/test_ckm_tpmi_43321_full_gwas_qc.py; .github/workflows/ckm-tpmi-43321-full-gwas-qc.yml
- scripts/ckm_stroke/verify_tpmi_43321_source_obsct_discrepancy.py and tests/test_ckm_tpmi_43321_obsct_discrepancy.py; .github/workflows/ckm-tpmi-43321-obsct-source-gate.yml


## 2026-11-21 SOURCE-VERIFIED TPMI sentinel rs187759247 and four-original GIGASTROKE EUR locus check

**Completed successful source-backed GitHub Actions:**

- [#44873979328 — one original TPMI source lead + bidirectional official Ensembl GRCh37 mapping](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873979328) SUCCESS, 8+8 offline tests PASS; public original TPMI lead [artifact #14695301070](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873979328/artifacts/14695301070), Ensembl exact map [artifact #14695296126](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873979328/artifacts/14695296126).
- [#44874466311 — four original publisher-MD5 GIGASTROKE EUR association files and independent four-trait review](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44874466311) SUCCESS, 10/10 offline test cases PASS, all four publisher original full compressed MD5 + original outcome whole-row count PASS. Combined comparison [artifact #14695641965](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44874466311/artifacts/14695641965). Shards read 1,014 MB of original GWAS gzip source in total, but uploaded only compact public variant presence/absence and original checksum metadata. The original files were not redistributed.

**TPMI 433.21 original PLINK2 GWAS lead variant (one genomewide SNP with P <= 5e-8, among 1,739,048 original records):**

| Field | Original association |
| --- | --- |
| rsID | **rs187759247** |
| GRCh38 | **chr19:15,462,568 T>C** |
| Original ID | chr19_15462568_T_C |
| A1 tested effect allele | **C** (original REF T, ALT C) |
| A1_FREQ | **0.0212101** |
| OR per C | **1.68467** |
| log(OR) | **0.5215696989** |
| log(OR) SE | **0.0662893** |
| Z | **7.86808** |
| P | **3.60115×10^-15** |
| PLINK2 OBS_CT | **204,099** (binary effective N not derived) |
| GRCh37 mapped coordinate | **chr19:15,573,379 T>C** |
| GRCh38→GRCh37→GRCh38 | **PASS unique bidirectional coordinate** |
| Official GRCh38 and GRCh37 Ensembl reference T | **PASS both genome builds** |
| A1 C preserved | **PASS same forward orientation** |

Coordinate conversion used official Ensembl REST mapping https://rest.ensembl.org/documentation/info/assembly_map with a reciprocal https://grch37.rest.ensembl.org map and both reference sequence bases. This establishes only the genotype reference locus and effect-allele orientation; it does not establish genomic or clinical relevance of the SNP.

**Original published GIGASTROKE EUR exact-position check — ALL original compressed MD5 validated:**

| GIGASTROKE phenotype | GWAS accession | Original source association rows scanned | chr19:15,573,379 rows with original REF/ALT | Clinical inference |
| --- | --- | ---: | ---: | --- |
| AIS | GCST90104540 | 7,482,032 | **0** | NOT OBSERVED IN PUBLISHED EUR SOURCE |
| LAA | GCST90104542 | 5,784,788 | **0** | NOT OBSERVED IN PUBLISHED EUR SOURCE |
| CE | GCST90104541 | 6,659,793 | **0** | NOT OBSERVED IN PUBLISHED EUR SOURCE |
| SVD | GCST90104543 | 5,974,028 | **0** | NOT OBSERVED IN PUBLISHED EUR SOURCE |

**Interpretation:** TPMI 433.21 has a genuine, single strong public association in the validated Han Chinese GRCh38 PLINK2 summary statistics; it has no *direct same-coordinate* comparison in the four public GIGASTROKE EUR original GWAS sources. “Not present in a particular published summary file” **does NOT mean the allele is absent from all Europeans, does NOT establish a null effect, and must not be coded as beta 0, P=1 or nonreplication**. Original EUR GWAS variant QC/MAF/imputation/participating-cohort coverage could explain missingness; these reasons are NOT verified. rs187759247 was not shown to be causal or the relevant biological driver, and the one lead SNP cannot support MR/colocalization by itself.

**Other unresolved gates:**

1. TPMI PheWeb 433.21 has 9,249 cases and 304,660 controls (total 313,909); every original PLINK2 variant has OBS_CT 204,099. Difference 109,810 (34.98148%) is still **not reconciled by actual analysis-case/control source counts**. Neither denominator can be silently used to assign effective case/control N.
2. TPMI phenotype 433.21 is clinically coded cerebral artery occlusion with infarction, not strictly adjudicated GIGASTROKE TOAST AIS. Exact clinical case and exclusion definitions plus the 2026 Neurology early-onset analysis lineage remain unverified.
3. TPMI individual-level sample disjointness from all GIGASTROKE contributors has NOT been demonstrated. Comparison across EAS vs EUR additionally requires suitable ancestry-specific LD and a pre-registered harmonization/replication design.
4. Original GIGASTROKE EUR meta-analysis per-SNP N_eff, actual LDSC post-munge regression variants and published study correction provenance are still HOLD. No LDSC h2/intercept/rg, MR or colocalization was calculated in this step.

**Best next research gate:** assess whether a validated **LD proxy** exists for rs187759247 in EUR and/or EAS and whether the locus appears in another genuinely independent clinically defined EAS AIS GWAS. Do not choose an LD proxy only by physical distance, and do not impute a direct EUR original source association without appropriate reference LD. Separately, resolve 433.21 original analysis case/control counts. The scientific replication decision remains **HOLD**, although original source and coordinate technical QC are now PASS.

Reproducible scripts: scripts/ckm_stroke/identify_tpmi_43321_lead_variant.py, scripts/ckm_stroke/map_tpmi_43321_lead_ensembl.py, scripts/ckm_stroke/probe_gigastroke_tpmi_43321_lead.py, scripts/ckm_stroke/review_tpmi_gigastroke_four_targeted.py; CI workflows .github/workflows/ckm-tpmi-43321-lead-grch37.yml and .github/workflows/ckm-tpmi-gigastroke-rs187759247-targeted.yml. Existing Draft PR #284 must remain unmerged without review.
