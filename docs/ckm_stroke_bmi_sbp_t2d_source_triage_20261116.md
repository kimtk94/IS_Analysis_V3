# CKM x ischemic stroke — EUR BMI/SBP/T2D candidate source triage

Evidence date: 2026-11-16. This document complements `docs/ckm_stroke_feasibility_20261116.md`.
**All entries remain CANDIDATE_ONLY**; do not copy them into approved `public_sources.tsv` or unlock LDSC without exact original source QC.

## EUR candidate shortlist

| Trait | Candidate | Original/curated provenance | Decision |
| --- | --- | --- | --- |
| BMI | GIANT + UK Biobank (Yengo et al. 2018), corrected source | GIANT consortium [official data page](https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files) explicitly warns pre-2018-06-25 meta-analysis BMI files are not full-N and were corrected. ~700k EUR subjects. | SOURCE_CANDIDATE / HOLD original file version, full genome SNPs, actual header/GRCh build |
| SBP | UK Biobank + ICBP (Evangelou et al. 2018), OpenGWAS `ieu-b-38` | [IEU / OpenGWAS dataset](https://gwas.mrcieu.ac.uk/datasets/ieu-b-38/) reports EUR N=757,601, 7,088,083 variants, GRCh37. Cross-reference [knowledge portal](https://kp4cd.org/node/679). OpenGWAS web redirects to a dynamic token-limited site. | SOURCE_CANDIDATE / HOLD source entitlement, original header, cohort overlap |
| T2D | T2DGGI Suzuki et al. 2024 EUR | [DIAGRAM/T2DGGI consortium downloads](https://www.diagram-consortium.org/downloads.html) lists the EUR ancestry GWAS (~365 MB) and publisher MD5 `46f13bf073c5560d9b068b623229ff01`. | SOURCE_CANDIDATE / HOLD original binary, precise build, N, effect scale, terms |
| T2D sensitivity | DIAMANTE Mahajan et al. 2022 EUR | Same [consortium downloads](https://www.diagram-consortium.org/downloads.html), EUR ancestry summary file (~191 MB) publisher MD5 `157743a89874b91ed078270c97e22f5b`; [Knowledge Portal](https://kp4cd.org/index.php/node/872) reports 80,154 cases and 853,816 controls. | CONDITIONAL sensitivity / HOLD |

No N/effect scale from the general multi-ancestry T2D report should be substituted for the EUR-specific dataset. Published MD5s are **not yet independently matched to downloaded bytes**. Data use consent/download checkbox must be respected; no bypass.

## Prerequisites before statistical release

1. Pin each chosen original **specific file**: immutable direct URI or exact accession; published file metadata; download permission; SHA256/MD5 of actual bytes.
2. Confirm genome build (GRCh37 for outcome and LD where possible), alleles, effect direction, unit/scale, P, SE, per-variant/effective N, frequency and treatment of X/multiallelic variants.
3. Confirm sample overlap between GIANT/UKB, Evangelou/UKB+ICBP, T2D meta-analyses and GIGASTROKE; use LD Score Regression intercept/reference overlap strategy, not naive independent-cohort assumptions.
4. Separate binary T2D liability/log-odds scale and effective sample size from continuous BMI/SBP, including ascertainment, overlap and response transformations.
5. Compare 2024 T2DGGI vs 2022 DIAMANTE: ancestry-specific EUR releases only, no multi-ancestry pooling for EUR LD.
6. Re-run ancestry-matched full autosome QC, LDSC `h2` checks and intercept before `rg`, then formal subtype-contrast tests. Verify source-provenance ledger before any colocalization.
7. GLGC 22/21/20 partial QC is 9/66 trait-chromosome observations, not 22/22 or LDSC PASS.

`config/ckm_stroke/p1_gwas_candidates.tsv` is a **candidate-only manifest**.
`scripts/ckm_stroke/audit_p1_gwas_candidates.py` reads only bounded consortium landing-page text, verifies exact markers and redirects, and emits a JSON **HOLD** artifact. Page-marker PASS never implies an approved source. No full GWAS download.

## Operational choice

Until the saturated GitHub Actions queue clears, run metadata probes locally inside the Orca research terminal; use GitHub Actions for reproducible scientific jobs after inputs pass verification. Do not dispatch duplicate large-source jobs. Maintain Draft PR #284 unmerged.

## Orca live metadata probe, 2026-11-16

Executed in the **Orca-managed IS_Analysis_V3 worktree terminal**, not inferred from a web search:

| Trait/priority | Landing page probe | Interpretation |
| --- | --- | --- |
| BMI primary | LANDING_MARKER_FOUND | GIANT source page exists; corrected full-coverage file, genome build and original checksum remain unverified |
| SBP primary | HOLD_MARKER_NOT_FOUND | OpenGWAS redirects to dynamically rendered app; absence of server HTML marker is not dataset absence |
| T2D primary (2024) | LANDING_MARKER_FOUND | DIAGRAM/T2DGGI official download list exists; actual EUR binary, cohort details and published checksum remain unverified |
| T2D secondary (2022) | LANDING_MARKER_FOUND | Same consortium lists independent older EUR candidate; scientific sensitivity only |

Regression tests: 4/4 PASS in Orca. This probe does **not** change `public_sources.tsv`, does not use OpenGWAS authentication, and downloads no binary GWAS.

## 2026-11-16 follow-up — pinned original file and bounded schema evidence

This follow-up distinguishes **exactly resolved original source URI**, **actual compressed header verification**, **third-party curated metadata**, and **consent/authentication-gated access**. `config/ckm_stroke/p1_exact_source_endpoints.tsv` is an unapproved provenance ledger, not an approved input manifest.

| Trait | Source link result | Header / build status | Decisions |
| --- | --- | --- | --- |
| BMI | Resolved exact **corrected** GIANT 2018 binary: [Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz](https://giant-consortium.web.broadinstitute.org/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz) | Live Orca `HTTP 206`, server reports **46,754,937 total bytes**, read **16,384 compressed bytes** only. Actual decompressed header has all 10 columns: `CHR POS SNP Tested_Allele Other_Allele Freq_Tested_Allele_in_HRS BETA SE P N`. Official [Yengo 2018 README](https://giant-consortium.web.broadinstitute.org/images/0/01/README_summary_statistics_Yengo_et_al_2018.txt) confirms **hg19**, tested allele is effect allele, HRS frequency is from **8,552 unrelated HRS participants**, not the entire GWAS meta-cohort. | **PASS_HEADER_ONLY**, NOT FULL SOURCE / NOT LDSC APPROVED |
| SBP | [NHLBI GRASP Full Results](https://www.nhlbi.nih.gov/grasp/full-results) lists `UKB-ICBPmeta750k_SBPsummaryResults.txt.gz`; a [peer-reviewed article](https://doi.org/10.1186/s12916-022-02483-1) cites the old [GRASP file URL](https://grasp.nhlbi.nih.gov/downloads/ResultsAug2019/2018/Evangelou/UKB-ICBPmeta750k_SBPsummaryResults.txt.gz). The current URL redirects to the GRASP overview instead of delivering an original gzip in external inspection; local Orca source-page probe additionally had a TLS trust error (do **not** disable SSL verification). | Separate [OpenGWAS ieu-b-38](https://opengwas.io/datasets/ieu-b-38) confirms EUR, `HG19/GRCh37`, N=757,601, 7,088,083 SNPs, but downloading VCF requires JWT / ephemeral link. | **HOLD_ACCESS** — exact historical URL not a verified currently retrievable binary; metadata is insufficient |
| T2D | [T2DGGI 2024 EUR consortium file](https://www.diagram-consortium.org/downloads.html) listed as a ~365MB ZIP with official published MD5 `46f13bf073c5560d9b068b623229ff01`. | Web form includes mandatory agreement to the no-reidentification terms and data documents; no direct file URL or header was inspected. | **HOLD_TERMS** — no bypass, no downloaded archive |

The BMI corrected source is the June 25, 2018 replacement; the older, incomplete-N GIANT upload must never be substituted. The checked BMI gzip is **not a genome-wide full integrity check**. The absence of EAF in GIANT cohorts (HRS frequency is supplied) limits downstream allele/frequency QC, and 2018 meta-analysis includes UK Biobank, with potential overlapping participants in outcome GWAS. Do not mix HRS frequency with cohort EAF, GRCh37 with GRCh38, or SBP in mmHg with BMI SD-scaled effects.

Orca evidence: `scripts/ckm_stroke/probe_primary_file_links.py` found uniquely named GIANT source + README links; NHLBI source-page probe returned a **TLS certificate verification error**; T2D form-only page had no direct anchor and requires authorization. `scripts/ckm_stroke/probe_bmi_raw_header.py` independently requested 512 KiB range and stopped after **16 KiB**, with **10/10 header columns**, `technical_header_status=PASS_HEADER_ONLY`; offline tests passed **3/3**.

**Scientific release gate:** Complete verified full-file download and checksum (if a publisher checksum exists, compare it; otherwise record SHA256), complete allele/position/effect/N/P QC, define sample overlap handling, then run EUR-matched LDSC h2 and intercept diagnostics before rg. T2D needs explicit case/control counts, effective N, log-OR scale and permitted download. All 3 traits remain **CANDIDATE_ONLY**, and `public_sources.tsv` was deliberately not promoted.
