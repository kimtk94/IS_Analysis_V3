# CKM → ischemic stroke subtypes: Actions-first study design

Study question: Which CKM-related metabolic and kidney genetic components
are shared across or specific to LAA, CE, and SVD ischemic stroke?

This is not a CKM Stage 4 to stroke causal analysis: stroke is itself an
end-stage clinical CKM criterion. Genetic correlation, MR, protein effects,
colocalization and tissue expression are distinct levels of evidence.

| Stage | GitHub Actions task | Advancement requirement |
|---|---|---|
| P0 | Phenotype/GWAS manifest | Source accession, build, ancestry and TOAST/ICD phenotype review |
| P1 | GWAS ingestion and streaming QC | Actual full normalized data, checksum, overlap, allele provenance |
| P2 | LDSC and optional Genomic SEM | Matched LD, h2/rg/intercept, formal stroke subtype heterogeneity |
| P3 | MAGMA/fine-mapping | Reviewed genes/loci, annotations, LD and FDR |
| P4 | cis-pQTL/eQTL MR and colocalization | Instrument quality, sensitivity, PP.H4, FDR and replication |
| P5 | GTEx/brain eQTL/snRNA | Donor-aware cell annotation and independent tissue plausibility |
| P6 | EUR/EAS comparison | Harmonized effect alleles, independent cohorts and precision |

A technical PASS means only that that gate ran correctly. The separate
scientific_status field remains HOLD/NOT_RUN/REVIEW_REQUIRED until
publication-quality evidence is available. A missing file, a missing
approved analysis command or absent official source mapping is HOLD,
not a simulated finding. P0 retains REVIEW_REQUIRED placeholders.

P1 expects a real, tab-delimited file with columns
CHR POS EA NEA BETA SE P. Bounded scans (--max-rows) return HOLD.
P2-P5 read a structured config/ckm_stroke/tasks.json plan containing
argv, inputs and outputs arrays and require explicit --execute. Commands
are launched with no shell. Unconfigured plans do NOT run.
P6 requires a real data/ckm_stroke/replication.csv with columns
locus,beta_eur,se_eur,beta_eas,se_eas,harmonized,independent.
Only the synthetic fixtures are included in Git.

The workflow produces seven JSON evidence artifacts. GitHub-hosted
runners are ephemeral and will not retain original GWAS data. Do not
commit summary statistics, patient-level data or tokens to this repo.
Before enabling P1-P6 real science, implement source downloads plus
independent hash/build/allele checks and persistent artifact storage.


## 2026-11-11 P1 real-source evidence

Official bounded gzip-prefix probe verified:
- GLGC EUR LDL-C, HDL-C, TG: tab-delimited, CHROM/POS_b37/REF/ALT/EFFECT_SIZE/SE/pvalue. ALT is the effect allele; TG is log-transformed and inverse-normalized.
- CKDGen EUR eGFR and UACR: **whitespace-delimited** (not tab), Chr/Pos_b37/Allele1/Allele2/Effect/StdErr/P-value. Allele1 is the effect allele.
- Header observation is not file integrity or full-row QC.

The source-specific normalizer now reads the CKDGen whitespace layout.
Excluded non-autosomal/non-SNP/missing/P-underflow rows are explicitly counted,
and any exclusion leaves a QC review HOLD rather than silently promoting a complete
analysis. Source-level checksums are needed before interpreting any output.

Dedicated workflow: .github/workflows/ckm-stroke-uacr-ingest.yml
- Downloads the entire **public** UACR EUR CKDGen Teumer2019 file within a 650 MiB safety cap.
- Verifies official consortium MD5 51270db6886892199975b343452f1d90
  from https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/MD5SUM .
- Reads the entire gzip stream (CRC during normalization), emits source and
  normalized SHA256, raw/retained row counts and exclusion breakdown.
- Uploads JSON metadata only; both raw and normalized GWAS remain ephemeral.
- A successful source-integrity pilot is *one validated file*, not P1 cohort-wide
  completion, clinical phenotype resolution, independent replication, or MR evidence.

The next scientific gate remains the full primary EUR phenotype manifest, including
BMI/SBP/T2D, and robust cross-study build/allele/overlap validation.


## Candidate source decisions (not yet production-locked)

| Trait | Potential source | Caution / next proof |
|---|---|---|
| BMI EUR | GIANT/Yengo 2018 updated GIANT+UKB meta-analysis | Consortium corrected an earlier sample-size error; confirm corrected version, source build, allele semantics and sample overlap with GIGASTROKE |
| SBP EUR | ICBP Evangelou 2018 / OpenGWAS ieu-b-38 | Metadata lists European n=757,601 and GRCh37, but bulk access and complete harmonized columns are not yet verified |
| T2D EUR | DIAMANTE/T2DGGI Mahajan 2022 or Suzuki 2024 ancestry-specific EUR | Website has click-through terms. Do not bypass consent or infer link from a study page; select study version with clear phenotype, ancestry and overlap |
| CKDGen UACR EUR | Teumer 2019 European subset | Official MD5 available; dedicated full-source pilot handles checksum and gzip, not overall P1 completion |
| GLGC lipids EUR | Graham 2021 ancestry-specific | Each raw EUR compressed file ~2.1 GB; do not attempt full download in the small bounded-source probe |

References:
- https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files
- https://gwas.mrcieu.ac.uk/datasets/ieu-b-38/
- https://www.diagram-consortium.org/downloads.html
- https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/README.txt
- https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/MD5SUM

These rows are *discovery candidates*, not verified and locked manifest entries.
