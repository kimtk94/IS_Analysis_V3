# P7 — EMT gene-level robustness and microbiome metadata readiness

## Run provenance
- [GitHub Actions P7 3-way matrix #41179413138](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179413138). Started 2026-10-31.
- Code `scripts/ocular/ocular_p7_emt_and_microbiome_gates.py`.
- Workflow `.github/workflows/ocular-p7-emt-microbiome.yml` manual only after one-time PR branch launch. Mode-specific artifacts retain aggregate summaries/genes but no patient- or sample-identifying mappings.
- P7a EMT job: **result gate pending independent artifact review** at time of initial draft. Do not interpret it as PASS until the actual P7a result is appended in a follow-up commit.

## P7b — sample-specific human clinical/control metadata accessibility [PASS metadata audit, HOLD clinical analysis]

Job `independent-evidence-gates (samplemeta)` SUCCESS; [artifact #13034745113](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179413138/artifacts/13034745113). Live official ENA `read_run` fields `sample_title`, `sample_description`, sample/run IDs, layout and library strategy were accessible for every study. An explicit case/healthy/negative extraction blank keyword screen of **only sample-specific title and description** returned:

| Project | Runs / unique samples | Provisional disease classifiable | Provisional healthy classifiable | Explicit sequenced extraction/PCR blank named in ENA sample text |
|---|---:|---:|---:|---:|
| PRJNA1145861 (conjunctival 16S) | 139/139 | **0/139** | **0/139** | **0/139** |
| PRJEB74248 (ocular shotgun) | 62/62 | **0/62** | **0/62** | **0/62** |
| PRJNA1416399 (oral 16S) | 129/129 | **0/129** | **0/129** | **0/129** |

**Interpretation: critical metadata linkage gap.** ENA `sample_title` or `sample_description` alone could not reconstruct case/control assignments, exposure to ocular hypotensive drops, steroids, antibiotic exposure, or bilateral-eye/repeated-subject grouping. The absence of `blank` in these metadata fields **does not prove** negative controls were absent from wet-lab procedures. Published articles ([PMC12075858](https://pmc.ncbi.nlm.nih.gov/articles/PMC12075858/) and [PMC11172891](https://pmc.ncbi.nlm.nih.gov/articles/PMC11172891/)) discuss negative controls; the latter reports negative-control libraries did not pass sequencing QC. Need supplementary clinical participant tables and true run/sample crosswalks to reconstruct groups and adjust for medications/low biomass. Do not publish or invent private sample labels.

**Decision:** GO metadata audit execution, **HOLD** use for reliable glaucoma patient-vs-control differential taxa or microbiome–tear proteomics correlations.

## P7c — actual shotgun submission file-format manifest [PASS metadata formats, HOLD data accessibility]

Job `independent-evidence-gates (formats)` SUCCESS; [artifact #13034760063](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179413138/artifacts/13034760063).

For **PRJEB74248**:
- ENA **62/62** shotgun runs with original `submitted_ftp` file paths; of these **31** have `fastq_ftp`, the remaining **31** have `submitted_ftp` only.
- **155/155** originally submitted individual filename entries end with `.fastq.gz` (not BAM/CRAM/SRA).
- **155/155** individual files have an associated **32-character MD5 checksum entry** in ENA metadata.
- Sum of reported original submitted sizes: **528,925,176,960 bytes** (approximately 528.9 decimal GB; the figure is *submitted file bytes*, not processed read counts or unique sequenced biological patient samples).
- This pass **did not** download the data, inspect gzip/FASTQ integrity, perform even a limited binary retrieval, run host-depletion, verify paired R1/R2 completeness, classify taxa, or link tear-protein samples.

**Decision:** GO file-type metadata and checksum manifest feasibility. **HOLD** raw-file download QC (e.g. representative HTTP HEAD/Range, gzip checksum), correct R1/R2 pairs, read-quality checks, controls, host depletion, phenotype crosswalk, actual differential-microbiome inference. The full submitted 529GB-scale dataset is not appropriate to download in an unbounded CI run.

## Main research logic
The DEX–TM experiment tests a human cultured meshwork glucocorticoid response (paired donors, but treatment-specific 10x prep). Ocular surface and oral microbiome studies are *different samples from different populations*, not steroid-exposed matched TM donors. Without patient-level steroid exposure or IOP response, no microbial mediation or steroid-induced glaucoma patient classifier can be validated. Preserve **GO computational QC; HOLD dose/batch and microbial clinical links; STOP causal microbial glaucoma risk or clinical prediction**.
