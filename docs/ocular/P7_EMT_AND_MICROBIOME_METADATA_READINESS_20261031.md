# P7 — EMT gene-level robustness and microbiome metadata readiness

## Run provenance
- [GitHub Actions P7 3-way matrix #41179413138](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179413138). Started 2026-10-31.
- Code `scripts/ocular/ocular_p7_emt_and_microbiome_gates.py`.
- Workflow `.github/workflows/ocular-p7-emt-microbiome.yml` manual only after one-time PR branch launch. Mode-specific artifacts retain aggregate summaries/genes but no patient- or sample-identifying mappings.
- **P7a EMT job SUCCESS**, verified [artifact #13033914032](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179413138/artifacts/13033914032). Public artifact ZIP downloaded and 194 EMT gene rows, 8 leave-top-effect tests, 24 candidate/marker rows directly inspected.

## P7a — 194 EMT genes and original GC/TM marker panel [PASS computational, HOLD mechanism]

P7a directly re-read official 2026 100nM 48h TM bulk supplemental workbook and original 16,050-gene P5b 500nM BeamA edgeR aggregate; **13,313** external 48h genes, **12,251** shared Ensembl IDs, **12,250** shared symbols. Official MSigDB 2025.1.Hs Hallmark gene sets with pinned SHA256 were reused. **194** shared EMT genes, of which **141/194 (72.68%)** have matching 500nM vs 100nM signs. For context **6,810/12,056 (56.49%)** non-EMT genes also agree in sign (not independent gene observations or formal significance). Among EMT genes, **49** satisfy BH-FDR<0.05 in both studies; **46 agree in sign, 3 disagree** (**SDC1, ITGA2, FAS**). Post-hoc pathways were not preregistered; no formal independent pathway p-value was calculated.

| Leave strongest EMT genes (by max abs log2FC) | BeamA 500nM EMT ES | Bulk 100nM EMT ES | Direction |
|---:|---:|---:|---|
| 0 | -0.478209 | -0.480715 | same (-/-) |
| 5 | -0.472174 | -0.441968 | same (-/-) |
| 10 | -0.437361 | -0.399887 | same (-/-) |
| 20 | -0.400991 | -0.330167 | same (-/-) |

The **UPR contrast remains discordant** after removing its 0/5/10/20 highest-effect pathway genes: BeamA ES remains negative (-0.344 to -0.304), external bulk ES positive (+0.288 to +0.243). This strengthens the **pathway-priority interpretation**, not causal dose attribution.

Positive-control and target interpretation:

| Marker | 500nM BeamA log2FC | 100nM bulk TM log2FC | Gene-level significance |
|---|---:|---:|---|
| **MYOC** | **+4.025688** | **+1.859591** | FDR<0.05 both; known GC-responsive TM gene, but not a formal batch-control experiment |
| **TGFBI** | -1.869208 | -0.842948 | FDR<0.05 both |
| **MMP2** | -0.589067 | -0.357295 | FDR<0.05 both |
| **ACTA2** | +0.469496 | +0.373745 | FDR<0.05 both |
| **SDC1** | -0.553539 | +0.631467 | FDR<0.05 both, opposite direction |
| **FN1** | +0.506583 | +0.152135 | external FDR .43: NOT independently significant |
| **COL1A1** | -0.814562 | +0.029319 | external FDR .95: NOT independently significant |
| **CCN2** | -0.537661 | +0.065488 | external FDR .89: NOT independently significant |

The MYOC cross-study concordance is a useful empirical positive control that both assays detected a GC/TM response; it **does not** independently prove randomized treatment-vs-library batch assignment. Thus **GO** EMT/ECM ranked-gene hypothesis prioritization, **HOLD** coherent TM ECM function/aqueous outflow mechanism until randomized and orthogonal studies, **STOP** clinical steroid-induced IOP prediction.

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

## P7d — official publication supplement vs 62 ENA aliases: computational audit PASS, clinical linkage HOLD

The original [P7d run #41179721638](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41179721638) FAILED with `TimeoutError` while reading the Europe PMC published supplemental ZIP, **not** with a sample-discovery conclusion.

We corrected the download logic to attempt a direct PMC supplement path first, fall back to Europe PMC with bounded transfer time, and write an **explicit HOLD** result rather than falsely passing sample mapping. [Corrected run #41180050136](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41180050136) **SUCCESS as an audit execution**, [artifact #13034413247](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41180050136/artifacts/13034413247) directly inspected:

- `status=HOLD_OFFICIAL_SUPPLEMENT_REMOTE_ZIP_UNAVAILABLE`
- Direct PMC route returned **non-ZIP response**; Europe PMC fallback **timed out** reading the archive.
- All 62 ENA `sample_alias + sample_title` entries remain unclassifiable into cases and controls by explicit keywords (**0 case / 0 healthy / 62 unknown**).
- No verified sample-to-patient crosswalk was obtained. The lack of a retrieved ZIP **does not prove** no supplementary clinical crosswalk exists.
- The original publication describes **32 human donors (16 glaucoma, 16 control)**, compared with **62 ENA runs**. Never equate run count with clinical N.
- Decision: **HOLD** independent differential taxa, paired tear proteome correlations, and clinical drug covariate analyses until exact sample↔patient mapping and low-biomass control metadata are obtained by a verified public or authorized source.

## Main research logic
The DEX–TM experiment tests a human cultured meshwork glucocorticoid response (paired donors, but treatment-specific 10x prep). Ocular surface and oral microbiome studies are *different samples from different populations*, not steroid-exposed matched TM donors. Without patient-level steroid exposure or IOP response, no microbial mediation or steroid-induced glaucoma patient classifier can be validated. Preserve **GO computational QC; HOLD dose/batch and microbial clinical links; STOP causal microbial glaucoma risk or clinical prediction**.
