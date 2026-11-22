# MASTER DEGREE — 8-track open evaluation before final 1–2 topic selection

**Decision revision:** 2026-11-02 KST. **No additional research project is disqualified yet.** Supersedes the earlier blanket "pause the other six" recommendation. Glaucoma has the highest **current manuscript readiness**; this does **not** prove all other tracks have lower scientific potential. Select 1–2 topics **only after** comparable bounded real-data and novelty experiments.

## What was unfair in the earlier shortlist?

- We contrasted a completed glaucoma manuscript with **other projects at P0 source acquisition** and treated lower engineering maturity as lower scientific value.
- HCC/Radiomics has actual original MRI and omics analyses; OSCC has a public 86 MB processed expression matrix; skin has a publicly deposited ~1 GiB full GWAS file. These are legitimate feasible source paths, not reasons to stop.
- GitHub CI success on a study's metadata or website HTML is not equivalent to actual scientific source access or new inference. Conversely, one HOLD metadata gate is not a negative biological result.
- IS's prior-art collision is important and is an originality constraint, not proof that a new locus/ancestry/technical-bias question cannot be formulated.

## Equal evaluation criteria (no fake scientific %)

Use 5 ordinal gates for **each** project, with verifiable outputs and no automatic ratio-to-final score until all are audited:

1. **Data:** direct proof of public original source, executable read, assembly/schema/metadata, sample units and access permissions.
2. **Identity:** unique original biological units, treatment/clinical labels and across-assay joins, without counting technical replicates as independent patients.
3. **Inference:** one prespecified falsifiable primary analysis with real data, appropriate independent N, sign/effect/FDR and leakage/batch/LD controls.
4. **External:** actually independent validation data or a defensible methodological reproducibility endpoint. Neither a same-cohort parameter perturbation nor cohort overlap counts as external replication.
5. **Novelty+thesis-fit:** claim-by-claim checked against 2022–2026 prior literature, clinical relevance and realistic manuscript scope; advisor approval is a distinct real-world gate.

Gate statuses: **PASS real scientific task**, **PARTIAL source/engineering**, **HOLD not verified**, **STOP proven impossible or substantially duplicate with no narrower contribution**. Do not promote a project to STOP just because it has fewer GitHub Actions or lacks a specific dataset at the first attempt.

## Proposed discriminatory experiment for each active topic

| Project | Next minimal **real-data** experiment | Pass criterion (not automatically thesis-ready) | Continue/HOLD/STOP decision trigger |
|---|---|---|---|
| Glaucoma DEX/TM | Complete existing 74/59/57/55 source-unit flow; integrate effect cutoffs and the 18-library rank diagram in manuscript | Actual figures/references checked; maintain strict drug-vs-library nonidentifiability warning | Already **GO for bounded computational thesis draft**; compare its novelty ceiling with other fields later |
| IS stroke pQTL/eQTL | From official EAS AIS original 6,725,216 rows and authentic exposure QTL, reproduce a **pre-specified locus-specific MR/coloc or allele/LD counterexample** and compare with 2025–26 proteome-wide MR papers | True original allele/SE/N/LD model and external-sample independence; one candidate finding not simply repeating named prior targets | **HOLD scientific inference** until QC passes, not STOP because prior art exists |
| HCC tumor RNA–Protein | Obtain **one verified official same-patient RNA–Proteome aliquot-linked cohort** (Gao PDC+OEP or another independent cohort) and one patient-level gene/pathway effect, not JS route hits | Quant original file, unique case/sample/aliquot matching, tissue/assay and original patient IDs verified; independent validation relative to Suter | If legitimate source inaccessible, consider a **methods/reproducibility HCC thesis** using authenticated Suter matrices before STOP |
| Radiomics ± RNA | Verify independent external MRI 59/59 actual image-mask extraction and patient independence from 14 registered PV test; then pre-specify reader/acquisition/segmentation robustness primary end point | Patient-level CI, independent external source and consistent PyRadiomics/IBSI parameters; *MRI–RNA remains separate HOLD* until matching | Could be a **standalone radiomics methodological thesis** even if MRI–RNA fails. Do NOT force merge with HCC unless identical patients/modalities |
| Ophthalmology P0 (outside original GSE313736) | Public PRJNA1243143 46-run accession–strain–treatment–timepoint map for matched 48h TM; then read original count data only where attribution verified | Complete unique run pairing and independent biological strain labels; genuine external inferential gene effects | Could strengthen glaucoma or become distinct ocular multiomics only if new disease/assay/exposure question justified |
| Muscle strength GWAS | Find **authenticated original 20 resistance-training loci with study IDs, variant IDs, trait/N and p**; validate genome-wide vs suggestive level and 1000G EAS/GTEx overlap for one locus | Valid study-to-SNP input, chromosome/build/alleles, outcome phenotype, ancestry LD | Missing upstream 20-study input -> HOLD rather than falsely call all variants GTEx-negative |
| OSCC oral cancer | **Download actual 86MB GSE103322 processed matrix**, hash and inspect matrix orientation/shape; connect original patient/cell annotations, then patient-aware one predeclared hypothesis and independently accessible bulk cohort | Verified patient pseudobulk ≥2 biological cohorts with real cohort IDs and independent test; raw FASTQs aren't needed for all analyses | If no reliable patient mapping, HOLD. If patient mapping real, upgrade beyond prior simple GEO accession-only gate |
| Skin aging MR | **Inspect original Bristol ~1GB perceived-age GWAS header and numeric rows**, then full allele/build/phenotype QC; identify independent exposure instruments and distinct outcome GWAS | Exact effect alleles, LOG_OR/LOG_SE, true coded phenotype and N, independent validated MR data, pleiotropy sensitivity | Perceived age may not validly measure *skin-specific photoaging*; absence of appropriate phenotype -> change question, not force a false cosmetic claim |

## Actual first discrimination CI result (2026-11-02)

[Run #41440909058](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41440909058) completed **3/3 jobs SUCCESS** (fixture + real source skin + real source oral). Note that job success is **not** scientific PASS, since source-level HOLD is recorded inside result JSON when blocked.

| Track | Actual observed outcome | Next decisive gate |
|---|---|---|
| Skin | **35,497 complete authentic GWAS prefix records sampled; 35,497/35,497 pass basic allele/frequency/LOG_OR/LOG_SE/P numerical checks** (4MB bounded prefix; 4,590 palindromic allele pairs, orientation uncertainty; header has SNP, CHR, BP, ALLELE1/0, A1FREQ, INFO, P_BOLT_LMM_INF, LOG_OR, LOG_SE, N_0/N_0.5/N_1, N_EFFECTIVE). This confirms an independently retrievable **real numeric dataset**, not just a readme; whole 1GiB file not fetched. | Original whole-file hash, genome build and logged perceived-age case/intermediate phenotype coding, MR target phenotype validation, valid independent genetic instruments and pleiotropy controls. |
| Oral | Direct **GEO FTP GSE103322_HNSCC_all_data.txt.gz returns HTTP 403** from this GitHub runner; scientific status stays HOLD. This is an access-route failure, **not evidence that processed tumor scRNA is absent**. | Try the officially documented Bioconductor ExperimentHub `EH5419` / GSE103322 SingleCellExperiment data through R/Bioconductor, authenticate 18-patient/cell biological-unit mapping, then choose a patient-aware test; raw FASTQ need not be required. |

**Decision consequence:** both remain candidates. In particular, skin's measured source status is materially stronger than last week's readme-only label, while OSCC requires a different public access route and must not be retired due to one runner's 403. The actual 2026-11-02 job results can be retrieved from per-track artifact JSON; the job fixture does not justify biological or MR promotion.

## Actual underway

- [GitHub CI real oral/skin source decision](https://github.com/kimtk94/IS_Analysis_V3/actions/workflows/research-decision-oral-skin-real-data-20261102.yml): distinct real-data jobs for full GEO processed matrix stream (bounded 125 MB compressed), and 4 MB authentic Bristol numeric GWAS prefix (Range, no large full download). All external access failures are **HOLD** with cause, not biological failure. No raw cell identities, genotypes, variant IDs or full expression matrices are exported to CI artifacts.
- Prior genuine external-data evidence for radiomics: [14 original registered-patient PV factorial CI #41208617950](https://github.com/kimtk94/Codex/actions/runs/41208617950). Claims about external 59/59 image success must be independently reconciled to actual workflow/run and source identity before being called 59 truly independent eligible patients.
- Official oral original processed expression: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322 ; GEO notes *raw sequencing* withheld due to privacy but **processed GSE103322_HNSCC_all_data.txt.gz 86 MB is public**. Bioconductor ExperimentHub independently provides a 5,902-cell/18-patient SingleCellExperiment wrapper: https://bioconductor.org/packages/release/data/experiment/html/GSE103322.html
- Official skin dataset/readme: https://data.bris.ac.uk/data/dataset/21crwsnj4xwjm2g4qi8chathha and https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt. Study effect scale is described as `LOG_OR`; perceived-age questionnaire does **not** automatically validate a skin-aging clinical trait.

## How to select when evidence is comparable

- **No final 1–2 choice today.** Preserve current 8-track daily scheduled checks; do not disable or rewrite automation without approval.
- Run one decisive actual analysis for every understudied track before ranking; judge HCC and Radiomics separately because either can qualify for a **methods thesis** without genuine MRI–RNA matched samples.
- After gate results, rank **defensible result novelty × independent source integrity × completion feasibility**, not number of commits/CI jobs or work already invested.
- To avoid excessive work, focus on one newly falsifiable experiment per candidate; stop open-ended API/source scraping if access is repeatedly blocked and no independent question remains.
- If any new project becomes more defensible than glaucoma, switch priority; glaucoma is currently an available fallback, not irrevocably the final thesis.

**Previous shortlist:** [2026-11-02 first pass](MASTER_DEGREE_FINAL_8_TO_2_SCIENTIFIC_SHORTLIST_20261102.md) is superseded for PAUSE/STOP recommendations; its numerical evidence and true prior-art overlap notes remain useful. No scientific inference was invalidated by this priority reversal.
