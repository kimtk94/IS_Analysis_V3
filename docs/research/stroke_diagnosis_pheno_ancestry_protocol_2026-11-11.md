# Stroke disease definition, phenotyping, and EUR–EAS study direction (2026-11-11)

**Status:** literature/diagnostic-source audit; NOT a re-adjudication of GIGASTROKE cases and NOT a course-submission draft. The course template forbids generative-AI-written submissions.

## 1. Three levels that must not be conflated

| Level | Meaning | Evidence | What we do NOT know |
|---|---|---|---|
| Clinical disease definition | AHA/ASA tissue-based definition: overt neurological dysfunction attributable to focal cerebral/spinal/retinal infarction. For objective infarction evidence, symptoms need not last 24 h; the classical >=24 h clinical criterion is a fallback if imaging/pathology is inadequate, after excluding alternate causes. | Sacco et al. *Stroke* 2013 doi:10.1161/STR.0b013e318296aeca | Actual patient diagnoses in the GWAS |
| Etiological classification | TOAST assigns LAA, cardioembolism, small-vessel occlusion, other determined, and undetermined based on clinical examination, CT/MRI, vessel imaging, cardiac assessment and exclusions; lesion diameter and stenosis criteria are *context-dependent*, not sufficient alone. | Adams et al. *Stroke* 1993 doi:10.1161/01.STR.24.1.35; ISGC doi:10.1161/STROKEAHA.114.006839, AHA/ASA 2021 guidance | Whether all contributing GWAS cohorts used identical TOAST/CCS adjudication |
| Operational GWAS ascertainment | Consortium-wide summary statistics combine multiple recruiting designs/criteria. GIGASTROKE had 29 population-based studies/biobanks and 25 clinic-based case-control studies, with cohort-specific case ascertainment; paper delegates definitions to Supplementary Information. | Mishra et al., *Nature* 2022 doi:10.1038/s41586-022-05165-3 | Exact EAS per-cohort diagnostic algorithm, ICD lists, imaging confirmation, incident/prevalent, TOAST vs other; not derivable from the study accession alone |

A physician's clinical ischemic-stroke diagnosis, an ICD-10 I63 discharge code, and the label "ischemic stroke" on a GWAS Catalog record are **different evidentiary strengths**.

### Clinical diagnosis – review checklist (not our GWAS source inclusion claim)

1. Sudden focal neurological symptoms or signs consistent with a vascular territory; record onset time and vascular territory; distinguish stroke mimics (seizure, migraine, hypoglycemia, neoplasm, etc.).
2. Brain CT (commonly non-contrast) to detect hemorrhage/exclude mimics; MRI with DWI as appropriate to identify acute infarction. A normal early CT does not exclude ischemic stroke.
3. Objective focal infarction with symptoms supports IS regardless of whether symptoms last >=24 hours; if imaging/pathology inadequate, the >=24-hour clinical criterion may establish CNS infarction when competing etiologies are excluded.
4. If ischemic, investigate etiology with vascular imaging (CTA/MRA/ultrasound), rhythm monitoring/ECG and cardiac evaluation where clinically appropriate.
5. Do **not** count isolated TIA without infarction, silent infarction without symptoms, intracerebral hemorrhage, subarachnoid hemorrhage, or nonstroke mimics as clinically symptomatic IS by default.
6. NIHSS = acute severity, mRS = disability/outcome; neither determines IS case status alone.

References: AHA/ASA 2013 https://pmc.ncbi.nlm.nih.gov/articles/PMC11078537/ ; Korean Stroke Society clinical initial assessment https://stroke.or.kr/guidelines/view.php?dep1=2&page=1&sid=33 .

### Stroke subtype definitions – source-qualified

- **LAA/LAS:** infarction related to large-artery atherosclerosis with compatible clinical/vascular imaging; >50% relevant arterial stenosis supports classic TOAST, alternative sources considered. Not any patient who has an arterial plaque.
- **CE/CES:** plausible cardioembolic source (e.g., AF/intracardiac thrombus), appropriately evaluated competing large-artery sources. AF presence **alone** does not automatically establish causal CE for each event.
- **SVO/SVS (often called SVD):** compatible lacunar clinical syndrome and small deep infarct, excluding competing major embolic/large vessel causes. Classic TOAST imaging cutoff <1.5 cm is not a universal replacement for clinical adjudication.
- **Other determined** vs **undetermined:** distinguish established rare causes, multiple plausible causes, incomplete workup, and truly cryptogenic after adequate workup.
- **TIA:** transient focal ischemic symptoms **without evidence of acute infarction** on appropriate evaluation; do not use a pure 24-hour clock to override positive MRI evidence.

References: original TOAST https://doi.org/10.1161/01.STR.24.1.35 ; ISGC phenotype standard https://pmc.ncbi.nlm.nih.gov/articles/PMC4465378/ ; 2021 AHA/ASA https://doi.org/10.1161/STR.0000000000000375 .

### Administrative phenotyping

- ICD-10 I63 = cerebral infarction; I60 = subarachnoid hemorrhage; I61 = intracerebral hemorrhage; I64 = unspecified stroke; G45 = TIA.
- ICD codes alone do **not** guarantee imaging-confirmed adjudicated clinical stroke. UK Biobank administrative stroke definitions show ascertainment-dependent differences (https://pubmed.ncbi.nlm.nih.gov/35185750/).
- Record separately: imaging adjudicated clinical case, specialist adjudication, ICD/hospital record, self-report, other. No assumption that every GIGASTROKE contributor is imaging confirmed.

## 2. GIGASTROKE labels and current verification

**Meta-data verified by GitHub Actions #42620703387 (5/5 catalog traits), NOT source GWAS bytes or cohort-level diagnosis.**

| EAS phenotype | GWAS Catalog accession | Literature cases | Clinical meaning | Proposed analysis |
|---|---|---:|---|---|
| Any stroke (AS) | GCST90104544 | 27,413 | May include ischemic, hemorrhagic and unspecified stroke; not interchangeable with IS | AS-only sensitivity, not primary |
| Ischemic stroke (IS/AIS) | GCST90104545 | 19,032 | Clinically ischemic phenotype as operationalized across contributing studies; modality/criteria not uniformly proven | **PRIMARY** |
| Cardioembolic stroke (CES) | GCST90104546 | 926 | Etiologic subtype; low power | Exploration, not EAS primary |
| Large artery stroke (LAS) | GCST90104547 | 1,735 | Etiologic subtype; low power | Exploration |
| Small vessel stroke (SVS) | GCST90104548 | 5,532 | Etiologic subtype; phenotyping uncertainty must be audited | Secondary only |

**Avoid assumptions:** main consortium paper states study-specific diagnostic definitions are in its Supplementary Information. The appendix and per-cohort mapping to EAS analyses have **not** been systematically verified here. Do not claim all 19,032 EAS IS cases were MRI-confirmed or diagnosed using identical TOAST rules; do not claim that all controls had imaging. Labels are aggregate GWAS phenotypes. TOAST misclassification and case ascertainment heterogeneity are plausible contributors to observed population differences.

**EAS accession mapping:** catalog metadata verified, but literature case counts are not file-verified. Source https://doi.org/10.1038/s41586-022-05165-3 ; https://www.ebi.ac.uk/gwas/ .

## 2A. Newly located cohort-specific ascertainment evidence (preprint supplementary appendix)

Source: [Mishra et al. GIGASTROKE supplementary appendix, Research Square preprint version](https://assets-eu.researchsquare.com/files/rs-1175817/v1/609a86c10d4db8cca4fc7ca5.pdf). This is an **earlier preprint supplement**, not a reconciled final version-of-record table. It verifies that multiple ascertainment methods were described in the consortium, but the accession-specific membership/weight and final publication amendments remain to be verified.

| Cohort example | Stroke case identification | Subtyping | Phenotype classification implication |
|---|---|---|---|
| BioBank Japan (BBJ), Japan | Physician diagnosis at 66 collaborating hospitals; patients >40 years; medical-record survey. Appendix reports 16,256 IS registrants, including 1,256 LAS, 710 CE, 4,613 SVS | TOAST from medical records | Clinician-ascertained with etiologic subtypes; **do not conflate** appendix cohort count with EAS meta-analysis sample count |
| Hisayama / Fukuoka Stroke Registry, Japan | Sudden focal nonconvulsive neurological deficit >24 h due to ischemia; stroke neurologist diagnosis using CT/MRI, vascular and cardiac evaluation | Cerebrovascular Disease III + TOAST | Classical duration-defined phenotype; differs from tissue-based <24 h MRI-positive IS inclusion |
| China Kadoorie Biobank (CKB), China | Incident records through death/disease registries and health insurance; baseline stroke/TIA-free; recorded I60, I61, I63, I64 | No uniform TOAST subgroup confirmation in viewed description | AS encompasses multiple etiologies; cannot treat I60/I61/I64 as IS |
| Korean Cancer Prevention Study-II (KCPS2), Korea | Prospective incident stroke ascertainment using National Health Insurance claims and linked mortality/hospitalization data; appendix states broad ICD-10 I60–I69 | No individual clinical TOAST adjudication demonstrated by appendix | Administrative coding source; **do not** assume I60–I69 are all infarctions or all imaging-confirmed |
| Tohoku Medical Megabank (TMM), Japan | Appendix describes population cohorts but gives no detailed incident stroke diagnostic algorithm in viewed section | Unverified | Do not infer clinical case definition from cohort name |

Specific PDF appendix evidence: BBJ pp.14–15; Hisayama pp.24–25; CKB and KCPS2 pp.39–40; TMM p.43. A definitive answer requires reconciling these with the final Nature supplementary DOCX and mapping each contributing dataset to GCST90104544–48.

### Revised phenotype interpretation

- **Partial source review PASS:** verified preprint narrative examples show both clinical/TOAST and administrative/ICD ascertainment among contributing cohorts.
- **P0 cohort-to-accession phenotype completeness HOLD:** we have NOT established the exact participant contribution of BBJ, Hisayama, CKB, KCPS2 or TMM to each EAS accession, nor whether final supplemental definitions changed.
- **Causal interpretation limitation:** Differences between EAS and EUR genetic associations can reflect variation in phenotype ascertainment, case mix, linkage disequilibrium, allele frequency and sampling power. Do not attribute a difference solely to ancestry.
- **Design implication:** Primary outcome remains IS, but add a documented phenotype ascertainment table and prespecified sensitivity to clinical versus administrative cohorts **only if source-stratified GWAS summary statistics can actually be obtained**. Without cohort-stratified data, disclose unresolvable heterogeneity rather than promising correction.

## 3. Proposed scientific question and design

**Question:** *Among independently preselected ischemic-stroke susceptibility loci, how much does ancestry-matched blood cis-eQTL evidence alter locus-to-gene prioritization in EUR versus EAS, after accounting for phenotype ascertainment, statistical power, variant coverage and ancestry-matched LD?*

- **Primary outcome:** binary genetic susceptibility to **ischemic stroke**. Not response, prognosis, recurrence, mortality or 6-month progression.
- **Discovery:** EUR IS GWAS, preselect loci without inspecting EAS results. Provenance/independent recruitment and sample overlap must be recorded.
- **EAS evaluation:** EAS IS GWAS (GCST90104545); compare matched variant effect size and SE, sign concordance and uncertainty. A different P value is *not* proof of ancestry difference.
- **Regulatory interpretation:** whole-blood cis-eQTL summary results in EUR and EAS, at identical genomic loci with harmonized SNP set and appropriate matched LD. Within each ancestry evaluate evidence of shared GWAS/eQTL variant(s) by methods compatible with number of independent signals (e.g. SuSiE + coloc); report PP(H3), PP(H4), credible set and sensitivity to priors and SNP missingness. Gene rank change must be estimated on a **coverage-matched** subset.
- **Optional biology:** independent single-cell transcript expression of prioritized genes in vascular/endothelial/immune lineages. Expression presence alone **cannot** show that specific GWAS variants regulate those cells.
- **Planned negative/sensitivity checks:** exclude weaker ascertainment cohorts *only if cohort-stratified summary GWAS are available*; compare AS versus IS **as phenotype-sensitivity only** recognizing shared samples and hemorrhagic inclusion; assess power, INFO/MAF, imputation platform, sample overlap and prevalence differences.
- **Independent replication criterion:** not merely querying another phenotype from the *same* GIGASTROKE individuals. Requires documented distinct subject sets, or explicitly label the result 'within-consortium cross-ancestry assessment', not independent external validation.

## 4. Novelty falsification and scope discipline

A 2026 EAS blood eQTL paper **already reports** that ancestry-matched QTL resources improve gene prioritization for EAS complex traits (MAEEA, n=2,024, https://doi.org/10.1093/nsr/nwag566). Thus that **generic** claim is not novel. A specific candidate gap remains **stroke-locus-level robustness under phenotyping and coverage matching**, provided a targeted literature matrix shows earlier stroke work did not already do it. Earlier stroke subtype GWAS/sc-eQTL colocalization (2024; https://pubmed.ncbi.nlm.nih.gov/38744227/) is directly competing prior work.

Additional directly competing publication: Surakka et al., *Cell Genomics* 2023, https://doi.org/10.1016/j.xgen.2023.100345, combined 16 biobanks (GBMI) and MEGASTROKE, used GBMI PheCode 433.21 (cerebral artery occlusion with cerebral infarction), identified heterogeneity at ABO across ancestries and a sex/EAS-specific ALDH2 rs671 signal, and conducted colocalization at several IS loci. Therefore **ABO / ALDH2 sex × ancestry, generic multi-ancestry IS, and generic GWAS–eQTL colocalization are pre-existing results, NOT novel hypotheses.** The biobank-based PheCode ascertainment also illustrates why claims of phenotype equivalence across clinical vs EHR studies require evidence.

This is a **conditional GO**, not a novel result. If exact cohort-level phenotype ascertainment cannot be reconstructed, explicitly limit inference to *published GWAS-level phenotype comparison* and do NOT attribute contrast to ethnicity or biology.

## 5. Gates and priority next actions

| Gate | Concrete evidence | Status |
|---|---|---|
| P0 Phenotype map | Preprint appendix checked for BBJ, Hisayama, CKB, KCPS2 and TMM; verify final Supplementary DOCX and map each contributing cohort/ascertainment to EAS accession/outcome | **PARTIAL / HOLD** |
| P1 Source files | Full IS EUR and EAS GWAS download, headers/alleles/build/INFO, hashes, actual sample counts and study provenance | **HOLD** |
| P2 EAS eQTL | MAEEA downloadable full cis-eQTL summary variants/effect/SE/MAF/N and distribution rights, matched EUR comparator | **HOLD** |
| P3 Statistical feasibility | Allele harmonization/coverage and ancestry-appropriate LD; credible-set overlap and detectable-effect analysis, prospectively selected discovery loci | **HOLD** |
| P4 Prior-work novelty | Locus-by-locus prior colocalization gene/variant matrix; quantify genuinely new interpretation | **HOLD** |

Decision tree: P0/P1/P2 pass + adequate P3 => conduct primary EUR/EAS IS analysis; otherwise limit to robust **EUR IS regulatory interpretation + EAS directional follow-up** or seek another phenotype-specific angle. EAS CES n=926 **never** serve as headline powered EAS confirmation without a prespecified power analysis.

**Source status:** Guidelines and consortium publication verified; individual-level GWAS diagnoses not re-adjudicated. GWAS Catalog accession is administrative evidence of a reported phenotype, not a clinical diagnostic audit.
