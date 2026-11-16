# CKM → IS feasibility decision | 2026-11-16

**Scope:** Scientific readiness triage of the existing Draft PR #284, **not** a claim of completed LDSC/MR.
This audit is deliberately metadata-only. It does not download original GWAS data or patient-level data.

## Verified project state

- P0–P6 pipeline is already implemented as staged evidence gates in PR #284; issue #285 is the central queue.
- eGFR and UACR full original-source normalization/variant technical QC have been reported in existing Actions evidence. UACR official checksum matched; eGFR official independent checksum still requires review.
- GLGC Graham2021 LDL-C, HDL-C and logTG **chromosome 22 only** technical QC passed (run #43431051168). Do not infer autosome-wide completion.
- A strictly gated chr21, chr20, chr19 expansion was launched as run #43848671279 on 2026-11-16; at observation it was QUEUED. This audit does not assume it passed.
- BMI, SBP, and T2D EUR primary sources remain unselected. Available source URLs are not verified full normalized GWAS files.
- tasks.json leaves P2–P5 argv/inputs/outputs empty. Thus **LDSC, MAGMA, MR/coloc and tissue-level science are NOT_RUN**.

## P0 accession/build source metadata gate: RESOLVED 2026-11-16

Directly downloaded EBI GWAS Catalog original metadata YAML and md5sum.txt for the four European outcomes. The originally recorded phenotype labels were correct. The secondary Song et al 2026 supplement reverses LAA/CE, and must not override the EBI publisher record.

| Outcome | Official accession | Original reported phenotype | Assembly | EUR total N | Published raw MD5 |
| --- | --- | --- | --- | ---: | --- |
| AIS | GCST90104540 | Ischemic stroke | GRCh37 | 1296908 | 1bfe4ae8a5042fb24cb0a562b05b0d2f |
| LAA | GCST90104542 | Large artery stroke | GRCh37 | 1241207 | 71ea9fd7c06791a831861053b6a1cfbe |
| CE | GCST90104541 | Cardioembolic stroke | GRCh37 | 1245612 | c252262586a6ea84b985f30dd7ab32be |
| SVD | GCST90104543 | Small vessel stroke | GRCh37 | 1241619 | f41bc5f55b3cab238e305e3e86363f69 |

Machine-readable evidence: config/ckm_stroke/gigastroke_eur_official.tsv. EBI source URL is recorded per study.
Original GWAS binary has not been downloaded or hash-checked in this gate: checksum/source-integrity scientific HOLD remains. Do not infer LDSC or MR result PASS.

## Novelty audit

Lei et al., 2026, *Journal of Cerebral Blood Flow & Metabolism*,
DOI: https://doi.org/10.1177/0271678X261462718 (online 2026-06-15)
already analyzed CKM metabolic signatures, incident stroke, MR and colocalization (LDLR, CELSR2).
Therefore 'CKM→any stroke with MR+coloc' alone is **not defensibly novel**.
The differential hypothesis should be **LAA/CE/SVD subtype heterogeneity, ancestry matched EUR/EAS genetic architecture, and independently replicated tissue/cell evidence**. Novelty is conditional on a formal up-to-date competitor table.

## KoGES independent parallel gate

KDCA KoGES official distribution lists Ansan–Ansung baseline **10,030** and F10 (2021–22) **5,511** observations. These are survey counts, **not** measured CKM transition event counts or genotype–metabolomics overlap. Data are available through CODA/OPEN KoGES following relevant review/approval.
2026 cohort update additionally references F11 data collection through 2024: verify distribution release availability before claiming 20/22/24-year coverage.

Sources:
- https://www.nih.go.kr/ko/main/contents.do?menuNo=300919
- https://coda.nih.go.kr/opkgs/
- https://doi.org/10.1093/ije/dyag192

KoGES scientific-status **HOLD** pending access rights, precise data use conditions, individual-level transition counts, attrition, genotype–metabolomics intersection and event power.

## Next actions in order

1. [RESOLVED metadata] Verify downloaded original GWAS file hash and variant columns separately; accession and GRCh37 labels locked.
2. Confirm chr21/20/19 job outcomes, review uploaded QC artifacts; do not dispatch additional chromosomes on queued evidence.
3. Lock source candidates for EUR BMI/SBP/T2D; independently verify full release provenance, allele orientation and cohort overlap.
4. Build durable normalized summary-stat artifact handoff/ancestry matched reference LD. Only then begin LDSC h2/intercept/rg.
5. Ask KDCA/CODA for authorized phenotype/genotype/metabolomics overlap counts and access instructions, not unrestricted individual-level data.

The executable metadata audit in scripts/ckm_stroke/feasibility_gate.py produces a JSON status artifact from current repository manifests, with fail-closed scientific HOLD. It does not claim to resolve unavailable scientific evidence.

## Additional near-identical prior art found in the 2026-11-16 audit

1. **Proteome-wide Genetic Study in East Asians and Europeans Identified Multiple Therapeutic Targets for Ischemic Stroke** (2025; PubMed 40304040): EAS/EUR MR and colocalization for FGF5, F11, ALDH2, ABO, GRK5, PROCR, MMP12 etc. This overlaps earlier user IS candidate screening and refutes novelty from EAS+EUR proteomics alone.
   https://pubmed.ncbi.nlm.nih.gov/40304040/
2. **Song et al., Proteome-Wide Mendelian Randomization Identifies Therapeutic Targets for Ischemic Stroke** (*Translational Stroke Research*, 2026-10-22, DOI 10.1007/s12975-026-01509-7): 589 INTERVAL pQTL instruments versus AIS, LAA, CES and SVD; deCODE replication, colocalization and 65,479 brain vascular nuclei. Highlights SPATA20/KDELC2. **Subtype-wide protein MR+coloc is already published**.
   https://pubmed.ncbi.nlm.nih.gov/42911303/
3. **Lu et al., Genomic Structural Equation Modeling Reveals CKM Syndrome Genetic Architecture** (*Journal of Diabetes*, 2026-04-13, DOI 10.1111/1753-0407.70225): existing CKM genomic SEM, LDSC, GWAS-by-subtraction, fine-mapping, TWAS and cell annotations.
   https://doi.org/10.1111/1753-0407.70225
4. **Multivariate genomic analysis of 5 million people elucidates the genetic architecture of shared components of metabolic syndrome** (*Nature Genetics*, 2024): a large pre-existing multivariate MetS genomic architecture baseline.
   https://www.nature.com/articles/s41588-024-01933-1

**Revised novelty decision:** Merely extending earlier IS proteomics to subtype MR, merely EAS/EUR triangulation, and merely CKM Genomic SEM are insufficient. The key potential contribution must be a *formally tested interaction between CKM component-specific liability and stroke subtype heterogeneity*, with population-specific genetic interpretation and independent biological validation. If the formal contrast is weak, STOP the claims and consider KoGES longitudinal transition analysis (subject to data access).
