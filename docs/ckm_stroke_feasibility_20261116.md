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

## Critical P0 discrepancy

The PR phenotype manifest maps **LAA → GCST90104542**, **CE → GCST90104541**, **SVD → GCST90104543**.
A separate 2026 published supplement (Table S2, https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs12975-026-01509-7/MediaObjects/12975_2026_1509_MOESM2_ESM.pdf)
reports **LAS → GCST90104541** and **CES → GCST90104542**. Its validity has not been independently adjudicated.
**Do not swap labels automatically**: cross-check against the *original* GIGASTROKE consortium release, GWAS Catalog authoritative study records, original summary-stat file headers and sample size. P0 must remain HOLD until resolved.

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

1. Resolve **GCST90104541/42 label conflict** against official original GWAS source; lock GRCh37, ancestry, phenotype/sample counts.
2. Confirm chr21/20/19 job outcomes, review uploaded QC artifacts; do not dispatch additional chromosomes on queued evidence.
3. Lock source candidates for EUR BMI/SBP/T2D; independently verify full release provenance, allele orientation and cohort overlap.
4. Build durable normalized summary-stat artifact handoff/ancestry matched reference LD. Only then begin LDSC h2/intercept/rg.
5. Ask KDCA/CODA for authorized phenotype/genotype/metabolomics overlap counts and access instructions, not unrestricted individual-level data.

The executable metadata audit in scripts/ckm_stroke/feasibility_gate.py produces a JSON status artifact from current repository manifests, with fail-closed scientific HOLD. It does not claim to resolve unavailable scientific evidence.
