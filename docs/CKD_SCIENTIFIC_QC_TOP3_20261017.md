# CKD SCIENTIFIC QC GATES — SDCCAG8 / INHBC / UMOD

**Date:** 2026-10-17. **Status:** DOCUMENTED-SOURCE SCIENTIFIC QC, NOT INDEPENDENT RAW-LEVEL RE-EXECUTION.

## Evidence provenance

All numbers below were read directly from the dated (2026-10-05) Google Drive candidate snapshots. These snapshots explicitly cite server source tables `STAGE2B_MATCH_QC.tsv`, `STAGE2C_LD_QC.tsv`, `STAGE2B_COLOC_SUMMARY.tsv`, `STAGE2C_SUSIE_DEFAULT.tsv` and `STAGE2C_SUSIE_FAILURES.tsv`. **They have not been independently recalculated from raw signed-allele matrices.** The previously SHA-pinned integrated nine-gene server file was verified and Figures 1–4 reproducibly rendered, but exact priors/alleles and LD were not reassessed in that workflow.

Source snapshots:
- [SDCCAG8](https://drive.google.com/file/d/1yZTqWjNpRHsqiqYMjhW8EEnuQ0-wQldR/view)
- [INHBC](https://drive.google.com/file/d/1h3lu4PKqWmQWoOVhUURVNAUE0mQCc6dV/view)
- [UMOD](https://drive.google.com/file/d/1DVDG2RM7TRQ7fEK1gFicggQhkOEvxgwH/view)

## 1. Primary source metrics

| Source metric | SDCCAG8 | INHBC | UMOD |
|---|---:|---:|---:|
| Stage2B shared GWAS–pQTL variants | 6,269 | 4,777 | 6,048 |
| Stage2B pQTL match fraction | 0.5664 | 0.5990 | 0.6141 |
| Stage2C signed EUR LD SNPs | 3,854 | 2,679 | 3,090 |
| EUR reference donor size | 503 | 503 | 503 |
| Stage2C pQTL credible sets | 9 | 9 | 10 |
| Stage2C eGFR credible sets | 1 | 1 | **NA — nonconverged** |
| Stage2B ABF H4, default p12=1e-5 | 0.9794495841 | 0.9512705477 | 5.5397846e-41 |
| Stage2C SuSiE max H4 | 0.9794116279 | 3.4359021e-14 | **NA** |
| EUR classification | core_shared_signal | abf_susie_discordant | susie_unresolved |
| EAS eGFR stage1 anchor P | 0.08717 | 0.06445 | 2.2599e-47 |

**Major interpretive distinction:** The Stage2B overlap SNP count and Stage2C 1KG EUR LD SNP count are *different analysis sets*. Do not claim that coloc.abf and coloc.susie used exactly identical variants simply because their H4 estimates are displayed side by side.

## 2. ABF prior sensitivity (same 2026-10-05 source)

The ABF uses p1 = 1e-4, p2 = 1e-4 and changes the shared-causal variant prior p12:

| p12 | SDCCAG8 H4 | INHBC H4 | UMOD H4 |
|---|---:|---:|---:|
| 1e-6 | 0.8265720 | 0.6612635 | 5.5398e-42 |
| 1e-5 | 0.9794496 | 0.9512705 | 5.5398e-41 |
| 1e-4 | 0.9979062 | 0.9949035 | 5.5398e-40 |

**Conclusion:** INHBC's H4 is meaningfully prior-sensitive and its large ABF H4 does not survive multi-signal SuSiE corroboration. The UMOD single-variant ABF provides essentially no H4 support throughout the tested grid and shows H3≈1, but this should be described as **these molecular-region models favor distinct signals**, not a denial of established UMOD kidney biology.

## 3. Paper-specific claim gates

### SDCCAG8 — supportive shared EUR locus, not established causal protein
- ABF and reference-EUR LD SuSiE are numerically concordant (H4≈0.979) and the source describes `ABF_and_SuSiE_shared_signal`.
- Stage2C uses 1KG EUR n=503, not the actual molecular donor matched LD. Independent SNP/allele, priors, conditional signal and locus coverage should be reviewed before saying unequivocal shared causal variant.
- The EAS stage1 anchor P≈0.087 supports, at most, same-SNP effect-direction annotation. It does not establish independent East Asian replication.
- HPA macrophage-enriched expression context is not a macrophage-specific pQTL mechanism.
- **Allowed tentative wording:** "SDCCAG8 shows convergent EUR regional colocalization support under the evaluated models, subject to reference-LD and allele-harmonization limitations."

### INHBC — ABF vs SuSiE discordance is the finding
- Source ABF PP.H4≈0.951 at p12=1e-5 but SuSiE max H4≈3.44e-14 despite 9 pQTL and 1 eGFR credible sets.
- H4 at p12=1e-6 falls to ~0.661; this prior sensitivity must be visible in Supplementary.
- Independent kidney pQTL evidence was previously published by Hirohama 2025 and does NOT resolve internal molecular-region discrepant causal architectures. Novelty must be qualified.
- EAS stage1 anchor P≈0.064 is not replicated at nominal 0.05.
- **Allowed wording:** "Single-signal ABF suggested colocalization at INHBC, but LD-aware multi-signal analysis did not support a shared credible-set signal in the evaluated data."

### UMOD — strong locus association ≠ current tested plasma protein mediation
- EUR eGFR MR FDR≈1.62e-218 and EAS anchor association P≈2.26e-47 in dated snapshot; these do NOT override the molecular-region results.
- ABF H3≈1 and H4≈5.54e-41 at p12=1e-5; the pQTL signal and eGFR GWAS association may be distinct.
- Source `susie_rss` eGFR model did not converge in 1,000 iterations. **Do not enter H4=0 or call SuSiE negative.** Report `NA (nonconvergence)`.
- Strong kidney Loop of Henle localization from HPA is biological context, not proof of plasma protein causal mediation.
- **Allowed wording:** "UMOD has robust kidney-trait genetic associations; the tested plasma pQTL–eGFR shared-signal hypothesis was not supported by ABF under tested priors and multi-signal inference remained unresolved due to nonconvergence."

## 4. Mandatory scientific QC before paper-ready PASS

1. **Allele/genome build:** For each locus, verify hg19/hg38 transformations, reference/alternate and effect alleles in the original pQTL and eGFR files, strand/ambiguous exclusions and signed effect direction. Summarized direct/swapped counts are not enough.
2. **Region SNP set identity:** Export variant ID lists underlying ABF (6269/4777/6048) and SuSiE LD (3854/2679/3090), quantify intersection and material lead/sentinel coverage. Similar H4 across non-identical SNP sets is not direct method replication.
3. **LD QC:** Evaluate actual signed LD ordering/allele dosage correspondence and near-positive-semidefinite conditioning. EUR 503-reference sample size limits fine mapping reliability, especially complex loci.
4. **SuSiE convergence:** Capture warnings, credible sets, posterior inclusion probabilities and diagnostic reasons. Mark UMOD not estimable until a scientifically justified data/model change; no force-run to obtain a desired H4.
5. **Prior sensitivity:** Preserve three p12 values; expose H3 and H4 together. High H4 at p12=1e-4 alone must not count as independent confirmation.
6. **Sample overlap and molecular specificity:** Check EUR kidney outcome GWAS vs UKB-PPP overlap and effects of aptamer assay artifacts, LD and horizontal pleiotropy.
7. **EAS:** Separate nominal single-anchor results from EAS regional ABF/dual-LD and KoGES longitudinal validation. Controlled individual genotype is still unverified.
8. **External precedent:** Separate candidate support from novelty, especially INHBC/UMOD; compare with Si 2024 and Hirohama 2025.

## 5. Figure / manuscript instructions

- Figure 3: retain all nine source panels and the NA UMOD result. Existing Figure 3 is a **method-comparison descriptive plot** with distinct SNP sets, not demonstrated replication.
- Add **Supplementary prior-sensitivity Table** with the 3x3 values above, plus per-locus H3, overlap SNP count, LD SNP count, credible-set convergence.
- Figure 2: annotate effect units, instrument strength and effect allele, plus source overlap limitations.
- Figure 4: expression and kidney significant-only QTL counts do not validate causal cell-type mediation.
- The finished paper must never claim any candidate `causal gene VERIFIED` solely from figures or pass of synthetic CI.

**Outcome:** Documented-source QC produces distinct interpretation tiers. **RAW_SNP_LD_ALLELE_VERIFICATION = PENDING** for all three; this is not the same as failing scientific inference.
