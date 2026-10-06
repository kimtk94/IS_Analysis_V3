# SKIN_BEAUTY_MR_MASTER

**Status:** Research-design + executable Stage 0 draft v0.2
**Date:** 2026-09-27
**Primary theme:** Proteogenomic causal inference for objective facial skin-aging phenotypes in East Asian populations

---

## 1. Executive decision

The most promising skin/beauty MR direction is **not another generic UK Biobank perceived-facial-aging MR**. That outcome has already been reused across obesity, glycemic traits, metabolites, micronutrients, hormones, microbiome, telomere, lipid multi-omics, and 2026 proteogenomic studies.

The strongest immediately executable opportunity identified in this review is:

> **Proteogenomic dissection of objective site-specific facial wrinkling in Korean women using cis-pQTL Mendelian randomization, ancestry-aware colocalization/fine-mapping, and cross-platform / cross-ancestry validation.**

Primary outcomes are Korean facial-wrinkle GWAS accessions from the 2022 JID study:

- **Under-eye wrinkles:** GCST90094903 / OpenGWAS `ebi-a-GCST90094903`
- **Crow's feet:** GCST90094904 / OpenGWAS `ebi-a-GCST90094904`

**2026-09-27 availability re-audit:** the paper states that full summary statistics were deposited in the NHGRI-EBI GWAS Catalog, but the current EBI study page for GCST90094903 displays “Full Summary Statistics: Not available.” OpenGWAS nevertheless currently indexes both studies as complete `ebi-a` datasets with **N=11,079 and 4,939,316 variants each**. Therefore the operational source for Stage 0/1 is OpenGWAS, with the EBI accession retained as provenance. The total paper cohort of 17,019 includes discovery + replication; the OpenGWAS complete GWAS dataset corresponds to the discovery GWAS N=11,079.

These traits were objectively quantified and come from a study co-led by Hong-Hee Won. This gives a substantially stronger thesis narrative than using only self-/socially-perceived facial age.

A targeted literature search for the exact GWAS accessions and the terms “under-eye wrinkles Mendelian randomization” and “crow's feet Mendelian randomization” did **not identify a direct published proteome-wide MR/colocalization analysis using these two Korean GWAS outcomes**. This is a promising novelty signal, but a final systematic/preregistered literature check is required before claiming novelty in a manuscript.

---

## 2. Why generic facial-aging MR is becoming saturated

The most commonly reused outcome is UK Biobank Field 1757 / OpenGWAS `ukb-b-2148`, which captures whether participants are told they look younger, older, or about their age. It is therefore a **perceived facial-age phenotype**, not a direct quantitative wrinkle, pigmentation, elasticity, or barrier measurement.

Published MR studies already cover, among others:

- obesity / glycemic / lifestyle traits → facial ageing
- circulating metabolites → facial skin ageing
- micronutrients → facial ageing
- sex and growth hormones → facial ageing
- gut / skin microbiome → facial ageing
- telomere length → facial ageing
- lipid multi-omics → facial ageing
- transcriptome/proteome multi-omics → facial ageing

By 2026, the space advanced further: a SEMA7A study integrated eQTL, pQTL, MR and colocalization, while a September 2026 BMC Medicine study integrated blood eQTLs, mQTLs, 2,923 pQTLs, NMR metabolites and hematologic traits in >400k UK Biobank participants.

**Implication:** UKB perceived age remains useful as a *secondary cross-phenotype comparator*, but is weak as the sole primary outcome for a new thesis.

---

## 3. Objective East-Asian skin phenotype resources

| Resource / phenotype | Population | Approx. N | Full summary statistics | Candidate role |
|---|---:|---:|---|---|
| Under-eye wrinkles | Korean women | Discovery GWAS 11,079; paper total 17,019 | **OpenGWAS ebi-a-GCST90094903; EBI accession provenance** | **Primary outcome** |
| Crow's feet | Korean women | Discovery GWAS 11,079; paper total 17,019 | **OpenGWAS ebi-a-GCST90094904; EBI accession provenance** | **Primary outcome** |
| Facial skin color L* | East Asian | 48,433 | **GCST90320257** | Powered pigmentation comparator |
| Facial skin color a* | East Asian | 48,433 | **GCST90320258** | Powered pigmentation comparator |
| Facial skin color b* | East Asian | 48,433 | **GCST90320259** | Powered pigmentation comparator |
| X-chr L*/a*/b* | East-Asian females | 42,770 | GCST90320260–62 | Optional sex-specific extension |
| Asian skin brightness | KR/VN/TH/ID women | ~2,762 | GCST90502777 | Validation / specificity |
| Elasticity R2/R5/R7 | KR/VN/TH/ID women | ~2,762 | GCST90502786–88 | Validation / specificity |
| Moisture, 3 sites | KR/VN/TH/ID women | ~2,762 | GCST90502789–91 | Validation / specificity |
| pH, 3 sites | KR/VN/TH/ID women | ~2,762 | GCST90502792–94 | Validation / specificity |
| Pore size | KR/VN/TH/ID women | ~2,762 | GCST90502795 | Validation / specificity |
| Sebum, 3 sites | KR/VN/TH/ID women | ~2,762 | GCST90502796–98 | Validation / specificity |
| TEWL, 3 sites | KR/VN/TH/ID women | ~2,762 | GCST90502799–801 | Validation / specificity |
| Wrinkle | KR/VN/TH/ID women | ~2,762 | **GCST90502802** | Independent-ish Asian wrinkle validation |
| SCINEXA wrinkle/pigment traits | Han Chinese women | 1,534 | Full stats not readily public | Locus-level biological comparison |
| Korean facial pigmented spots | Korean women | 17,019 | GWAS published; full-stat access needs re-audit | High-interest future outcome |

### Key interpretation

The two Korean wrinkle GWAS are the best starting point because they combine **objective phenotype + reasonable discovery N + OpenGWAS-accessible complete GWAS data + Korean/East-Asian relevance**. The 2025 four-country Asian cohort gives unusually rich cosmetic phenotypes but is considerably smaller, so it should generally be treated as replication/triangulation rather than the discovery dataset.

---

## 4. Exposure resources: proteomics / molecular QTL

| Resource | Population / platform | Scale | Role |
|---|---|---:|---|
| UKB-PPP | Mostly European; Olink Explore 3072 | 54,219; 2,923 proteins | High-power pQTL discovery |
| deCODE | Icelandic; SomaScan | 35,559; 4,907 aptamers | High-power cross-platform pQTL discovery/replication |
| CKB Olink/SomaScan | Chinese; Olink + SomaScan | ~3,976 | **EAS ancestry-matched validation; platform concordance** |
| CKB SomaScan atlas 2026 | Chinese; SomaScan v4.1 | 3,965; 7,289 proteins | Broad EAS proteome; 1,092 proteins with cis-pQTL |
| BioBank Japan / Japanese Olink | East Asian; Olink Explore 3072 | smaller cohorts | Additional EAS pQTL replication where accessible |
| GTEx v8 skin eQTL | Multiple tissues | skin: sun-exposed lower leg + not-sun-exposed suprapubic | Tissue mechanism / colocalization |
| eQTLGen | Blood | large | Secondary gene-expression triangulation |

### Critical QC point

Olink and SomaScan do not measure every protein equivalently. In ~3,976 Chinese CKB participants, matched proteins had only modest cross-platform abundance correlation (median rho ≈ 0.29), and only a subset of cis-pQTL signals colocalized across platforms. Protein-altering variants can also create binding/epitope artifacts.

Therefore a high-confidence protein should ideally satisfy several of:

- cis-pQTL rather than trans-only evidence
- MR signal with adequate F-statistic
- outcome/pQTL colocalization
- consistency in Olink and SomaScan when both exist
- replication or direction consistency in EAS pQTL data
- no obvious assay-binding artifact / protein-altering-variant concern, or explicit sensitivity analysis
- mechanistic support in skin tissue / relevant cell type

---

## 5. Exposure × outcome gap matrix

| Research combination | Literature saturation | Data readiness | Novelty signal | Thesis suitability |
|---|---:|---:|---:|---:|
| EUR pQTL → UKB perceived facial ageing | High | Very high | Low | Low–medium |
| Metabolites → UKB facial ageing | High | High | Low | Low |
| Acne → proteome-wide MR | High | High | Low | Low |
| Rosacea → microbiome/metabolite/protein MR | Medium | High | Medium | Medium |
| AGA → metabolite / molecular MR | Medium | High | Medium | Medium |
| **pQTL → Korean under-eye wrinkles** | **Low / no direct study found in targeted search** | **Very high** | **High** | **Very high** |
| **pQTL → Korean crow's feet** | **Low / no direct study found in targeted search** | **Very high** | **High** | **Very high** |
| **Joint under-eye vs crow's-feet proteogenomics** | **Very low** | **Very high** | **Very high** | **Very high** |
| pQTL → EAS skin color L*/a*/b* | Low–medium | Very high | High | High, but less “aging”-specific |
| pQTL → Asian elasticity/moisture/TEWL/pore/sebum | Very low | Medium | Very high | Medium due to N/power |
| Objective Korean wrinkles ↔ UKB perceived-age cross-phenotype comparison | Low | High | High | High as secondary aim |

---

## 6. Candidate thesis/project designs

| Rank | Candidate design | Novelty | Executability | Main risk | Recommendation |
|---:|---|---:|---:|---|---|
| **1** | **Proteome-wide cis-pQTL MR + coloc of Korean under-eye & crow's-feet wrinkles** | Very high | Very high | EAS-vs-EUR pQTL ancestry mismatch if discovery uses European pQTLs | **Core thesis candidate** |
| **2** | **Site-specific causal architecture: common vs discordant proteins for under-eye vs crow's feet** | Very high | Very high | Outcomes are correlated / same cohort | **Integrate as Aim 2 of #1** |
| **3** | **East-Asian skin-trait atlas: wrinkles + skin color + elasticity/moisture/TEWL** | Very high | Medium | Small N for 2025 Asian traits; multiple testing | Strong extension |
| **4** | **Objective Korean wrinkles vs European perceived facial age** | High | High | Cross-ancestry and phenotype-definition heterogeneity | Strong secondary / validation aim |
| **5** | **Proteogenomic causal mapping of EAS objective pigmentation** | High | High | Skin color is partly constitutive pigmentation, not strictly aging | Good separate paper / backup |

**Preferred architecture:** combine #1 + #2 as the thesis core, then use #3/#4 selectively for external triangulation rather than making the project too broad.

---

## 7. Recommended thesis question

### Primary question

**Which genetically proxied circulating proteins influence objectively measured periorbital wrinkling in Korean women, and which signals are supported by shared causal variants and skin-tissue biology?**

### Secondary questions

1. Which protein effects are shared between under-eye wrinkles and crow's feet, and which are site-specific?
2. Do high-confidence protein signals replicate or show concordant direction in East-Asian pQTL resources?
3. Are the same signals observed for UKB perceived facial ageing, or are they specific to objective wrinkling?
4. Do the implicated genes/proteins localize to fibroblasts, keratinocytes, melanocytes, immune cells, vascular/endothelial cells, or other skin compartments?
5. Are any high-confidence targets druggable or linked to existing therapeutic/cosmetic mechanisms?

---

## 8. Proposed title

**Proteogenomic dissection of site-specific facial wrinkling in Korean women: Mendelian randomization, colocalization, and cross-ancestry validation**

Alternative conservative title:

**Genetically proxied plasma proteins associated with objective periorbital wrinkling in Korean women: Mendelian randomization and colocalization analyses**

---

## 9. Analysis architecture

### Stage 0 — Data audit and harmonization

**Outcome**
- Use OpenGWAS IDs `ebi-a-GCST90094903` and `ebi-a-GCST90094904` as the operational data source; retain GCST accession metadata for provenance.
- Stage 1 should query only the cis-pQTL instrument SNPs rather than downloading all ~4.94M outcome variants.
- Retrieve full VCF / regional outcome data only when required for locus-level colocalization or fine-mapping.
- Confirm genome build, variant IDs, effect/reference alleles, beta/SE/P, EAF, N, imputation/QC fields.
- Audit allele-frequency distribution and duplicated/multiallelic variants.
- Standardize to one build; use explicit liftover only when needed.

**Exposure**
- Build cis-pQTL catalog from UKB-PPP and/or deCODE.
- Prioritize ±500 kb to ±1 Mb cis window around protein-coding gene.
- Select genome-wide significant instruments.
- Use ancestry-matched LD for clumping / conditional independence.
- Record F-statistics; exclude weak instruments.
- Annotate PAV/missense/epitope-binding risk.

**EAS exposure validation**
- Query CKB Olink/SomaScan pQTLs and BBJ/Japanese resources where available.
- Record whether the same protein has an EAS cis-pQTL and whether the signal/direction is compatible.

### Stage 1 — Proteome-wide cis-pQTL MR

For each protein × wrinkle phenotype:

- single instrument: Wald ratio
- multiple independent cis instruments: IVW as primary
- report beta, SE, 95% CI, P, F, number of instruments
- harmonize strands; handle palindromic SNPs using EAF
- control multiple testing across tested proteins × outcomes (FDR, with Bonferroni as stringent reference)

Do **not** mechanically require MR-Egger/weighted median for proteins with too few cis instruments. Those methods are useful only when instrument count permits meaningful estimation.

### Stage 2 — Locus-level causal validation

For MR-positive or prioritized proteins:

- regional pQTL + wrinkle GWAS extraction
- ancestry-aware LD
- coloc.abf first-pass
- prior sensitivity analysis (including p12)
- PP.H4 threshold pre-specification, e.g. >0.8 as strong evidence
- SuSiE / conditional colocalization when multiple independent signals are present
- report credible sets and top posterior variants

For Korean outcomes, **do not use a European LD matrix for SuSiE merely because the pQTL discovery cohort is European**. Separate the ancestry issue explicitly and prioritize EAS LD / EAS pQTL replication for interpretation.

### Stage 3 — Site-specific wrinkle architecture

Classify proteins as:

- shared: significant/colocalized in both under-eye and crow's feet
- under-eye predominant
- crow's-feet predominant
- directionally discordant
- MR-only without colocalization

Compare effect sizes rather than merely comparing significance status. Because both traits derive from overlapping/same participants, correlated-outcome inference must account for covariance where formal effect-difference testing is attempted.

### Stage 4 — Skin-tissue mechanism

Integrate:

- GTEx v8 sun-exposed lower-leg skin eQTL
- GTEx v8 non-sun-exposed suprapubic skin eQTL
- eQTL colocalization / SMR where appropriate
- Human Protein Atlas tissue/cell-type expression
- public human skin scRNA-seq datasets
- fibroblast / keratinocyte / melanocyte / immune / endothelial annotation
- ECM, collagen, elastin, senescence, oxidative-stress, pigmentation and inflammatory pathway enrichment

### Stage 5 — Cross-phenotype / cross-ancestry triangulation

Priority order:

1. CKB / Japanese EAS pQTL validation
2. Korean under-eye vs crow's feet concordance
3. 2025 Asian wrinkle GCST90502802
4. EAS skin-color L*/a*/b* (GCST90320257–59)
5. UKB perceived facial ageing as a secondary comparator
6. elasticity/moisture/TEWL/pore/sebum only for sufficiently instrumented prioritized proteins

### Stage 6 — Translational annotation

For Tier-1 targets:

- Open Targets / drug-target evidence
- approved or investigational drugs
- target directionality
- skin expression / localization
- known collagen/ECM, melanogenesis, barrier, inflammation, senescence mechanism
- adverse phenotype PheWAS

Avoid equating a genetically supported circulating-protein effect with efficacy of a topical cosmetic ingredient. That requires separate pharmacologic and tissue-delivery evidence.

---

## 10. Candidate prioritization framework

### Tier 1

- cis-pQTL MR passes multiplicity threshold
- strong outcome/pQTL colocalization
- effect is directionally compatible in an EAS pQTL resource or across platforms
- relevant skin tissue/cell expression
- no major assay/PAV artifact concern

### Tier 2

- significant MR + moderate colocalization or strong skin mechanism
- not fully replicated across ancestry/platform

### Tier 3

- MR-only signal, trans-pQTL-only signal, or signal with poor colocalization
- exploratory until additional support is obtained

The thesis should be written around Tier-1 targets, with Tier-2 as supportive and Tier-3 in supplement/exploratory sections.

---

## 11. Main threats to validity

### Ancestry mismatch

European pQTL → Korean GWAS is analytically possible but can be biased or weakened by ancestry differences in LD, allele frequency and pQTL architecture. The solution is not to discard high-powered European pQTL resources, but to treat them as discovery and seek **EAS pQTL confirmation**.

### LD for colocalization / fine-mapping

LD must match the analyzed population as closely as possible. Use 1000G EAS as a minimum public reference; a Korean reference is preferable if available and legally accessible.

### Assay-binding artifacts

Affinity-proteomics cis-pQTLs can represent altered reagent binding rather than true concentration change, especially near protein-altering variants. Cross-platform Olink/SomaScan comparison and PAV annotation are required.

### Outcome covariate adjustment

Review exact wrinkle-GWAS model adjustment. Outcome GWAS covariate adjustment can affect causal interpretation if the exposure is upstream of an adjusted covariate.

### Multiple testing

A proteome-wide screen across two correlated outcomes requires explicit multiplicity control and a pre-specified primary family of hypotheses.

### Phenotype correlation

Under-eye and crow's-feet wrinkles are related traits. “Significant in one but not the other” is not evidence of a true effect difference. Effect-difference claims require formal modeling or covariance-aware inference.

---

## 12. Proposed project directory

```text
/srv/is-analysis/IS_Analysis_V3/
├── data/skin_beauty/
│   ├── outcomes/
│   │   ├── GCST90094903_under_eye/
│   │   ├── GCST90094904_crows_feet/
│   │   ├── GCST90502802_asian_wrinkle/
│   │   └── skin_color_eas/
│   ├── pqtl/
│   │   ├── ukb_ppp/
│   │   ├── decode/
│   │   ├── ckb/
│   │   └── bbj_japan/
│   ├── eqtl/
│   └── ld_reference/
│       └── eas/
├── results/skin_beauty/
│   ├── stage0_audit/
│   ├── stage1_cis_mr/
│   ├── stage2_coloc/
│   ├── stage2b_susie/
│   ├── stage3_site_specific/
│   ├── stage4_skin_tissue/
│   ├── stage5_replication/
│   └── stage6_targets/
├── scripts/
│   └── skin_beauty_*.py|R
└── docs/
    └── SKIN_BEAUTY_MR_MASTER.md
```

---

## 13. Immediate execution sequence

1. **Download/audit GCST90094903 and GCST90094904.**
2. Confirm genome build and harmonize variant schema.
3. Build an initial cis-pQTL exposure table from resources already used in the CKD pipeline.
4. Audit EAS overlap and allele frequencies against 1000G EAS / available CKB pQTL.
5. Run a minimal proteome-wide cis-MR prototype on both wrinkle outcomes.
6. Select top candidates by FDR + effect consistency.
7. Run ancestry-aware regional colocalization for those candidates.
8. Only after the above, expand to SuSiE, cross-platform pQTL replication and skin tissue/cell annotation.

This sequence deliberately tests whether the project has a real signal before investing in the more expensive multi-omics stages.

---

## 14. Go / no-go criteria after the prototype

**GO** if at least one of the following occurs:

- one or more proteins pass FDR and show strong coloc support;
- several biologically coherent proteins show moderate coloc plus EAS/platform support;
- distinct but reproducible protein architectures emerge for under-eye vs crow's feet.

**PIVOT / EXPAND** if proteome-wide signal is weak:

- add eQTL/TWAS/SMR as a second molecular layer;
- use the 48,433-person EAS skin-color GWAS as a powered complementary outcome;
- broaden to the Asian skin-trait panel for phenotype triangulation rather than treating a null wrinkle MR as project failure.

---

## 15. Key references / datasets

1. **Lee SG et al.** Identification of Genetic Loci Associated with Facial Wrinkles in a Large Korean Population. *Journal of Investigative Dermatology* (2022). GWAS Catalog: GCST90094903, GCST90094904.
2. **Kim B et al.** Mapping and annotating genomic loci to prioritize genes and implicate distinct polygenic adaptations for skin color. *Nature Communications* (2024). GWAS Catalog: GCST90320257–GCST90320262.
3. **Comprehensive Profiling of Genetic and Nongenetic Factors that Influence Skin Traits in Asian Women from 4 Countries.** *Journal of Investigative Dermatology* (2025). GWAS Catalog: GCST90502777–GCST90502802.
4. **Sun BB et al.** Plasma proteomic associations with genetics and health in the UK Biobank. *Nature* (2023). UKB-PPP: 54,219 participants; 2,923 proteins.
5. **Ferkingstad E et al.** Large-scale integration of the plasma proteome with genetics and disease. *Nature Genetics* (2021). deCODE SomaScan: 35,559 Icelanders; 4,907 aptamers.
6. **Wang B et al.** Comparative studies of 2168 plasma proteins measured by two affinity-based platforms in 4000 Chinese adults. *Nature Communications* (2025). CKB Olink/SomaScan cross-platform pQTL resource.
7. **Pozarickij A et al.** Genomic atlas of 7,000 plasma proteins and their associations with diseases and traits in East Asian populations. *medRxiv* (2026). CKB SomaScan: 7,289 proteins in 3,965 Chinese adults; 1,092 proteins with cis-pQTL.
8. **Multi-omics Mendelian randomization revealing SEMA7A as potential drug target for facial skin aging.** *Journal of Big Data* (2026).
9. **Causal inference maps the sex-dimorphic genetic architecture and inflammatory-metabolic networks underpinning perceived aging in UK Biobank.** *BMC Medicine* (2026-09-24).

---

## 16. Current recommendation

Proceed with **Korean objective facial wrinkle proteogenomics** as the skin/beauty MR prototype.

The distinguishing features are:

- objective cosmetic phenotype rather than perceived age;
- full public GWAS statistics;
- East-Asian/Korean population;
- two anatomically distinct wrinkle outcomes from the same research program;
- direct methodological reuse of the existing pQTL → MR → coloc → SuSiE → tissue/cell pipeline;
- strong opportunity to add EAS pQTL validation and cross-platform QC;
- close continuity with Hong-Hee Won's prior skin-genetics work.

The first computational milestone is **Stage 0 outcome download + schema/QC audit for GCST90094903 and GCST90094904**.


---

## 15. Stage 0 executable implementation (v0.2)

Generated companion scripts:

- `prepare_skin_stage0.py`: inventories existing UKB-PPP/deCODE cis-pQTL files and audits local outcome tables.
- `download_skin_wrinkle_opengwas.R`: authenticates to current OpenGWAS with `OPENGWAS_JWT`, records metadata/file links, and validates queryability through top hits.
- `run_skin_stage0.sh`: server launcher. It deliberately does **not** use `set -e`; component exit codes are printed explicitly.

Recommended server locations:

```text
/srv/is-analysis/IS_Analysis_V3/scripts/skin_beauty/
/srv/is-analysis/IS_Analysis_V3/data/skin_beauty/stage0_gwas/
/srv/is-analysis/IS_Analysis_V3/results/skin_beauty/stage0_audit/
```

OpenGWAS now requires JWT authentication for protected API endpoints. Tokens are currently valid for 14 days. Keep the token only in the environment (e.g. `OPENGWAS_JWT`); never commit it to Git or Drive.
