# IS G0022 / ALDH2 rs671 — Stroke subtype comparison and mediation feasibility
**Analysis date:** 2026-10-10 KST
**Research branch:** `research/is-broad-discovery-20261009`
**Verdict:** **EAS/Japanese subtype association and molecular-QTL evidence established; ALDH2→stroke direct effect or mediated effect NOT established.**

## 1. Source and allele contract

GIGASTROKE EAS five phenotypes `GCST90104544` (Any Stroke, AS), `GCST90104545` (Any Ischemic Stroke, AIS), `GCST90104546` (Cardioembolic Stroke, CES), `GCST90104547` (Large Artery Stroke, LAS), `GCST90104548` (Small Vessel Stroke, SVS), plus BioBank Japan ischemic stroke GWAS. Sources taken from **existing canonical GRCh37 reference-harmonized original GWAS data**, no regenerated effect estimates.

All associations were oriented to **rs671 GRCh37 chr12:112241766 G>A, alternative allele A**, checking original effect allele, sign of beta, inverted EAF for reference-effect rows, 95% CI from original SE, and actual original p. Confirmed **24/24** of four G0022 AIS credible-set SNPs × six phenotype analyses had ALT protective direction. This is NOT 24 independent signals: the four SNPs are strongly correlated (`rs671 EAS r² ≈0.87–0.99`) and studies/phenotypes overlap.

## 2. Original rs671 stroke association by phenotype

| Source and phenotype | ALT A beta (log OR) | SE | OR per A (95% CI) | Original p |
|---|---:|---:|---|---:|
| BBJ Japanese ischemic stroke | −0.10975 | 0.01253 | **0.896 (0.874–0.918)** | 1.99e-18 |
| GIGASTROKE EAS AS | −0.09830 | 0.01400 | **0.906 (0.882–0.932)** | 2.226e-12 |
| GIGASTROKE EAS AIS | −0.15140 | 0.01760 | **0.860 (0.830–0.890)** | 8.426e-18 |
| GIGASTROKE EAS SVS | −0.15590 | 0.03120 | **0.856 (0.805–0.910)** | 5.87e-7 |
| GIGASTROKE EAS CES | −0.23400 | 0.07400 | **0.791 (0.685–0.915)** | 0.00156 |
| GIGASTROKE EAS LAS | −0.09250 | 0.05310 | **0.912 (0.822–1.012)** | 0.08141 |

**Interpretation**:
- Strong AIS/SVS associations and nominal CES association observed; LAS CI crosses 1 and p>0.05, **not evidence for true difference between subtypes**. CES has larger effect point estimate but much larger uncertainty (SE=0.074).
- These GIGASTROKE phenotypes are nested and share or may share cases/controls; BBJ also has possible cohort reuse within meta-analyses. They are **not six independent confirmations**, must not run naïve fixed-effect pooling, Cochran Q or formal subtype-difference p-values without covariance or demonstrably disjoint cohorts.
- Forest figure is descriptive association only, no adjusted risk, no population-attributable inference and no clinical practice claim.

### Generated research figure

`G0022_RS671_STROKE_PHENOTYPE_FOREST.png` shows ALT A odds ratios and 95% CIs, with figure caption warning about shared controls/cohorts. The source R code regenerates this file from audited TSV. It is not a formal meta-analysis plot.

## 3. Coding variation, East Asian QTL, and ALL alternative nearby genes

Top AIS full-locus CS variant ALDH2 rs671 is **missense coding functional** (Glu/Lys amino-acid change, greatly reduced aldehyde dehydrogenase activity). This is an **independent functional annotation**, not evidence of direct protective protein mechanism for stroke.

JCTF Japan Omics Browser in Japanese COVID-19 patient samples:
- rs671 → ALDH2 eQTL nominal p=3.51e-12, QTL SuSiE PIP=0.354.
- rs671 → ALDH2 Olink pQTL nominal p=8.42e-49, QTL SuSiE PIP=0.0926.
- All four GWAS CS variants had ALDH2 eQTL and Olink pQTL associations.
- The complete four-page browser evidence ledger contains **72 QTL entries for 19 distinct genes**, including **17 genes with eQTL** and **3 with pQTL**, representing **20 distinct gene–assay pairs** (ALDH2 appears in both).
- Other reported signals include BRAP pQTL (maximum JCTF PIP=0.0767), PCNPP1 eQTL (0.0735), RPH3A eQTL (nominal p down to 1.61e-21 but fine-mapping PIP=0 among four GWAS variants), HECTD4, NAA25, PTPN11 and more.

A tiny marginal QTL p-value and a positive disease GWAS at the same correlated locus **do not prove colocalization**; GWAS and QTL SuSiE PIPs are computed on different SNP sets and cannot be multiplied for an informal posterior. The JCTF public ZIP is thresholded (nominal p<0.05 or PIP>0.001), so missing non-significant variants blocks valid full-cis coloc. All **causal-gene and mediation gates remain closed**.

## 4. rs671 stroke causal DAG (testable, not confirmed)

```text
Genotype: ALDH2 rs671
   ├── ALDH2 catalytic activity ──> aldehyde/oxidative stress ──> vascular injury ──> stroke
   ├── alcohol flush/tolerance ──> alcohol behavior ──> BP / metabolic effects ──> stroke
   ├── ALDH2 mRNA abundance ──> ALDH2 enzyme and other cell-state effects ──> stroke
   └── missense protein epitope/conformation ──> Olink binding signal (may not equal abundance)
Sex, ancestry, culture and COVID patient ascertainment modify observed effects/assays.
```

Key human epidemiology context:
- Large Han Chinese targeted sequencing study reported **genotype-by-drinking interaction** for ALDH2 rs671 and ischemic stroke. It provides a concrete warning against declaring the locus a drug target solely from cis-pQTL beta.
- Prospective Korean ALDH2-tagging variant study examined sex-specific stroke risk and alcohol consumption, illustrating why sex and environmental context can modify interpretation.
- 2026 Taiwanese first-ever ischemic stroke cohort evaluated genotype/alcohol interaction in stroke onset age, a *different outcome* from incidence GWAS.

These external studies are used as biological/methodological context, **not pooled with GIGASTROKE**.

## 5. Explicit Mendelian randomization / mediation gating

The server currently has:
- Signed GWAS effects/SE for six (potentially overlapping) stroke phenotype datasets at four CS SNPs.
- Japanese variant-level eQTL/pQTL browser associations and the source study metadata.
- Signed 1000G EAS LD reference n=504 and exploratory rs671-conditional GWAS residual z analysis (no remaining GWS under one-lead approximation).

It **does not have verified sex-specific EAS stroke effect summaries, harmonized EAS alcohol consumption exposure GWAS, EAS blood-pressure exposure GWAS, complete tested JCTF cis-QTL summary, independent instruments or disjoint-sample validation.**

Therefore:
- **One rs671 genetic instrument cannot support MR-Egger intercept, Cochran Q or weighted median/robust pleiotropy diagnostics.** Its direct catalytic, behavior, protein and possible binding mechanisms violate a simple single-mediator exclusion-restriction assumption.
- Do not compute a causal Wald ratio `beta_AIS / beta_JCTF_ALDH2_QTL` and call it mediated ALDH2 effect. Even if ALT directions seem aligned, it would misattribute pleiotropic coding and alcohol pathways.
- A sound design requires EAS/sex-aware alcohol and SBP/DBP GWAS sources, instrument harmonization, mutually distinguishable loci or multivariable/mediation identification, formal QTL multi-signal colocalization with complete variant coverage, and independent EAS stroke data.
- Do not infer an absence of alcohol-mediated stroke risk from low or absent eQTL in non-East-Asian panels.

Eight separate pathways and their missing causal inputs are in the machine-readable `G0022_RS671_MEDIATION_PATHWAY_READINESS.tsv`; **0** valid mediation effects were estimated.

## 6. Reproducibility / code artifacts

Under `/srv/is-analysis/worktrees/is-broad-discovery-20261009`:
```bash
python3 scripts/is/audit_g0022_rs671_stroke_subtypes.py
python3 scripts/is/audit_g0022_jctf_all_gene_qtl.py
python3 scripts/is/audit_g0022_rs671_mediation_readiness.py
Rscript scripts/is/plot_g0022_rs671_stroke_subtypes.R
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Outputs under `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/`:
- `G0022_FOUR_CS_SIX_GWAS_EFFECTS.tsv`
- `G0022_RS671_SIX_GWAS_SUBTYPE_EFFECTS.tsv`
- `G0022_RS671_STROKE_SUBTYPE_SUMMARY.json`
- `G0022_RS671_STROKE_PHENOTYPE_FOREST.png`
- `G0022_JCTF_ALL_GENE_QTL_SIGNAL_AUDIT.tsv`
- `G0022_JCTF_ALL_GENE_QTL_SUMMARY.json`
- `G0022_RS671_MEDIATION_PATHWAY_READINESS.tsv`
- `G0022_RS671_MEDIATION_READINESS_SUMMARY.json`

All full IS gene candidates (2,225 unique positional Ensembl genes) are preserved; G0022 is a local mechanistic case study within broad IS discovery.

## 7. Next true identification priorities

1. **Sex-stratified EAS GWAS** for ischemic stroke, alcohol intake and SBP/DBP (allele/variant/case definition per exposure) to test proposed genotype–sex–alcohol links; not currently available in production data directory.
2. Source verified, full *unfiltered* cis-QTL variant universe around rs671 from JCTF or MAEEA (currently Zenodo restricted), allowing joint LD fine mapping and `coloc.susie`.
3. Independent EAS AIS replication (exclude BBJ reuse) and larger ancestry-matched LD panel with covariance-aware subtype analysis.
4. Orthogonal ALDH2 protein abundance measurements (LC-MS) to test Olink assay-binding effects caused by protein coding variation.
5. Vascular/brain tissue cell-state evidence for all gene alternatives. Avoid selecting ALDH2 solely due large pQTL magnitude.

## Primary study sources
- GIGASTROKE: https://www.nature.com/articles/s41586-022-05165-3
- Han Chinese ALDH2 rs671 interaction in stroke: https://pmc.ncbi.nlm.nih.gov/articles/PMC9673712/
- Korean sex-specific ALDH2 and stroke: https://pmc.ncbi.nlm.nih.gov/articles/PMC4317484/
- Taiwan onset and drinking interaction (2026): https://pmc.ncbi.nlm.nih.gov/articles/PMC12828095/
- rs671 ALDH2 biochemical function: https://www.nature.com/articles/s41467-024-46899-0
- JCTF rs671 variant QTL: https://japan-omics.jp/variant?input_value=rs671
- Japan Omics Browser help: https://www.japan-omics.jp/help
- NBDC JCTF summary metadata: https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000193
