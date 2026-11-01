# Ocular P6d/P6e — independent 100nM TM pathway comparison and shotgun raw route recovery (2026-10-31)

## Real verified GitHub Actions runs
- **P6d** [run #41173246350](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173246350) **SUCCESS**, [aggregate artifact #13031791978](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173246350/artifacts/13031791978).
- **P6e** [run #41173368834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173368834) **SUCCESS**, [aggregate artifact #13031452582](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41173368834/artifacts/13031452582).
- Artifacts opened and numerical content checked; no full FASTQs, donor metadata or sample-level patient information uploaded.
- Repeatable code: `scripts/ocular/ocular_p6d_crossstudy_pathway_concordance.py`, `scripts/ocular/ocular_p6e_shotgun_alternative_access.py`. Workflows are **manual dispatch only after first verified run**.

## P6d — identical gene background, opposite dose/platform studies

Compared original P5b **16,050** human TM BeamA edgeR genes (nominal **500nM DEX 48h**) with Mehrotra et al. IOVS 2026 official 48h **100nM bulk TM** supplementary DESeq2 contrasts, loaded live (original ZIP **23,188,342 bytes**). The external 48h table contains **13,313** Ensembl-gene rows, **12,251** of which map to primary Ensembl IDs, and **12,250** unique gene symbols after deterministic duplicate collapse. Official MSigDB `2025.1.Hs` Hallmark GMT SHA256 is held fixed (`f22066af72e215ccb7b89d88e492c07e1eef17534c2ca7b0f9902cfecbbdd8e9`).

Unlike P6b's P-value-weighted exploratory rank tests, P6d **uses the absolute log2FC magnitude to weight each ranked list**, so the numbers are **different estimands**. Both cohorts are evaluated using *the same* mapped gene-symbol background, avoiding a misleading platform-specific gene-universe difference. These `ES` values do **not** have cross-study P-values, are not a pooled effect, and do not test dose×study interactions.

| MSigDB Hallmark | 500nM BeamA ES | 100nM bulk TM ES | Same ES direction | Same log2FC sign among matched pathway genes |
|---|---:|---:|---|---|
| EPITHELIAL_MESENCHYMAL_TRANSITION | **-0.478209** | **-0.480715** | YES | **141/194 (72.68%)** |
| INFLAMMATORY_RESPONSE | -0.362696 | -0.386645 | YES | 70/110 (63.64%) |
| REACTIVE_OXYGEN_SPECIES_PATHWAY | +0.333068 | +0.505666 | YES | 27/44 (61.36%) |
| IL6_JAK_STAT3_SIGNALING | -0.351400 | -0.269417 | YES | 40/53 (75.47%) |
| APOPTOSIS | -0.315524 | -0.320401 | YES | 84/135 (62.22%) |
| UNFOLDED_PROTEIN_RESPONSE | **-0.343666** | **+0.287631** | **NO** | 55/104 (52.88%) |
| TGF_BETA_SIGNALING | +0.248864 | -0.322218 | NO | 30/52 (57.69%) |
| TNFA_SIGNALING_VIA_NFKB | -0.233868 | +0.307208 | NO | 114/164 (69.51%) |
| HYPOXIA | -0.300960 | +0.337237 | NO | 114/166 (68.67%) |

**5/9** examined Hallmark pathway ES signs agree in these two different culture experiments. The **EMT pathway is the most compelling follow-up among P6b's two exploratory FDR-positive pathways**, because its negative ranked ES is consistent and **49** EMT genes reach each cohort's gene-level BH-FDR<.05, of which **3** have opposite signed log2FC; however both studies contain substantial covariate and protocol differences, and these are not formal independent pathway replication p-values.

**UPR is specifically discordant (−/+), so P6b UPR FDR=.009 must not be taken as external validation or dose-specific mechanistic proof.** TGF-β and NF-κB also disagree at pathway ranking level. All **12** previously selected opposite-sign genes remain opposite-sign by construction: a selected 12-gene list cannot be counted as an independent held-out confirmation.

**GO** numerical cross-study descriptive directional ranking on a shared feature universe. **HOLD** true pharmacological specificity, protein-level ECM change, dose reversal, individual cell-state confounding, and clinical steroid-induced IOP.

## P6e — no missing raw archive ROUTES for 31 shotgun runs previously lacking FASTQ metadata

For the **62 WGS runs** associated with `PRJEB74248`, a fresh official ENA `read_run` report was retrieved including `fastq_ftp`, `submitted_ftp`, `sra_ftp`. Actual counts:

| Raw file metadata route | Number of runs |
|---|---:|
| `fastq_ftp` | **31/62** |
| `submitted_ftp` | **62/62** |
| `sra_ftp` | **0/62** |
| `fastq_ftp` AND `submitted_ftp` | **31/62** |
| `submitted_ftp` only | **31/62** |
| ALL routes absent | **0/62** |

This **corrects the preliminary P6c concern**. All 62 study runs have **a reported raw submission route**; the 31 runs without `fastq_ftp` are **not automatically missing data**. But `submitted_ftp` paths may lead to different raw-file formats; this pass **did not test any file downloads, file readability, host DNA depletion, sample metadata, or tear proteomic sample join**. Therefore the appropriate decision is:

- **GO** metadata-discoverable file routes **62/62**.
- **HOLD** confirmation of actual binary access, format, paired layout, sample case/control annotation, low biomass negative controls, matched tear proteomics, and clinical drug/IOP phenotype.
- **STOP** causal microbial effects on steroid-induced IOP based on this cross-sectional ocular study alone.

## Combined next focus
1. If expanding to **cell-specific DEX pathophysiology**, prioritize experimentally controlled ECM/stiffness/aqueous outflow validation and gene-level agreement under randomized 0/100/500nM donor-matched conditions, including MYOC as a glucocorticoid TM assay marker.
2. Treat 16S/shotgun findings as a **separate glaucoma ocular-surface and oral microbiome project**, until joint human patient sample IDs, negative controls, clinical eye-drop exposure and a steroid response outcome can be independently verified.
3. Any future shotgun pipeline requires downloaded raw file type identification, archive checksums, removal of host reads, negative-control decontamination and clear BioProject sample-level covariates. This PR contains no read processing for these microbiome runs.
