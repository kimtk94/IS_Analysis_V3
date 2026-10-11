# IS — COL4A1/COL4A2 paired-cis colocalization diagnostic and prior sensitivity (2026-10-11)

**Project:** broad discovery of ischemic stroke (IS), BBJ Japanese ischemic stroke and GIGASTROKE EAS acute ischemic stroke (AIS); eQTL Catalogue complete cis inputs.
**Study status:** official input integrity and 270-model coloc.abf calculations **EXECUTED / TECHNICALLY VERIFIED**; human causal gene identification **EXPLORATORY / NOT ESTABLISHED**. No original data, no existing 646 coloc posterior, and no 2,225-gene candidate master have been changed.

## Executive conclusion

The previously discovered full-cis COL4A2–ROSMAP DLPFC (QTD000434) tag-SNP coincidence is real as an **allele-harmonized marginal-association record**. When the immediate neighboring gene **COL4A1** is given exactly the *same SNP panel*, however, a single-causal-signal coloc.abf gives nearly identical strong H4 posterior for both genes:

- **COL4A1:** PP.H3=0.022883, PP.H4=**0.948198** at p12=1×10⁻⁵.
- **COL4A2:** PP.H3=0.034158, PP.H4=**0.957119** at p12=1×10⁻⁵.

Both genes have the same top **H4-conditional tag SNP `13:110954510:T:C` (GRCh37)**, and the adjacent gene eQTL Z-score vectors have Pearson r=**0.823636** in ROSMAP DLPFC. *Neither high H4 nor this shared tag independently identifies the causal transcript.* The analysis prioritizes a **13q34 regulatory locus**, not an exclusive COL4A2 therapeutic or causal target.

The **molecular eGene study-wide permutation FDR** nonetheless differs: COL4A2 ROSMAP p_beta=0.00627568, full-study BH q=**0.0171972** (PASS); COL4A1 p_beta=0.078202, q=**0.156907** (NOT PASS). This provides a molecular-QTL evidence distinction but does not overcome shared LD-tagging or biological cell-type ambiguity.

**Why not claim a validated causal gene?** The 3,946 common non-palindromic SNPs are only **~45.0% of the 8,764 COL4A2 official cis-tested variants**. The remaining variants include source INDELs, allele-count quality exclusions, and GWAS missingness. The source ROSMAP molecular QTL donors are not the GIGASTROKE EAS AIS cohort; ancestry-matched LD and multiple causal-signal tests were not available/verified. The EAS AIS sample case fraction was cross-checked against secondary publications, not the original cohort-by-cohort source, and GIGASTROKE contains or may contain BBJ contributors. Both datasets are **not independent stroke replication**.

## 1. Source inputs and exact SNP overlap

We compared COL4A1 and COL4A2 in each *identical physical SNP set* rather than prioritizing one gene's own best-matched SNPs. The public original eQTL Catalogue `QTD000131`, `QTD000171`, `QTD000434` complete cis slices were indexed against the separately audited BBJ and GIGASTROKE AIS summary statistics, with alleles independently reconstructed against the GRCh37 reference FASTA (not canonical GIGASTROKE other:effect pseudo-ID).

| Study | Gene | Official cis genetic variants | Common BBJ/EAS AIS non-palindromic SNPs used | QTL N from AN/2 median | Same physical SNP: QTL p≤1e−4 and EAS AIS p<5e−8 |
|---|---|---:|---:|---:|---:|
| GTEx Artery Aorta | COL4A1 | 9,025 | 4,730 | 387 | 0 |
| GTEx Artery Aorta | COL4A2 | 9,026 | 4,730 | 387 | 0 |
| GTEx Brain Cortex | COL4A1 | 8,218 | 4,648 | 205 | 0 |
| GTEx Brain Cortex | COL4A2 | 8,218 | 4,648 | 205 | 0 |
| ROSMAP DLPFC | COL4A1 | 8,763 | 3,946 | 560 | **1** |
| ROSMAP DLPFC | COL4A2 | 8,764 | 3,946 | 560 | **2** |

Input TSV rows: 26,648 gene–tissue–SNP records (double counting the same physical SNP when measured for both genes). All six assays include beta, SE, p, ALT AF for QTL, BBJ and GIGASTROKE. No variants from alternative sample sets or unverified sources were imputed.

For ROSMAP (QTD000434), **1,191 genomic SNP records per gene** had discrepancies between nominal source `MAF` and `AC/AN` and were conservatively excluded before effect comparison. For each gene, SNP/INDEL categories and the tested `n_variants` remain separately documented. Never describe the 3,946 matches as the 8,764 full original SNP universe.

## 2. Coloc model and preregistered sensitivity variables

Package **coloc 5.2.3**, function `coloc.abf` and normal-approximation GWAS/QTL beta/SE. Quantitative QTL N is the median of original per-variant `AN/2`, 387/205/560; the quantitative trait sdY is internally inferred by coloc from source N, MAF, var(beta), **not observed from individual donor expression records**. Binary AIS data N=256,274, cases=19,032, case fraction=0.074264264 (secondary reported source arithmetic); original per-SNP N is not provided.

- Primary variant set: all matched non-palindromic variants in both genes and both stroke GWAS datasets, identical for the two genes per QTL study.
- Sensitivity variant sets: the same matched SNPs with QTL vs EAS AIS ALT allele-frequency differences ≤0.2 and ≤0.1 *in both genes*. These are QC stress tests, **not validated ancestry matching** and post hoc AF filters can change inferred posterior.
- Prior p1=p2=1e−4, p12 grid **1e−6, 3e−6, 1e−5, 3e−5, 1e−4**.
- Molecular QTL N perturbation: **0.8×, 1.0×, 1.2×** solely as a sensitivity analysis, without altering the published input file.
- Total **2 genes × 3 tissues × 3 SNP selections × 5 priors × 3 N settings = 270 successful coloc.abf fits**. All five posterior components H0–H4 are numerically finite and sum to 1 per fit.
- `coloc::check_dataset` warns for some eQTL datasets whose strongest P is weaker than its default threshold; this is expected when genes have no source-significant eQTL, and is **not ignored as evidence**. No data rows were invented to satisfy it.

**Critical method restriction:** coloc.abf assumes up to one causal signal per trait. It can produce a high H4 for neighboring genes that share strongly LD-correlated marginal tag signals. This is not the same as a multi-signal credible-set comparison; individual donor genotype-derived matched LD is not available for this locus. See: <https://chr1swallace.github.io/coloc/reference/coloc.abf.html> and <https://chr1swallace.github.io/coloc/articles/a02_data.html>.

## 3. Main result: exact SNP-matched tissue-by-gene posterior

| Public QTL study | COL4A1 PP.H3 | COL4A1 PP.H4 | COL4A2 PP.H3 | COL4A2 PP.H4 |
|---|---:|---:|---:|---:|
| GTEx Aorta | 0.324387 | 0.049896 | 0.260651 | 0.147165 |
| GTEx Cortex | 0.271206 | 0.100936 | 0.274500 | 0.148056 |
| ROSMAP DLPFC | **0.022883** | **0.948198** | **0.034158** | **0.957119** |

These rows use full available common non-palindromic SNP matches, original AN-derived N, and `p12=1×10⁻⁵`. The QTL eGene evidence is source-level FDR, not coloc H4. In all three tissues the H4-conditional top SNP is the same `13:110954510:T:C` for both genes (but **do not interpret the conditional SNP posterior as a confirmed causal variant**).

### ROSMAP H4 changes sharply with p12

| p12 | COL4A1 H4 | COL4A2 H4 |
|---|---:|---:|
| 1e−6 | 0.646697 | 0.690596 |
| 3e−6 | 0.845947 | 0.870063 |
| 1e−5 | 0.948198 | 0.957119 |
| 3e−5 | 0.982115 | 0.985286 |
| 1e−4 | 0.994566 | 0.995540 |

ROSMAP at `p12=1e-5`, source N and AF-band sensitivity: COL4A1 H4=0.9482 (full 3,946), 0.9537 (AF diff≤0.2, 3,095), 0.9565 (AF diff≤0.1, 2,053). COL4A2 H4=0.9571, 0.9577, 0.9588, respectively. QTL N 0.8× and 1.2× runs are preserved in the `270`-row source table. These stability checks do **not** cure the missing SNP universe, LD, cohort ancestry or single-causal assumptions.

### Gene-to-gene correlation and locus-gene resolution

| Molecular QTL source | COL4A1 versus COL4A2 marginal cis Z-score r |
|---|---:|
| GTEx Aorta | 0.807454 |
| GTEx Cortex | 0.868490 |
| ROSMAP DLPFC | 0.823636 |

These are *across-SNP marginal-effect* Pearson correlations, not gene-expression correlation across individuals. They reflect substantial shared locus structure/LD; no direct regulatory direction is inferred. In the selected ROSMAP SNP universe the strongest cis QTL SNP for both genes is `13:110954486:T:C`. The H4-conditional top shared tag is 24 bp downstream at `13:110954510:T:C`.

## 4. Relation to prior 646 internal ABF and published cerebrovascular studies

The existing **16 GTEx tissue assays per gene** in prior BBJ coloc reruns show p12=1e−5 maximum H4 of **0.1171 for COL4A1** and **0.3308 for COL4A2**, while the new paired ROSMAP-DLPFC–EAS-AIS analysis gives 0.9482 and 0.9571. These **are not direct replication statistics**: tissue, donors, case phenotype, sample ancestry, quantified genotype/QTL regions and source SNP inclusion differ. Do not claim that newly recovered data "corrected" the previous 646 posterior.

**Prior literature, explicitly separated from our current AIS analysis:**

1. **Rannikmäe et al., Neurology 2017**, *COL4A2 is associated with lacunar ischemic stroke and deep ICH* (DOI 10.1212/WNL.0000000000004560). In primarily **European ancestry lacunar stroke**, COL4A2 intronic **rs9515201** was associated with lacunar IS (OR 1.17; p=6.62×10⁻⁸); distinct phenotype from all-etiology EAS AIS. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5664302/>
2. **2019 Brain GWAS of cerebral small vessel disease**, rs9515201 at 13q34 relates to *non-lobar ICH*, with nominal eQTL signatures for COL4A1 and COL4A2 in putamen. This reinforces that two gene expressions in the locus may track a common regulatory background and that ICH is **not** AIS. <https://pmc.ncbi.nlm.nih.gov/articles/PMC6763741/>
3. **2025 Human brain vascular multi-omics**, a plausible COL4A2 regulatory element at rs9515201 is active in **smooth muscle, pericyte and perivascular fibroblast** populations with chromatin and peak-to-gene support. This is stronger cell-type functional context than whole-tissue DLPFC but **does not validate the present H4 as a COL4A2-specific causal effect**. <https://pmc.ncbi.nlm.nih.gov/articles/PMC12321221/>
4. **Collagen IV pathobiology:** rare COL4A1/COL4A2 variants can disrupt secretion and vascular basement membrane integrity, with strong mechanistic evidence predominantly in **hemorrhage/small-vessel disease**, not automatically broad ischemic stroke. <https://pmc.ncbi.nlm.nih.gov/articles/PMC3257894/>

**Important position-based comparison:** `rs9515201` maps to **GRCh37 chr13:111040798:A:C** in the same paired data, but in the available GIGASTROKE EAS AIS cohort its p=**0.03702**, not genome-wide significant. In ROSMAP DLPFC its COL4A1 QTL p=**0.636974** and COL4A2 QTL p=**0.945829**. This is **not** the same physical SNP as our two independently FASTA-harmonized joint association variants at `13:110954237:A:G` and `13:110954510:T:C`. Any LD tagging relationship between them remains **unknown** in an appropriate ancestry- and cohort-matched panel. Neither rs9515201 nor two other SNPs should be silently described as the same variant.

## 5. Status, blockers and decision

**VERIFIED (technical):**
- Official source QTL 19-column nominal schema / `permuted n_variants` and gene ID, original AN-derived N.
- BBJ and GIGASTROKE historical OA:EA to independent GRCh37 FASTA biological REF:ALT.
- Paired gene-by-study common non-palindromic SNP universe, full source and signed effects/allele frequency.
- 270-coloc.abf model runs with finite posterior, matched input SNP and 5 priors, 3 N and 3 AF filters.
- The association status and source-level BH q values, retained separately from posteriors.

**EXPLORATORY:** shared 13q34 AIS + molecular QTL locus; ROSMAP H4 signal for BOTH genes; cell-type hypothesis leaning to mural/perivascular compartments from external published data.

**NOT VERIFIED / BLOCKED:** cohort-matched EAS and ROSMAP DLPFC genotype LD, multi-signal SuSiE credible sets, independent non-BBJ EAS AIS cohort, original AIS case-control breakdown by contributing cohort, 55% missing original cis variants in GWAS-matched paired panel, transcript expression perturbation / patient tissue causal mechanism, **unique causal COL4A1 or COL4A2 designation**.

**Decision:** retain BOTH COL4A1 and COL4A2 within the original unrestricted **2,225 broad IS candidate universe**. Mark ROSMAP shared locus `EXPLORATORY_SHARED_REGULATORY_REGION`, not new unequivocal COL4A2 causal gene; do not replace previous 646 assay results. For highest return from locally available datasets, pursue small-vessel/vascular cell-type annotations and cross-ancestry LD proxy sensitivity separately, without claiming matched LD when absent.

## 6. Reproducibility and figures

Python inputs:
- `scripts/is/prepare_is_col4a1_col4a2_paired_cis_coloc.py`
- `scripts/is/summarize_is_col4a1_col4a2_paired_abf.py`

R analysis:
- `scripts/is/run_is_col4a1_col4a2_paired_abf_sensitivity.R`
- `scripts/is/render_is_col4a1_col4a2_abf_prior_grid.R`
- `scripts/is/render_is_col4a1_col4a2_locus_associations.R`

Original derivative data (only on server):
- `/srv/is-analysis/results/is/audits/is_col4a1_col4a2_paired_coloc_20261011_v1/` — 26,648 matched SNP rows with exact original QTL / BBJ / EAS AIS signed association statistics, 6 input-QC rows, SHA256.
- `/srv/is-analysis/results/is/audits/is_col4a1_col4a2_paired_abf_sensitivity_20261011_v1/` — all 270 posterior sensitivity grid rows, QC and two figures.
- `/srv/is-analysis/results/is/audits/is_col4a1_col4a2_paired_abf_science_gate_20261011_v1/` — 6 gene×study summary, 30 prior rows, 54 AF-N QC rows and source hashes.

GitHub small derivative tables and Figures:
- [Single causal posterior grid: paired H3 and H4](../figures/is/IS_COL4A1_COL4A2_THREE_TISSUE_H3_H4_PRIOR_SENSITIVITY.svg)
- [ROSMAP DLPFC QTL and GIGASTROKE EAS AIS locus](../figures/is/IS_COL4A1_COL4A2_ROSMAP_EAS_AIS_GWAS_CIS_LOCUS_PLOT.svg)
- `docs/data/is/IS_COL4A1_COL4A2_6_TISSUE_PAIRED_ABF_EXPLORATORY_CONCLUSIONS.tsv`
- `docs/data/is/IS_COL4A1_COL4A2_THREE_QTL_P12_AF_N_ABF_SENSITIVITY.tsv`
- `docs/data/is/IS_COL4A1_COL4A2_PAIRED_ABF_SCIENCE_GATE_MANIFEST.json`

Heavy source QTL and full genotype variant matrices remain on the server only. This is a source-anchored **diagnostic analysis**, not a finished colocalization gene-discovery manuscript result.
