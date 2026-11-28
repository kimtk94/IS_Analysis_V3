# Ocular frozen-panel ex-vivo tissue results (draft, 2026-11-25)

## Results text for manuscript

In GSE37474, we evaluated the frozen 241-gene glucocorticoid-response panel and its 16-gene ECM organization subset in five paired donor-eye samples exposed to 100 nM dexamethasone for 10 days versus untreated controls. Public sample titles and treatment, age and sex fields supported the pairing. Strict annotation and probe filters retained 27,920 eligible probe sets, with 199/241 (82.6%) and 14/16 (87.5%) candidate genes mapped. The median donor-level mean direction-aligned percentile shift was +5.21 percentage points for the 241-gene panel and +8.11 for ECM16. Both panels had positive scores in all five donor pairs (two-sided exact sign P = 0.0625 each; Holm-adjusted P = 0.125 each). The mean score remained positive in every leave-one-pair-out analysis (+4.93 to +5.78 and +8.02 to +9.45 percentage points, respectively). A post-readout sensitivity analysis collapsing eligible probes to 17,686 unique Entrez genes yielded median shifts of +5.49 and +8.27 percentage points, respectively; all five donor pairs remained positive for each panel. These results describe directional consistency in a mixed TM/corneoscleral tissue context, not statistically significant independent validation or clinical response prediction.

## Evidence and claim boundaries

| Analysis | Observed | Permitted interpretation |
|---|---|---|
| GSE313736 donor/pool pseudobulk | 57 paired discovery donors; panel also uses 48h published TM statistics | Discovery only; donor/pool models are not independent cohorts |
| GSE301517 seven-day | 241/241 mapped; 234/241 equal-strain mean directions; 204/241 all three strains agree | Descriptive time-point corroboration; prior examined study, not blind holdout |
| GSE16643 cross-drug | 229/241 and 16/16 mapped; positive score in both strains for DEX, FA, TA | Descriptive drug-context consistency, not clinical drug equivalence |
| PRJNA729873 paired | 5/8 complete pairs; both panels 5/5 positive; Holm .125 | Incomplete data, missingness concentrated in nonresponders; not confirmatory |
| GSE37474 primary | 5/5 paired scores positive; Holm .125 each | Descriptive mixed-tissue consistency; donor independence not certified |
| GSE37474 sensitivity | Same five pairs; 5/5 positive each | Post-readout same-data sensitivity; no second independent test |
| GSE65240 fourteen-day | 220/241 and 16/16 mapped; 3/3 conditional lot scores positive for both; lot 4973 alone positive | Metadata-conflict sensitivity only; source_name unresolved, biological independence unverified |

A common effect scale across studies is not established. Do not pool p values, patients, genes or ranks across studies to claim significance. Five pairs cannot achieve two-sided exact sign P < .05. No claim that an ECM pathway, a causal glaucoma mechanism, clinical IOP prediction, or an independent statistically significant replication has been newly established. Prior GSE37474 reanalyses and overlap with original selection must be checked; the manuscript's candidate contribution is the fixed-panel transportability assessment across platforms and conditions, not novel DEG discovery.

## Reproducible figure and evidence

- Primary: [run 45535697612](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45535697612); [aggregate JSON](evidence/ocular_exvivo_paired_20261125.json).
- Sensitivity: [run 45536989022](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45536989022); [aggregate JSON](evidence/ocular_exvivo_gene_rank_sensitivity_20261125.json).
- Figure 3: [run 45537145972](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45537145972), [artifact 14987506456](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45537145972/artifacts/14987506456). Actual job stdout reported SVG/PDF/PNG/legend and PNG dimensions 3260 x 1580. Upload-reported ZIP SHA256 3d9238eaf91255b055e2666ae67573c6468ea548b2f63f4e59f2c2dbd757c86e. File generation verified; independent ZIP digest and visual review pending.

## Before submission

1. Audit selection lineage, cohort relationships and exact prior fixed-signature work; use descriptive language throughout until independent evidence is established.
2. Seek legitimately public uniformly quantified expression for missing PRJNA729873 arms and PRJNA781014. Do not mix incompatible quantification pipelines or change the frozen panels. If unavailable, document this as a compute/data limitation.
3. Review actual Figure 2, Figure 3 and legends, units and denominators; preserve null/uncertain findings. No main merge or submission decision implied.

## Fourteen-day conditional results (GSE65240)

The frozen 241/16 panels were evaluated in three reported commercial primary human TM cell lots exposed to 100 nM DEX for 14 days. The original title, explicit lot/treatment characteristics and supplementary filenames supported the six intended identities (6/6); the original source_name field conflicts for four samples (lots 2584 and 3423). We retain the primary identity HOLD and report a conditional analysis of the depositor labels, plus the conflict-free lot 4973 alone. The two Agilent platforms were analyzed separately with 29,831 eligible probes per platform and exact common Entrez panel mappings (220/241, 16/16). Direction-aligned percentile shifts for the 241-gene panel were +7.86, +7.57 and +11.11 percentage points in lots 2584, 3423 and 4973, respectively; for ECM16, +9.92, +9.43 and +11.65. Median conditional shifts were +7.86 and +9.92 percentage points. The three lots cannot be assumed to be three independent donors. No inferential P is assigned. Concordant conditional results in the metadata-conflict-free lot are not independent validation, and no outcome-based label change was made.

## Cross-study evidence, not a time course

| Context | Frozen mapping (all; ECM) | Directional evidence (all; ECM) | Effect measure and boundary |
|---|---|---|---|
| 48h published TM context used in panel construction | selection inputs | Not an independent test | No new comparator estimate |
| GSE301517 seven-day, three strains | 241/241; 16/16 | 234/241; 16/16 equal-strain mean direction, and 204/241; 13/16 in every strain | Normalized CPM log2 ratios; technical repeats averaged, not independent donor confirmation |
| GSE37474 ten-day, five pairs | 199/241; 14/16 | 5/5; 5/5 positive | Median signed percentile shift +5.21; +8.11 pp, Holm .125 for both |
| GSE65240 fourteen-day, three commercial lots | 220/241; 16/16 | 3/3; 3/3 positive conditional on depositor labels | Median signed percentile shift +7.86; +9.92 pp; identity HOLD and donor independence not certified |
| GSE16643 cross-drug, two strains | 229/241; 16/16 | 6/6 strain-drug contrasts positive for both | Descriptive percentile shift; distinct drugs and repeats not six donors |
| PRJNA729873 incomplete paired study | 5/8 complete pairs | 5/5; 5/5 positive | Median signed log2 difference +0.8158; +0.9163, Holm .125; three missing nonresponders |

These are different studies with distinct tissues/cultures, exposure conditions, rank universes and mapping. The larger 14-day median cannot establish a temporal increase, stronger ECM effect or dose response. Gene-level directional counts cannot be pooled with donor-level score counts. The common panel lineage and prior dataset use remain limitations. The primary claim is descriptive portability of a fixed response signature across contexts, NOT a newly discovered glaucoma causal mechanism or clinical predictor. The 48-hour comparison is part of selection. Seven-day mapping counts are 241/241, not 234/241; the latter is the number of concordant equal-strain mean effects.

Verified [14-day score](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45563105992/artifacts/14998946580), aggregate [JSON](evidence/ocular_14d_conditional_score_20261125.json). Prior readbacks are archived as `ocular_7d_frozen_directional_readback_20261125.json`, `ocular_crossdrug_directional_readback_20261125.json`, and `ocular_incomplete_paired_readback_20261125.json` in evidence/. No main merge.
