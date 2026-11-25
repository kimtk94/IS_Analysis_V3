# Ocular frozen-panel ex-vivo tissue results (draft, 2026-11-25)

## Results text for manuscript

In GSE37474, we evaluated the frozen 241-gene glucocorticoid-response panel and its 16-gene ECM organization subset in five paired donor-eye samples exposed to 100 nM dexamethasone for 10 days versus untreated controls. Public sample titles and treatment, age and sex fields supported the pairing. Strict annotation and probe filters retained 27,920 eligible probe sets, with 199/241 (82.6%) and 14/16 (87.5%) candidate genes mapped. The median donor-level mean direction-aligned percentile shift was +5.21 percentage points for the 241-gene panel and +8.11 for ECM16. Both panels had positive scores in all five donor pairs (two-sided exact sign P = 0.0625 each; Holm-adjusted P = 0.125 each). The mean score remained positive in every leave-one-pair-out analysis (+4.93 to +5.78 and +8.02 to +9.45 percentage points, respectively). A post-readout sensitivity analysis collapsing eligible probes to 17,686 unique Entrez genes yielded median shifts of +5.49 and +8.27 percentage points, respectively; all five donor pairs remained positive for each panel. These results describe directional consistency in a mixed TM/corneoscleral tissue context, not statistically significant independent validation or clinical response prediction.

## Evidence and claim boundaries

| Analysis | Observed | Permitted interpretation |
|---|---|---|
| GSE313736 donor/pool pseudobulk | 57 paired discovery donors; panel also uses 48h published TM statistics | Discovery only; donor/pool models are not independent cohorts |
| GSE301517 seven-day | 234/241 equal-strain mean directions; 204/241 all three strains agree | Descriptive time-point corroboration; prior examined study, not blind holdout |
| GSE16643 cross-drug | 229/241 and 16/16 mapped; positive score in both strains for DEX, FA, TA | Descriptive drug-context consistency, not clinical drug equivalence |
| PRJNA729873 paired | 5/8 complete pairs; both panels 5/5 positive; Holm .125 | Incomplete data, missingness concentrated in nonresponders; not confirmatory |
| GSE37474 primary | 5/5 paired scores positive; Holm .125 each | Descriptive mixed-tissue consistency; donor independence not certified |
| GSE37474 sensitivity | Same five pairs; 5/5 positive each | Post-readout same-data sensitivity; no second independent test |

A common effect scale across studies is not established. Do not pool p values, patients, genes or ranks across studies to claim significance. Five pairs cannot achieve two-sided exact sign P < .05. No claim that an ECM pathway, a causal glaucoma mechanism, clinical IOP prediction, or an independent statistically significant replication has been newly established. Prior GSE37474 reanalyses and overlap with original selection must be checked; the manuscript's candidate contribution is the fixed-panel transportability assessment across platforms and conditions, not novel DEG discovery.

## Reproducible figure and evidence

- Primary: [run 45535697612](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45535697612); [aggregate JSON](evidence/ocular_exvivo_paired_20261125.json).
- Sensitivity: [run 45536989022](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45536989022); [aggregate JSON](evidence/ocular_exvivo_gene_rank_sensitivity_20261125.json).
- Figure 3: [run 45537145972](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45537145972), [artifact 14987506456](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45537145972/artifacts/14987506456). Actual job stdout reported SVG/PDF/PNG/legend and PNG dimensions 3260 x 1580. Upload-reported ZIP SHA256 3d9238eaf91255b055e2666ae67573c6468ea548b2f63f4e59f2c2dbd757c86e. File generation verified; independent ZIP digest and visual review pending.

## Before submission

1. Audit selection lineage, cohort relationships and exact prior fixed-signature work; use descriptive language throughout until independent evidence is established.
2. Seek legitimately public uniformly quantified expression for missing PRJNA729873 arms and PRJNA781014. Do not mix incompatible quantification pipelines or change the frozen panels. If unavailable, document this as a compute/data limitation.
3. Review actual Figure 2, Figure 3 and legends, units and denominators; preserve null/uncertain findings. No main merge or submission decision implied.
