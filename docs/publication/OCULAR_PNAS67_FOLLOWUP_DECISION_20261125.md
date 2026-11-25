# Ocular separate disease-contact follow-up — decision, not panel refit

## Verified current evidence
All 67 exact Ensembl IDs from Singh et al. PNAS 2026 S2 were evaluated, without effect-based additions. Zero overlap with 241/16 does not establish novelty or non-regulation. Of 67, 48 tested in both original internal models, 37 both FDR<.05 and concordant, four with both absolute log2FC >=1. Current official 48h workbook reproduced 305 archived rows exactly but historical archive identity remains HOLD; 42/67 present in the current 48h table, 2/67 FDR<.05. Four internally strong rows:

| Gene | Original donor log2FC | Pool log2FC | Current 48h log2FC / FDR | Seven-day equal-strain mean log2 CPM ratio | All strain sign | Minimum arm mean CPM |
|---|---:|---:|---|---:|---:|---:|
| ANXA13 | +1.6639 | +1.5801 | absent | +0.1210 | 1/3 | 0.376 |
| TGFB3 | -1.1556 | -1.3011 | -0.02385 / .98959 | -1.5882 | 3/3 | 23.702 |
| LRIG1 | -1.0500 | -1.1392 | -0.27792 / .11950 | -0.4764 | 2/3 | 67.914 |
| TPPP3 | -1.0906 | -1.0700 | -0.56924 / .70285 | -0.3101 | 3/3 | 0.397 |

Prioritize TGFB3 for hypothesis generation because direction and effect magnitude hold across the three reported strains with detectable expression in both arms. LRIG1 has higher expression but strain heterogeneity; TPPP3 and ANXA13 have low minimum arm mean CPM. These ranks are post-readout triage only, not predeclared statistical tests. No proof of functional importance or clinical outcome. Neither 7-day nor 48-hour comparison is a new independent replication. No change to frozen 241/16.

## Next design
1. Keep all 67 publisher IDs as a separate fixed exploratory family. Retain the 48 known donor directions and the already observed 37/4 subgroups; do not impute unavailable directions. 
2. Inspect source coverage of all 67 on separately verified GSE37474, GSE65240 and GSE16643 platforms with a mapping/metadata gate before computing descriptive gene-level changes. Use exact Entrez/Ensembl mapping and source-specific denominators. No new P from same-study or technical replicate count; unresolved lot labels remain conditional. Report all four candidates together, including missing and discordant results. Avoid optimizing expression cutoffs.
3. For a genuinely new functional study, lock TGFB3 and a comparator or all four before readout; confirm mRNA/protein changes across independently sourced donor cultures, matched vehicle/DEX, defined dose/time, technical repeats averaged; randomize/blind and preregister analysis/power. Direct perturbation/rescue and functional outflow assays are separate experiments. No mechanistic claim based on association alone. No outreach or laboratory work performed.

Evidence: [current 48h](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45594940382/artifacts/15012727333), [7d context](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45595015809/artifacts/15012343987), corresponding JSON in evidence/. No main merge.
