# Draft results and decision — ischemic stroke fixed 24 subtype extension

## Results

Of 192 marker-study cells from 24 historically selected EAS regional markers and eight GIGASTROKE ancestry/phenotype summaries, 170 were matched and 158/192 fixed contrasts were evaluable. Eleven contrasts passed Bonferroni correction across 192 under a zero-covariance assumption and four passed a maximum-variance bound. The latter included the known PITX2-region marker rs2129981 in EUR and EAS SVS-versus-CES contrasts (delta -0.2983 and -0.3118; adjusted P 3.49e-6 and 8.43e-4), rs2767060 in EAS-versus-EUR LAS (delta +0.3672, adjusted P 0.03283), and rs7203740 in EAS SVS-versus-CES (delta -0.2945, adjusted P 0.04376). These are conditional comparisons of reported coefficients from the same release, not direct case-case or causal estimates. Retrospective EAS marker selection, phenotype ascertainment, population LD and sample imbalance, model-scale comparability and independent replication remain unresolved.

At rs2767060, the exact C/T T-allele association in EAS LAS (beta 0.3582, SE 0.0680, P 1.38e-7) did not meet conventional genome-wide significance, while the EUR LAS association was near zero (beta -0.0090, SE 0.0297, P 0.7604). The 1000 Genomes Phase 3 reference comparison returned 78 unique high-LD proxies across CHB/JPT/CEU, 15 shared. No exact index or proxy token occurred in the pinned GIGASTROKE supplementary workbook or in the 51 Catalog regional variants retrieved from the 2026-11-17 release. The variant remains intergenic; no causal gene is established. Regional and independent literature coverage remains incomplete.

At rs7203740, the exact C/T (not C/A) T allele showed nominal PLCG2 association in ROSMAP DLPFC (beta 0.279006, SE 0.068636, P 5.51034e-5); the other four tested datasets did not have nominal P < 0.05 for PLCG2. The gene-level permutation BH P was 4.5327492e-6 among 23,076 tested genes. Its official lead rs11640696 T/G at chr16:81878856 differs from the GWAS candidate. Nominal and permutation files agree on lead variant identity and direction but not the exact estimates (nominal beta -0.249511 and P 1.04157e-10; permutation beta -0.246995 and nominal P 7.86674e-11). The mismatch remains a model/source-equivalence HOLD. Full hash-verified ROSMAP SuSiE output (377,180 rows, 7,083 gene IDs) contained no PLCG2 gene row or rs7203740 variant; no PIP can be assigned. The Catalog region contained 342 curated variants but no exact rs7203740 match. This is exploratory molecular context, not evidence of tissue specificity, colocalization, or mediation.

## Publication decision

- Primary framing: genetic effect heterogeneity across ancestry and ischemic stroke subtypes, with known positive controls, full fixed-family denominators, explicit missingness and conservative covariance bounds.
- rs2767060: exploratory cross-ancestry LAS candidate; require independent subtype-level evidence and molecular anchoring.
- rs7203740/PLCG2: exploratory subtype-associated molecular candidate; nominal-only regional comparisons with independently verified nominal anchor estimates may proceed, while cross-file estimate equivalence remains HOLD.
- No formal MR/colocalization, matched population LD, causal gene, independent validation or claim of novelty. No main merge.

## Verified evidence
- [Fixed family](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45534891070)
- [rs2767060 LD](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45549485788)
- [QTL extension](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45549445784)
- [Gene-level result and CS](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45558846009)
- [Exact identifier audit](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45559096023)
- [Lead diagnostic](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45559352081)
- [Complete curated Catalog regional audit](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45559287128)
