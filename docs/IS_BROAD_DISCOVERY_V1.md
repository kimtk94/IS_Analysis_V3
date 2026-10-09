# IS Broad Discovery v1 — expansion contract

## Rationale
The previously selected BBJ loci L001-L004, 9 manually selected genes, and positional nearest-five rule are a **deep-validation anchor**, not an exhaustive ischemic stroke discovery universe.

## Data currently accessible
- BBJ ischemic stroke: canonical GRCh37 variant-level GWAS
- GIGASTROKE EAS: AS, AIS, CES, LAS, SVS canonical GRCh37 GWAS
- GIGASTROKE full cross-ancestry, GBMI, independent replication: **not claimed loaded** until exact files are inventoried and QC verified.

## Phase 0 exploratory screen
Run: `python3 scripts/is/build_broad_discovery_v1.py`

- Include quality-filtered p<=5e-8 and an independent suggestive analysis p<=1e-6.
- Do not drop any genome-wide significant loci because they are absent from the historic four-locus study.
- Exclude variants with INFO<0.7 or MAF<0.01 from the primary screen; retain separately in a QC-exclusion ledger in next revision.
- Current one-megabase positional clustering is **provisional**, not LD clumping.
- Output is under `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1`.
- Existing upstream canonical GWAS and Phase9 results remain untouched.

## Required Phase 1 before claiming unbiased locus-to-gene mapping
1. Source manifest with study, phenotype/subtype, population, sample sizes, genome build, overlap and QC.
2. LD-independent loci, fine-mapping and all credible-set variants; ancestry-matched LD and allele checks.
3. Include all protein-coding and noncoding genes in the credible-set regions plus distal targets linked by molecular QTL/ATAC/chromatin interactions; use distance as a score, not an exclusion.
4. Colocalization (multi-signal when supported); human and animal evidence separately labeled.
5. Audit universe→tested→excluded→unresolved→ranked counts with exclusion reasons.
6. Rank whole candidate universe. Deep validation capacity may be limited, discovery is not.
7. Separate discovery and replication sources and assess cohort overlap.
8. Record all caveats: absence of expression or a missing resource does **not** prove a locus has no causal effect.

## Stable anchor
FGF5, ALDH2, SH3PXD2A, COL4A2, COL4A1, C4orf22, CALHM2, NEURL1 and INA remain benchmark/anchor genes. Do not grant them automatic ranking bonuses.

## Claim discipline
A provisional region is not a new independent locus. A gene annotated by position is not a causal gene. A suggestive signal is not genome-wide significant. Existing discovery cohorts must not be falsely presented as independent replication.
