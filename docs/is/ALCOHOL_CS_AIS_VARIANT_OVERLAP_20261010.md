# Exploratory alcohol fine-mapped variants vs. AIS GWAS — 2026-10-10

## Status / boundaries

**DESCRIPTIVE ASSOCIATION OVERLAP ONLY.** Not MR, mediation, colocalization, shared causality, or independent replication. The original Koyanagi alcohol GWAS official assembly declaration has not yet been independently confirmed. Eleven deduplicated coordinate:REF:ALT variants matched the on-server GRCh37 FASTA REF allele (11/11). The Ensembl GRCh37 VEP colocated rsID plus positional allele-set audit matched 11/11; it is not independent dbSNP evidence. EAS 1000G (n=504) is proxy LD, not matched GWAS LD; SuSiE R reliability/sensitivity warnings remain. ALDH2 fine-mapping BLOCKED, ADH1B EXPLORATORY.

## Direct AIS summary statistics (aligned to verified-input ALT)

| Gene | rsID | BBJ Japanese IS p | GIGASTROKE EAS AIS p | GIGASTROKE variant allele QC |
|---|---|---:|---:|---|
| ALDH2 | rs671 | 1.99193648e-18 | 8.426e-18 | EXACT |
| ALDH2 | rs184590119 | MISSING | MISSING | MISSING |
| ADH1B | rs1229984 | 0.517197 | 0.3659 | REF_ALT_SWAP_REVIEW |
| ADH1B | rs34144181 | 0.754749 | 0.5823 | EXACT |
| ADH1B | rs12502498 | 0.577723 | 0.1665 | REF_ALT_SWAP_REVIEW |
| ADH1B | rs9997653 | 0.586082 | 0.1656 | EXACT |
| ADH1B | rs1826906 | 0.591041 | 0.1618 | REF_ALT_SWAP_REVIEW |
| ADH1B | rs1442485 | 0.612576 | 0.1721 | EXACT |
| ADH1B | rs12502290 | 0.587270 | 0.1850 | EXACT |
| ADH1B | rs2584461 | 0.897024 | 0.5094 | REF_ALT_SWAP_REVIEW |
| ADH1B | rs17028965 | 0.968565 | 0.5746 | EXACT |

**Totals:** BBJ exact 10/11, absent 1/11. GIGASTROKE exact 6/11, opposite REF/ALT orientation but identical allele sets 4/11, absent 1/11. A swapped canonical allele pair is **not** FASTA REF-verified, even when GWAS effect and other allele permit descriptive sign harmonization. All four GIGA discrepancies are T:C versus C:T at the same coordinate.

## Exposure–outcome allelic direction

- rs671 (12:112241766:G:A): original Japanese alcohol GWAS beta_ALT=-1.3215 for log2(grams/day+1); BBJ IS beta_ALT=-0.1097493721, p=1.99193648e-18; GIGASTROKE EAS AIS beta_ALT=-0.1514, p=8.426e-18. Consistent descriptive signs, no demonstrated mediation.
- rs1229984 (4:100239319:T:C): alcohol beta_ALT=+0.1680; BBJ IS beta_ALT=+0.0083996654, p=0.5172; GIGASTROKE AIS beta_ALT=+0.0124, p=0.3659, flagged REF_ALT_SWAP_REVIEW. Lack of single-SNP AIS association does not exclude any region or effect.
- rs2584461 (4:100321573:T:C): alcohol beta_ALT=-0.1048; BBJ IS beta_ALT=+0.00216995, p=0.897; GIGA beta_ALT=+0.0116, p=0.5094, flagged REF_ALT_SWAP_REVIEW. Null AIS associations do not support interpreting direction.
- Alcohol source p-values equal to literal 0 in input reflect numeric underflow/representation; use source beta and SE, not p=0, for effect-size interpretation.

## Sources / generated output

- BBJ input: `/srv/is-analysis/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz`
- GIGASTROKE input: `/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz`
- Alcohol source: `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1/{ADH1B,ALDH2}/variants.tsv`
- Relative output root: `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/susie_rss_finite_ref_sandbox_v1/`
- `ais_overlap_v2/ALCOHOL_CS_AIS_GWAS_DIRECT_OVERLAP.tsv`
- `ais_overlap_v2/ALCOHOL_CS_AIS_GWAS_DIRECT_OVERLAP_MANIFEST.json`
- `ais_effect_direction_v1/ALCOHOL_CS_AIS_ALIGNED_EFFECTS_EXPLORATORY.tsv`
- `ais_effect_direction_v1/ALCOHOL_CS_AIS_EFFECT_DIRECTION_MANIFEST.json`
- Scripts: `scripts/is/audit_alcohol_cs_ais_direct_overlap.py`, `scripts/is/summarize_alcohol_cs_ais_effect_directions.py`
- Four isolated regression tests passed for the AIS direct-overlap parser/harmonizer.

## Next required validations

1. Explain and repair, **in a derivative artifact only**, four GIGASTROKE EAS canonical REF/ALT swaps by comparing source records, provenance and GRCh37 FASTA; run a wider reference integrity audit before locus-level inference.
2. Independently verify Koyanagi original GWAS assembly/build; resolve study sample overlap between BBJ IS and GIGASTROKE EAS AIS before treating them as independent replication.
3. Reassess locus-wide ALDH2/ADH1B conditional signals with appropriately matched LD and sensitivity tests. Do not infer alcohol causality or distinguish mediation versus pleiotropy from present associations.
4. Maintain ADH1B EXPLORATORY and ALDH2 BLOCKED fine-mapping status pending LD and conditional reliability resolution.
