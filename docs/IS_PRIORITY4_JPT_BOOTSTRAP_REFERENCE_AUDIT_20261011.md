# JPT104 reference LD bootstrap audit — IS priority-four loci

**Study date:** 2026-10-11 KST. **Exploratory population-reference audit only.** No formal conditional BBJ GWAS, no in-cohort GTEx LD, no new SuSiE-coloc or inferred causal gene.

## Sources and population sampling

JPT: 104 Japanese in Tokyo, Japan individuals from the existing 1000 Genomes Phase 3 EAS504 panel. Non-JPT EAS: 400 disjoint EAS individuals. EAS504 contains JPT104 (not independent). EUR503: 503 separate reference participants.

Exact GRCh37 genomic keys `chr:pos:REF:ALT` were used and PLINK-counted REF dosages were converted to ALT dosage where needed. Existing local EAS PGEN genotypes and newly extracted full-cis EAS VCF ALT dosages were identical for the lead and candidate SNP in all 504 subjects.

Lead `10:105459834:T:C`, candidate `10:105694487:G:T`. BBJ lead z=−6.130469; candidate marginal z=+4.832637.

| Reference population | N | Signed r | r² | candidate ALT freq | approximate residual z (alpha=0.05) |
|---|---:|---:|---:|---:|---:|
| JPT104 | 104 | **−0.243875** | .059475 | .399038 | **3.5078** |
| non-JPT EAS400 | 400 | −0.124580 | .015520 | .410000 | 4.1362 |
| EAS504 | 504 | −0.147903 | .021875 | .407738 | 4.0111 |
| EUR503 | 503 | **+0.109206** | .011926 | .347913 | **5.4983** |

The approximate residual `z=(z_j−r_alpha*z_lead)/sqrt(1−r_alpha²)`, with `r_alpha=(1−alpha)*r_reference`, is **not a validly calibrated conditional association statistic**. We checked alpha=0,.01,.05,.1,.2 and the pre-existing exploratory |z|≥5.45 screen. Both chr4 and chr10 had **zero** JPT hits at all five alphas. EUR had one chr10 hit at alpha up to .1 and zero at .2.

## Reference panel bootstrap — 2,000 resamples

Paired individual genotypes were sampled with replacement, seed=20261011, with BBJ z values held fixed. Percentile intervals represent **only reference sample LD uncertainty**, not BBJ GWAS sampling variance.

| Population | 95% bootstrap signed r interval | 95% heuristic z interval, alpha=.05 | Resamples meeting |z|≥5.45 |
|---|---|---|---:|
| JPT104 | [−.4223,−.0562] | [2.591,4.511] | 0/2000 |
| Other EAS400 | [−.2203,−.0278] | [3.630,4.673] | 0/2000 |
| EAS504 | [−.2319,−.0626] | [3.570,4.476] | 0/2000 |
| EUR503 | [.0201,.1966] | [4.951,6.084] | 1136/2000 (56.8%) |

Reference LD difference `r_JPT−r_EUR=−.35308` has 95% bootstrap interval [−.5505,−.1527]; `r_JPT−r_otherEAS=−.11929` has interval [−.3243,.0989], which crosses zero. The latter does **not** establish a distinctive JPT-versus-other-EAS difference. 0/2000 JPT crossings is not zero true probability.

## Interpretation

The previously highlighted chr10 EUR-only residual is **not supported by JPT reference LD**. Do not promote it to an independent causal or conditionally genome-wide significant BBJ stroke signal. The comparison is strongly dependent on the reference ancestry and uses 104 JPT individuals only.

The complete original IS discovery universe (2,225 genes), the 646 historical coloc outputs, and the separate 120 EBI–BBJ ABF sensitivity runs remain unchanged. Definitive follow-up requires BBJ in-study conditional GWAS/correct cohort LD and separately validated GTEx QTL cohort LD and effective molecular N.

## Reproducibility

The actual server-side audit was run under `/srv/is-analysis/results/is/audits/is_priority4_jpt_bootstrap_20261011_v1/`; its files include `JPT_REFERENCE_LOCUS_QC.tsv`, `JPT_EAS_EUR_ALPHA_SCREEN.tsv`, `CHR10_JPT_EAS_EUR_CANDIDATE_ALPHA.tsv`, `CHR10_JPT_EAS_EUR_BOOTSTRAP.tsv`, and `JPT_REFERENCE_BOOTSTRAP_MANIFEST.json`. Experimental Python source and local tests were placed in the research worktree, but are **not attested by a successful full-suite test run or GitHub source commit in this step**.

Population background: https://www.internationalgenome.org/data-portal/data-collection/phase3 ; https://test.internationalgenome.org/data-portal/population/JPT .
