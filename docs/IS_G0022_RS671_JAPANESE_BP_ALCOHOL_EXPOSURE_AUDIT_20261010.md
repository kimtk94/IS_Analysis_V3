# IS G0022 rs671 — Japanese blood-pressure evidence and alcohol-GWAS source audit
**Date:** 2026-10-10 KST
**Branch:** `research/is-broad-discovery-20261009`
**Scope:** Continued ALDH2 rs671/chr12 stroke locus mechanistic analysis. **Genetic association shown; causal mediation NOT estimated.**

## 1. Scientific question and current evidence

Does the East-Asian ischemic stroke association at `rs671` (GRCh37 `12:112241766:G:A`) arise through reduced alcohol intake, modified blood pressure, ALDH2 enzyme activity, cis-molecular regulatory effects, or another pathway?

A validated combined GWAS–QTL causal/multivariable mediation estimate **cannot yet be obtained**. Rs671 affects alcohol metabolism directly and changes drinking behavior. The G0022 SuSiE-RSS four SNPs are highly correlated (`r²=0.87–0.99`), not four independent exposure instruments. EAS LD reference n=504; exact SNP-level GWAS effective N unavailable. JCTF molecular cis-eQTL/pQTL full variant universe unavailable; no valid molecular coloc.

Verified GWAS effect (per rs671 alternate A):
- **GIGASTROKE EAS AIS** beta −0.1514 (SE 0.0176), OR=**0.8595**, p=8.426e-18.
- BBJ Japanese ischemic stroke beta −0.10975 (SE 0.01253), OR=0.8961, p=1.99e-18.
- EAS SVS OR=0.856, p=5.87e-7; CES OR=0.791, p=.00156; LAS OR=.912, p=.0814.
These are source study effect estimates, not disjoint cohorts or independent genetic replications. Do not run naive meta-analysis across these overlapping datasets.

## 2. NEW direct public Japanese BP GWAS evidence

A raw official `https://pheweb.jp/variant/12-112241766-G-A` PheWeb HTML source containing the full `window.variant` structured payload was downloaded, saved under `/srv/is-analysis/data/is/gwas_exposure/japanese_bp_sakaue2021/rs671_BBJ_PheWeb.html`, and SHA256-stamped. The parser confirms `chr12:112241766 G>A`, effect allele **ALT A**, study ID `SakaueKanai2021`, and exact phenocode, n, beta, SE and original p.

| BBJ blood-pressure phenotype | Japanese N | Beta (ALT A, original reported scale) | SE | P |
| --- | ---: | ---: | ---: | ---: |
| SBP — systolic BP | **145,505** | **−0.062** | 0.0039 | 1.2e-56 |
| DBP — diastolic BP | **145,515** | **−0.063** | 0.0041 | 9.4e-54 |
| MAP — mean arterial pressure | **145,502** | **−0.067** | 0.0040 | 1.6e-63 |
| PP — pulse pressure | **145,445** | **−0.033** | 0.0039 | 1.3e-16 |

**NEVER interpret these β as raw mmHg until a phenotype-specific unit/normalization check is complete.** BBJ quantitative GWAS effects are presented on study-reported transformed scales. No absolute BP reduction or individual risk prediction is claimed.

Other PheWeb `SBP_GxE` and `DBP_GxE` entries are *combined GxE tests across environments* and contain no beta/SE. They are **not automatically alcohol-interaction tests**; these rows were excluded from directional/mediation analysis.

**Interpretation:** rs671 A is associated with lower BP measures, consistent in direction with lower odds of AIS, but shared samples and pleiotropy can explain this association structure. It is not proof that blood pressure mediates the ALDH2–stroke relationship.

## 3. NEW Japanese alcohol behavior GWAS sources

Official published study: Koyanagi et al., *Science Advances* 2024 (doi:10.1126/sciadv.ade2780). We independently obtained Zenodo record 10038152 metadata and original README, verified official per-file sizes and MD5 contracts, and stored both in the personal server's separate untracked source directory:

`/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024/`

The open v1 sources relevant to direct rs671 exposure-allele tests:

| Open file | Phenotype | Japanese sample N | Expected compressed size | Verified source MD5 |
| --- | --- | ---: | ---: | --- |
| `1_Alcohol_intake_Unstratified.tsv.gz` | Daily consumption, **log2(grams/day + 1)** | 154,570 | 193,465,188 bytes | `a37efa6044605f187d93866d162b2026` |
| `5_Drinking_Unstratified.tsv.gz` | Ever vs never drinking, logistic odds | 175,672 | 193,865,182 bytes | `9d15385620a15520a30de84c0393c0f1` |

**Actual acquisition status as of this analysis:** A source transfer was attempted for alcohol intake but the public endpoint closed the connection after **23,897,084 bytes**, curl code 18. The file remains explicitly `.part`, no MD5 verification, **NO allele-effect output from this incomplete source**. The second drinking-status file has not been acquired. Code `download_g0022_japanese_alcohol_unstratified.py` was changed to **PLAN_ONLY by default** and enforces exact official file sizes and MD5 on explicit `--execute`. `G0022_RS671_ALCOHOL_BP_EXPOSURE_READINESS.tsv` stores the missingness and source sizes. Do not claim that the files were fully downloaded.

The public 2024 manuscript confirms rs671 is a major determinant of alcohol consumption; in its validation dataset the rs671 genotype alone explained **19.22%** of residual variation of `log2(grams/day+1)`, but this is **not a summary-GWAS beta** and must never be used as a substitute for the allele-specific GWAS exposure effect.

Original six Japanese data-cohort meta-analysis includes **BioBank Japan**, J-CGE population cohorts and Nagahama. BBJ GWAS of stroke/BP may share underlying people or clinical selection. `BBJ included` is certain by the source README; *the exact overlap fraction* remains unknown.

Genotype-specific GA or GG GWAS files estimate effects of **other SNPs conditional on fixed rs671 genotype**, while `Interaction` files test SNP×rs671 interactions. These cannot replace rs671's own allele effect from `Unstratified` full-sample exposure GWAS. The source original analyses were sex-adjusted, **not separate male/female stroke-GWAS replication**.

J-CGE downloads also list these alcohol traits but require submitting researcher/affiliation/purpose/email and agreeing to conditions; the open Zenodo source is the compliant access path used here.

## 4. Mediation/MR assumptions and fail-closed decision

### Evidence supported
1. EAS rs671 ALT A is associated with lower AIS odds (GIGASTROKE EAS) and Japanese BBJ BP traits.
2. Japanese rs671-A is well-characterized as loss-of-function ALDH2 missense; original 2024 paper shows strong alcohol consumption impact.
3. JCTF Japanese blood eQTL and measured protein pQTL source results exist for rs671 and its three high-LD GWAS credible-set variants.
4. These are **separate marginal studies/associations**, not formal proof of one shared molecular causal signal or mediation.

### Missing identification
1. Unstratified alcohol / drinker GWAS *source-complete* rs671 β, SE, ALT coding, effect units, exact case definitions: **NOT VERIFIED**.
2. Independent instruments outside rs671's high-LD block and defensible exclusion restrictions. ADH1B/ALDH1B1 variants may be candidates but also directly influence alcohol metabolism/acetaldehyde; cannot assume absence of horizontal pleiotropy.
3. EAS BP GWAS phenotype normalization and sex-specific exposure/outcome effects, and shared BBJ participant covariance.
4. Fine-mapped QTL cis-variant universe and correctly ancestry-matched LD for signal-level colocalization.
5. Replication of GWAS effects in disjoint EAS AIS cohorts.
6. Correct causally identifiable pathway or multivariable MR design, sensitivity to selection, drinking patterns, and direct ALDH2 coding effects.

**Therefore:** no Wald-ratio mediated effect, no IVW/weighted median/MR-Egger/Cochran Q, no BP mediation percent, no therapeutic candidate clinical claim. Additional correlated variants in the same 4-SNP CS do not create independent genetic instruments.

## 5. Reproduce source and code audits

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
python3 scripts/is/audit_g0022_rs671_bbj_bp.py
python3 scripts/is/download_g0022_japanese_alcohol_unstratified.py  # PLAN ONLY
# Explicit consent-free open-data attempt, guarded; public server network reliability varies:
python3 scripts/is/download_g0022_japanese_alcohol_unstratified.py --file alcohol --execute
python3 scripts/is/audit_g0022_rs671_alcohol_bp_readiness.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Canonical GWAS and QTL data are read-only. Derived results:
- `G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv`
- `G0022_RS671_BBJ_BP_PHEWEB_SOURCE_AUDIT.json`
- `G0022_RS671_ALCOHOL_BP_EXPOSURE_READINESS.tsv`
- `G0022_RS671_ALCOHOL_BP_READINESS_SUMMARY.json`

Full locus previous outputs, inherited:
- `G0022_RS671_STROKE_PHENOTYPE_FOREST.png`
- `G0022_RS671_MEDIATION_PATHWAY_READINESS.tsv`
- `G0022_EAS_JCTF_INTEGRATED_VARIANT_EVIDENCE.tsv`

## 6. Next verifiable research milestone

First, acquire and checksum-complete two Japanese exposure GWAS files **then** inspect header_definition effect-allele schema, GRCh37 coordinate, p-values, phenotype units and rs671 measured associations. The original files are full genome gzip TSV, not GRCh37 indexed tabix, so `rs671` effect should be read by safe streaming rather than random access. After harmonization, **do not** compute causal mediation without a defensible IV, nonoverlap, outcome model and pleiotropy design. Pursue independent EAS stroke GWAS/sex-stratification and original QTL/LD source completeness in parallel.

## Primary sources
- BBJ PheWeb rs671: https://pheweb.jp/variant/12-112241766-G-A
- BBJ DBP GWAS source: https://www.pheweb.jp/pheno/DBP
- Koyanagi et al. 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC10816704/
- Official open Zenodo source: https://zenodo.org/records/10038152
- J-CGE conditional access policy: https://epi.ncc.go.jp/cgi-bin/cms/public/index.cgi/jcge/en/download/index
- Korean alcohol flush IV sensitivity: https://pmc.ncbi.nlm.nih.gov/articles/PMC5765011/

## 7. Additional verified input contract (2026-10-10, source-file schema audit)

Successfully fetched the official Zenodo v1 `header_definition_20230615.xlsx` (9,972 bytes; **MD5 `79ebeb63f12317f89b291477cf231a76`**), inspected source `Sheet1` rather than guessing the GWAS columns. Original headers **`SNP, CHR, POS, EA, NEA, EAF, BETA, SE, P, HetP, N`**. `POS` is hg19/GRCh37; `EA` is the effect allele and `NEA` its complement/other allele; `EAF` is frequency of EA. `N` is the per-row total meta-analysis sample size, which need not equal the maximum phenotype-wide N.

Implemented `scripts/is/extract_g0022_koyanagi_alcohol_rs671.py`. It validates **both original file length and MD5** before streaming the gzip TSV, checks the exact schema, matches rs671/three AIS CS variants by GRCh37 chromosome-position-REF/ALT, flips beta/EAF when EA is genomic REF instead of ALT, and reports missing tested variants. It does **not** infer any direct/indirect causal effects. For drinking status, the source binary outcome **case coding requires manual confirmation**; no unverified claim that its `BETA` measures ever-drinking odds is made.

When both source files are fully retrieved and source MD5 passes, the runnable extraction steps are:

```bash
python3 scripts/is/extract_g0022_koyanagi_alcohol_rs671.py --exposure alcohol_intake
python3 scripts/is/extract_g0022_koyanagi_alcohol_rs671.py --exposure drinking_status
python3 scripts/is/audit_g0022_rs671_alcohol_bp_readiness.py
```

Until then, the intended extraction must fail on the existing incomplete `.part` file; this is **correct fail-closed behavior**, not an absence of rs671 alcohol GWAS association.
