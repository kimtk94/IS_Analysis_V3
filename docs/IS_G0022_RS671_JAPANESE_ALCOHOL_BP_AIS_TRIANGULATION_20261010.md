# IS G0022: Japanese rs671 alcohol GWAS recovered, ALT harmonization and AIS/BP evidence triangulation
**Research snapshot:** 2026-10-10 KST
**Repository/worktree:** `research/is-broad-discovery-20261009` / `/srv/is-analysis/worktrees/is-broad-discovery-20261009`
**Status:** **REAL OPEN JAPANESE GWAS EXPOSURE ASSOCIATIONS VERIFIED; CAUSAL MEDIATION NOT IDENTIFIED**

## Executive result

Two full, official, public unstratified Japanese alcohol GWAS source files have now been acquired on the user-owned server after earlier non-resumable API downloads terminated. A new reproducible downloader targets Zenodo's **official record-file HTTP Range endpoint**, requests bounded 4-MiB chunks, validates **HTTP 206 Content-Range and byte counts**, assembles, then requires the full official **MD5 checksum** before exposing a source file. **Both source MD5 checks PASSED.**

- Alcohol amount GWAS, Koyanagi et al. 2024: **7,676,853 original genome-wide variant rows**, N=154,570 at rs671.
- Ever-vs-never drinking status GWAS: **7,679,425 original genome-wide variant rows**, N=175,672 at rs671.
- Combined **15,356,278 source rows** actually streamed and audited.
- All **four G0022 AIS GWAS 95% credible set variants** matched each source by **GRCh37 chromosome, position, REF/ALT**, with effect direction harmonized to the GWAS **ALT allele**.
- rs671 allele A is strongly associated with **lower log2 daily alcohol intake**. The drinking status model also has negative A-coded beta, but binary case/event coding is **not yet confirmed** in the official source definition; never infer an ever-drinking odds decrease until case encoding is verified.
- The Japanese BBJ rs671 A allele is associated with lower BP phenotypes and GIGASTROKE EAS rs671 A is associated with lower AIS odds. **This is descriptive association triangulation and not proof of causal mediation**. Potential BBJ cohort overlap and major ALDH2 coding pleiotropy remain.

## 1. Provenance: original Zenodo open study dataset

Official source: https://zenodo.org/records/10038152 (v1, DOI 10.5061/dryad.tmpg4f546). The associated 2024 *Science Advances* paper discusses rs671-stratified analysis; **ONLY the two unstratified sources** are appropriate for direct rs671 genotype effect.

| File | Official size (bytes) | Official MD5 | Download / QC |
|---|---:|---|---|
| `1_Alcohol_intake_Unstratified.tsv.gz` | **193,465,188** | `a37efa6044605f187d93866d162b2026` | **PASS** |
| `5_Drinking_Unstratified.tsv.gz` | **193,865,182** | `9d15385620a15520a30de84c0393c0f1` | **PASS** |
| `header_definition_20230615.xlsx` | 9,972 | `79ebeb63f12317f89b291477cf231a76` | **PASS** |

An earlier partial API transfer ended at 23,897,084 bytes with curl error 18, **never used for statistics**. The successful source was assembled using `https://zenodo.org/records/10038152/files/{filename}?download=1`; this endpoint supports HTTP 206 and validated `Content-Range`, unlike Zenodo API file-content download in the failed attempt. Both datasets required 47 finite byte-range blocks. No synthetic effect has been used.

The official `header_definition_20230615.xlsx` defines `SNP CHR POS EA NEA EAF BETA SE P HetP N`. `POS` is hg19 **GRCh37**, `EA` effect allele, `NEA` non-effect allele, `EAF` frequency of EA, `HetP` between-cohort heterogeneity test, `N` variant-specific meta-GWAS sample count.

**Version note:** Zenodo exposes a 2023 early v1 dataset; Dryad additionally displays later versions with much larger files and altered metadata. This analysis is **explicitly tied to the v1 source**. Verify any version differences before final manuscript publication; do not silently mix data between releases.

## 2. Real G0022 four-SNP ALT-coded exposure statistics

The source p-value is **numerical zero for all 8 strong association rows** (underflow in the released file). This **does not mean probability exactly zero**. Infer strength from `beta / SE`, not (p=0). Do not invent precise p-values where source rounds or truncates them.

### Daily alcohol consumption GWAS: log2(grams/day + 1)

| AIS CS variant GRCh37 | Source EA | ALT allele | β_ALT | SE | N | HetP |
|---|---|---|---:|---:|---:|---:|
| `12:112168009:G:A` | A | A | **−1.3168** | 0.0074 | 154,570 | 1.35e−68 |
| **`12:112241766:G:A` (ALDH2 rs671)** | A | A | **−1.3215** | 0.0074 | 154,570 | **1.15e−62** |
| `12:112468206:C:T` | T | T | **−1.4267** | 0.0081 | 154,570 | 4.68e−71 |
| `12:112736118:A:G` | A | **G** | **−1.4580** | 0.0086 | 154,570 | 3.16e−72 |

**Critical check:** The last SNP has `EA=A` but the GWAS ALT is G. The original exposure β was explicitly **negated** and EAF complemented (`1−EAF`) to align with the ALT effect allele. Without this step the apparent direction would be wrong.

These effect sizes use log2(grams/day+1), not raw grams/day. One A-allele association on this transformed scale must not be presented as a fixed gram/day absolute reduction. The heterogeneity test is very significant across study cohorts and precludes treating one Japanese meta summary as a homogeneous experimental perturbation.

### Drinking status GWAS: binary never/ever status

| AIS CS variant GRCh37 | Source EA | ALT allele | β_ALT | SE | N | HetP |
|---|---|---|---:|---:|---:|---:|
| `12:112168009:G:A` | A | A | −1.8169 | 0.0121 | 175,672 | 0.00518 |
| **`12:112241766:G:A` (rs671)** | A | A | **−1.8264** | 0.0121 | 175,672 | **0.04178** |
| `12:112468206:C:T` | T | T | −1.9636 | 0.0131 | 175,672 | 4.70e−5 |
| `12:112736118:A:G` | A | **G** | −1.9874 | 0.0136 | 175,672 | 1.69e−5 |

The status variable is documented as never versus ever, but the precise case coding of the binary logistic model is not explicitly resolved by the file header. `BETA_ALT<0` is **a model direction only** until event coding is confirmed. All four p-values are recorded as numerical zero (extreme p underflow), not interpretable as literal zero probability.

## 3. rs671-A exposure/BP/AIS association triangulation

All effects are aligned to **A at GRCh37 12:112241766 G>A**; source β units differ and are never combined by arithmetic:

| Japanese / EAS summary | Trait scale | β_A | SE | Interpretation |
|---|---|---:|---:|---|
| Koyanagi et al. Japanese unstratified amount | log2(grams/day+1) | **−1.3215** | 0.0074 | A associated with lower daily intake |
| Koyanagi et al. drinking status | Binary logistic, case code unresolved | **−1.8264** | 0.0121 | Negative model coefficient; event interpretation blocked |
| Sakaue/Kanai BBJ systolic blood pressure | Original GWAS standardized/transformed BP trait | **−0.0620** | 0.0039 | A associated with lower reported BP phenotype |
| GIGASTROKE EAS AIS | Case-control log odds | **−0.1514** | 0.0176 | A associated with **OR=0.8595** per allele |

Additional Japanese BBJ DBP β_A −0.063, MAP β_A −0.067, PP β_A −0.033. Original BP phenotype normalization has not been reverse-transformed into mmHg. The BBJ source has possible cohort overlap with the alcohol source, and GIGASTROKE EAS AIS may reuse BBJ individuals. **These are not four independent replications.**

## 4. Scientific gate: what is supported and what is not

**Supported:**
- ALDH2 rs671-A shows a robust (in-sample) negative **Japanese alcohol-intake** association and a negative BBJ BP association.
- GIGASTROKE EAS stroke shows rs671-A lower AIS odds, with ancestry-compatible allele direction.
- JCTF Japanese blood QTL shows rs671 and all three high-LD GWAS CS variants associated with ALDH2 transcript and Olink-measured protein.
- The 4-SNP AIS CS correlates strongly with rs671 in 1000G EAS n=504; external-LD residualization showed no additional GWS in a **single-lead approximate** scan, not a formal conditional GWAS.

**Not supported / blocked:**
- Ischemic stroke prevention from reduced alcohol consumption, or a `mediated percentage`. ALDH2 rs671 is a functional missense variant with strong direct acetaldehyde/aldehyde metabolism impacts, alcohol drinking behavior changes, and other potential pathways. Its use as an alcohol-only IV violates a naïve exclusion-restriction assumption.
- Formal IVW, MR-Egger, weighted median and Cochran Q using four highly correlated G0022 CS SNPs as if independent.
- `coloc.susie` without full, ancestry-matched, unfiltered QTL summary statistics and LD.
- Quantitative comparison between β values on different scales; BBJ blood-pressure GWAS effects may be standardized and not raw mmHg.
- Independent EAS AIS replication, alcohol-by-sex stroke effect modification, true direct enzyme→stroke effect.
- Any conclusion that nominally negative LAS result proves stroke subtype differences (overlapping samples; covariance unverified).

**Proposed DAG:** rs671 missense influences ALDH2 enzyme activity and drinking tolerance; these may influence alcohol intake and BP, while enzyme activity may independently change oxidative/aldehyde vascular stress; blood eQTL/pQTL could reflect cis regulation and/or epitope binding. Sex, environment and source case-control sampling affect the relationships. Without separate instruments and mediator identification, no path-specific causal estimate is defensible.

## 5. Reproducible scripts and integrity controls

The scripts are on the isolated research worktree; data are outside Git:

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
# PLAN only (no network by default):
python3 scripts/is/download_g0022_japanese_alcohol_ranges_parallel.py
# Execute only if originals absent:
python3 scripts/is/download_g0022_japanese_alcohol_ranges_parallel.py --file both --workers 8 --block-mib 4 --execute

# Complete-source exact MD5 and allele orientation, streaming >15 million rows:
python3 scripts/is/extract_g0022_koyanagi_alcohol_rs671.py --exposure alcohol_intake
python3 scripts/is/extract_g0022_koyanagi_alcohol_rs671.py --exposure drinking_status

python3 scripts/is/audit_g0022_rs671_alcohol_bp_readiness.py
python3 scripts/is/audit_g0022_rs671_allele_triangulation.py
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Original verified sources:
`/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024/`

Outputs under:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/`

- `G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv`
- `G0022_4CS_KOYANAGI2024_DRINKING_STATUS_EXPOSURE.tsv`
- `G0022_KOYANAGI2024_ALCOHOL_INTAKE_AUDIT.json`
- `G0022_KOYANAGI2024_DRINKING_STATUS_AUDIT.json`
- `G0022_RS671_ALCOHOL_BP_EXPOSURE_READINESS.tsv`
- `G0022_RS671_ALCOHOL_BP_AIS_ALLELE_TRIANGULATION.tsv`
- `G0022_RS671_ALCOHOL_BP_AIS_TRIANGULATION_SUMMARY.json`

### Next legitimate scientific milestones
1. Confirm the case coding for the binary drinking-status source. Do not infer a reduced probability of ever drinking from the negative β until this is resolved.
2. Query a truly independent Japanese/East-Asian ischemic stroke GWAS, identify potential sample reuse with BBJ alcohol and BP cohorts, and examine ancestry-specific imputation quality for rs671.
3. Investigate the **substantial between-cohort heterogeneity** of alcohol amount (rs671 HetP≈1.15e−62). Check study cohort drinking definitions, sex structure and covariates.
4. Acquire sex-stratified EAS alcohol / BP / stroke outcome summaries and independent *non-ALDH2* instruments with falsification controls for separate drinking vs acetaldehyde pathways.
5. Obtain full unfiltered QTL cis summary and harmonized LD to run multi-signal colocalization; test ALDH2 Olink protein binding artifact in a separate assay before inferring abundance mediation.
6. Preserve broad 2,225-gene IS discovery; G0022 is one mechanistic case study, not an exclusive shortlist.

## Primary references

- Official open Zenodo (v1): https://zenodo.org/records/10038152
- Dryad original and later versions: https://datadryad.org/dataset/doi%3A10.5061/dryad.tmpg4f546
- Science Advances original: https://www.science.org/doi/10.1126/sciadv.ade2780
- BBJ rs671 blood pressure: https://pheweb.jp/variant/12-112241766-G-A
- Japan Omics Browser rs671 molecular evidence: https://japan-omics.jp/variant?input_value=rs671
- GIGASTROKE: https://www.nature.com/articles/s41586-022-05165-3
