# IS G0022 AS/AIS LD 불일치 원인 감사 — V2
**일자:** 2026-10-09 (KST)  
**브랜치:** `research/is-broad-discovery-20261009`  
**상태:** READONLY molecular-causality gate. Canonical GWAS 및 2,225개 후보 유전자 universe 불변.

## Executive summary

기존 G0022(chr12) SuSiE-RSS pilot에서는 EAS AS `GCST90104544`의 두 영역에서 심한 GWAS–1000G EAS LD 불일치가 관찰됨(`estimate_s_rss` `s≈0.278/0.178`); EAS AIS `GCST90104545`의 세 영역에서는 `s≈0.006–0.011`. 이에 대해 SNP allele flip, phenotype effect 방향, population reference 차이, SNP coverage/QC를 단계적으로 검사했음.

**핵심 발견:** AS와 AIS에 공통으로 존재하는 변이만을 *동일한 ALT 방향·동일한 LD 행렬*에 적용하면, AS 두 영역의 `s`가 `0.0065/0.0100`으로 떨어짐. 반면 SNP 수를 동일하게 줄인 20회의 **무작위 제거 대조군**에서는 불일치가 지속됨(첫 영역 `s=0.248–0.295`, 둘째 `s=0.151–0.193`). AS pilot에만 수록된 35/43개, 총 **78개 변이는 AIS canonical GWAS에서 해당 염색체 위치가 아예 없음**. 이들은 모두 단독으로 GWS 수준을 충족하지 않으며, 최대 `|Z|=2.646`임.

이 결과는 **GWAS 간 marker coverage 차이가 mismatch에 연관됨**을 지지함. 다만 그 78개 각각이 오류라는 근거는 없으며, AS case mix, meta-study의 per-SNP 참여 코호트, imputation quality, phenotype 구성 등의 실제 원인까지는 규명하지 못함.

## 데이터 출처 및 보호 조건

- GIGASTROKE EAS AS: `GCST90104544`, 27,413 cases + 237,242 controls = 264,655명. Source: https://kp4cd.org/node/1138 .
- GIGASTROKE EAS AIS: `GCST90104545`, 19,032 cases + 237,242 controls = 256,274명. Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC13621555/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC12668099/ .
- Source total population verified against original GWAS Catalog accession GRCh37 YAML, but **variant-specific effective sample sizes are unavailable**.
- Analytical sensitivity working `N_eff=4×N_cases×N_controls/(N_cases+N_controls)`: AS `98,294`, AIS `70,474`. This is a **study-level rough effective N**, not a validated meta-GWAS per-SNP effective N.
- All pilot genotypes: 1000 Genomes EAS Phase 3, **504 reference individuals**, expanded chr12 G0022 stable-ID GRCh37 panel. Windows are **±250kb** around five clump indexes (NOT full joint locus fine mapping).
- The AS/AIS studies contain related/overlapping cohorts and phenotypes; their concordance cannot be called external independent replication.

## 단계별 실제 검증

### 1. GWAS Z 방향 및 빈도 일관성
G0022 region-wide ALT harmonized GWAS에서 **5,254 common variants** 검출. AS 2개 clump window의 AS/AIS Z-score Pearson `r=0.9758` / `r=0.9900`. 두 분석에서 `|Z|≥3`인 변이끼리 서로 반대 부호인 사례 0건. 단순 일괄 allele-flip 가설은 지지되지 않음.

### 2. SuSiE-RSS kriging 기반 outlier
`kriging_rss`에서 `logLR>2` **그리고** `|Z|>2`인 변이 0건(5개 window 모두). Kriging이 미검출한 non-isolated LD mismatch는 존재할 수 있으므로 **부정 증명 아님**.

### 3. MAF/EAF filter sensitivity
`estimate_s_rss` `s` (pilot local LD; effective N 근사):

| EAS AS window | Baseline MAF≥0.01 / ΔEAF≤0.15 | Strict MAF≥0.05 / ΔEAF≤0.05 | Very strict MAF≥0.10 / ΔEAF≤0.03 |
| --- | ---: | ---: | ---: |
| 111629389 | 0.27793 | 0.20937 | 0.15865 |
| 112930475 | 0.17794 | 0.22918 | 0.09641 |

AS 두 영역 모두 단순 EAF/MAF 강화만으로 불일치 해소되지 않음. Strict filter에 따른 variant count/LD 구조 변경은 별도의 과학적 sensitivity일 뿐 'correction'으로 간주하지 않음.

### 4. 동일 SNP 집합 + 동일 signed-LD 비교
AS windows에서 AIS 원본과 allele-matched인 SNP만 유지. **동일 1000G EAS signed-LD 행렬**을 AS와 AIS 요약통계에 공통 적용함:

| AS window center | Full AS SNP | Matched common SNP | Full AS `s` | Common-SNP AS `s` | Common-SNP AIS `s` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 111629389 | 459 | **424** | 0.27793 | **0.00651** | 0.00823 |
| 112930475 | 529 | **486** | 0.17794 | **0.01001** | 0.01197 |

이는 전체 AS phenotype이 무조건 잘못되었다는 주장과 반대되는 근거임. 다만 common-SNP selection 자체가 분석 대상과 LD 차원을 변경했으므로 해당 서브셋을 최종 결론으로 사용하면 안 됨.

### 5. 임의 변이 제거 대조 (20회/영역)
기존 AS SNP 집합에서 공통 SNP만큼 동일 개수로 무작위 추출. `set.seed(20261009)` fixed.

| AS window | Common-sample subset `s` | 20개 무작위 subset `s` range | 무작위가 common `s` 이하로 떨어진 횟수 |
| --- | ---: | ---: | ---: |
| 111629389 | 0.00651 | 0.24840–0.29529 | **0/20** |
| 112930475 | 0.01001 | 0.15072–0.19258 | **0/20** |

단순 SNP 수 감소가 아닌 **특정 변이 집합 포함 여부**의 영향을 지지함. 이 반복 횟수는 QC 대조이지 multiple-testing-adjusted formal p-value 검정이 아님.

### 6. AIS canonical GRCh37 전체 원본과 78개 AS-only 변이 대조
AIS의 6M+ canonical variant rows를 직접 다시 읽어, AS pilot에서 AIS matched summary에 없던 총 78개 변이의 chr12 position을 조회함:
- **78/78 모두 `NO_AIS_CANONICAL_POSITION`**
- 원본 위치는 있으나 allele이 달라진 변이: 0
- 원본 위치와 allele pair는 있으나 harmonization 실패한 변이: 0
- 해당 AS-only 78개 중 GWS: **0**, 최대 `|Z|≈2.646`

이것이 공개 AIS GWAS에서 SNP가 누락된 이유까지 증명하지는 않음. 메타연구 참여 수, imputation 정보, 연구별 QC/필터 자료가 필요.

## Effective-N SuSiE 재실행

Source-supported **근사 case/control effective N** 적용:
- AIS 3개 local windows: 3/3 수렴, `Xcorr=R`, purity≥0.5 적용 95% local CS 각 1개 유지.
- AS 2개 local windows: 1/2 수렴; 원래 높은 LD mismatch 해소되지 않음. 설령 CS/PIP가 크더라도 확정 근거로 사용 불가.

최종 safe gate: **AS 2개 BLOCKED**, **AIS 3개 EXPLORATORY**, **validated causal signals/genes 0**. 상위 positional genes ALDH2, PTPN11, IFT81 등은 아직 QTL·coloc 및 기능적 검증 필요.

## 재현 방법

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

python3 scripts/is/prepare_g0022_effective_n.py
Rscript scripts/is/run_g0022_effective_n_sensitivity.R
python3 scripts/is/audit_g0022_as_ais_z_concordance.py
Rscript scripts/is/diagnose_g0022_kriging_outliers.R
Rscript scripts/is/sensitivity_g0022_ld_variant_filters.R
python3 scripts/is/prepare_g0022_matched_as_ais_ld.py
Rscript scripts/is/diagnose_g0022_shared_snp_as_ais_ld.R
Rscript scripts/is/sensitivity_g0022_as_missing_vs_random.R
python3 scripts/is/audit_g0022_as_only_source.py
python3 scripts/is/audit_g0022_diagnostic_decision_v2.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

**Read-only canonical**; all derived files are under:

`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot/`

Principal outputs:
- `G0022_DIAGNOSTIC_DECISION_V2_SUMMARY.json`
- `G0022_DIAGNOSTIC_DECISION_V2.tsv`
- `G0022_AS_ONLY_SNP_CANONICAL_SOURCE_AUDIT.tsv`
- `G0022_AS_ONLY_SNP_CANONICAL_SOURCE_SUMMARY.json`
- `matched_as_ais_ld_sensitivity/G0022_SHARED_SNP_AS_AIS_LD_DIAGNOSTIC.tsv`
- `matched_as_ais_ld_sensitivity/G0022_AS_MATCHED_SNP_VS_RANDOM_DROPOUT.tsv`
- `G0022_KRIGING_SNP_DIAGNOSTICS.tsv`
- `effective_n_sensitivity/G0022_SUSIE_EFFECTIVE_N_SENSITIVITY.tsv`

### 향후 우선순위

1. **GIGASTROKE source cohort/meta-specific sample sizes / INFO 및 imputation QC** 확보. AS-only 78 SNP이 왜 AIS에서 빠졌는지 연구별 contributing cohorts와 비교, 원본 VCF/INFO/allele strand 검토. 별도 검증 전 AS 원자료 일부를 일괄 삭제하거나 오류로 선언하지 말 것.
2. **AIS 3개 locus pilot의 full-region 4Mb multi-signal fine mapping**을 수행하되 valid N, LD r~ reference and numerical convergence guardrails 선행. 504 EAS 참조 표본 수와 강한 LD로 overconfident PIP 위험 존재.
3. QTL cis e/s/pQTL, brain/vascular endothelial pericytes/human IS tissue single-cell, enhancer-to-gene 연결을 통해 chr12 signal별 causal gene 후보 대조.
4. 전체 80 ancestry/study provisional areas / 2,225 unique positional genes 유지. 지리적 중첩·clumping을 independent causal loci로 해석하지 않음.

## 방법론 참고

- SuSiE-RSS 공식 diagnostics: https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html
- LD mismatch/finite-reference issue: https://stephenslab.github.io/susieR/articles/rss_mismatch.html
- GIGASTROKE original paper: https://www.nature.com/articles/s41586-022-05165-3
