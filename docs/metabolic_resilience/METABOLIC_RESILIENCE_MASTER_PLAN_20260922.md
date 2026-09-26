# Metabolic Resilience 연구 설계 마스터 문서
## Proteogenomic determinants of metabolic resilience: cross-ancestry causal inference and longitudinal validation in a Korean population

**Version:** 0.1  
**Date:** 2026-09-22  
**Primary purpose:** 성균관대학교 Digital Health 석사 연구 후보 주제의 실제 분석 착수용 설계서  
**Status:** Feasibility / Stage 0 설계 완료, 데이터 audit 착수 전  
**Primary concept:** 질병 발생 이후가 아니라, 노화와 비만에도 불구하고 유리한 대사 상태를 장기간 유지하는 **metabolic resilience**의 유전적·단백질적 결정요인을 규명한다.

---

# 0. Executive summary

## 0.1 연구 질문

> **왜 어떤 사람은 나이가 들고 비만도가 증가해도 장기간 양호한 대사건강을 유지하는 반면, 다른 사람은 빠르게 metabolic syndrome, type 2 diabetes 및 cardiometabolic disease 방향으로 악화되는가?**

본 연구는 이 질문을 다음의 세 축으로 분해한다.

1. **Proteogenomic causal screening**
   - 대규모 plasma pQTL을 exposure로 사용한다.
   - glucose, lipids, adiposity, blood pressure, T2D 등의 대규모 GWAS를 outcome으로 사용한다.
   - cis-pQTL 기반 MR과 colocalization을 통해 대사건강에 일관되게 유리한 candidate proteins를 찾는다.

2. **Cross-ancestry replication**
   - European discovery에서 발견한 protein-trait 관계가 East Asian populations에서도 재현되는지 평가한다.
   - CKB East Asian pQTL 및 KoGES / BBJ / TWB summary statistics를 이용한다.
   - ancestry-specific LD와 allele frequency 차이를 반영한다.

3. **Korean longitudinal validation**
   - KoGES 안산·안성 지역사회 코호트의 반복 측정 데이터를 이용한다.
   - 단순한 baseline disease 여부가 아니라 약 20년 동안의 metabolic trajectory를 정의한다.
   - 후보 protein의 cis-pQTL 기반 genetically predicted protein score와 장기 metabolic resilience 간 연관성을 검증한다.
   - physical activity, obesity, sleep, diet 등을 modifier로 평가한다.

---

# 1. 왜 CKD 대신 이 주제를 고려하는가

기존 CKD 연구 설계의 장점은 명확하다.

- outcome이 명확하다.
- CKD/eGFR GWAS가 매우 크다.
- pQTL → MR → colocalization 구조가 잘 정립되어 있다.
- kidney single-cell annotation이 가능하다.
- KoGES에서 eGFR slope를 longitudinal하게 검증할 수 있다.

그러나 plasma proteomics와 CKD를 연결한 MR/colocalization 연구는 이미 상당수 존재한다. 따라서 CKD 연구에서 실제 novelty는 다음에 집중된다.

- East Asian replication
- Korean longitudinal eGFR decline
- ancestry heterogeneity
- 한국인 tissue/functional validation

반면 metabolic resilience는 질병 하나를 분석하는 대신 다음 개념을 다룬다.

> **질병 발생 이전 단계에서 건강 상태가 유지되거나 붕괴되는 과정**

따라서 Digital Health의 다음 키워드와 직접 연결된다.

- Wellness
- Prevention
- Healthy aging
- Precision health
- Lifestyle modification
- Longitudinal monitoring
- Personalized intervention

KoGES 역시 최근 cohort profile update에서 건강노화와 장기추적 연구 자원으로서의 방향을 강조하고 있어 연구 주제와 시점이 잘 맞는다.

---

# 2. 핵심 연구 가설

## H1. Proteomic causal architecture

일부 circulating proteins는 단순한 metabolic disease biomarker가 아니라 다음 여러 metabolic domains에 유전적으로 일관된 영향을 준다.

- adiposity
- glycemia
- dyslipidemia
- blood pressure
- type 2 diabetes

즉 특정 단백질은 하나의 질환이 아니라 **metabolic homeostasis 전체에 영향을 미치는 molecular determinant**일 수 있다.

---

## H2. Cross-ancestry reproducibility

European ancestry에서 발견되는 protein-metabolic trait causal signal 중 일부는 East Asian population에서도 방향과 효과가 보존된다.

반대로 일부 signal은 다음 이유로 ancestry-specific할 수 있다.

- allele frequency 차이
- LD structure 차이
- ancestry-enriched pQTL
- 환경 및 생활습관 차이
- proteomic assay/platform 차이

따라서 단순한 EUR 연구보다 cross-ancestry 검증이 필요하다.

---

## H3. Longitudinal metabolic resilience

baseline metabolic health가 유사하더라도 장기적으로 다음과 같은 이질적인 trajectory가 존재한다.

1. Stable healthy / resilient
2. Slowly deteriorating
3. Rapidly deteriorating
4. Persistent high metabolic burden

후보 protein의 유전적 proxy가 이러한 trajectory와 연관될 수 있다.

---

## H4. Lifestyle modification

candidate protein의 효과는 lifestyle과 독립적이지 않을 수 있다.

특히 다음 interaction을 가정한다.

\[
Protein_{genetic} \times PhysicalActivity
\rightarrow MetabolicTrajectory
\]

추가 exploratory modifier:

- diet quality
- smoking
- alcohol
- sleep
- BMI/central adiposity
- sex
- age
- menopause status where appropriate

---

# 3. 연구의 novelty를 어디에 둘 것인가

## 3.1 피해야 할 단순한 연구

다음 설계는 선행연구와 차별화가 약하다.

### A. Proteomics → prevalent metabolic syndrome

이미 plasma proteomics와 prevalent/incident metabolic syndrome을 분석하고 MR까지 수행한 연구가 있다.

### B. Proteomic BMI score → T2D/CVD

UK Biobank에서 BMI-associated proteomic score와 incident T2D 및 CVD를 연결한 연구가 이미 존재한다.

### C. Sleep → plasma proteomics → MR

2026년 UK Biobank 43,709명 연구에서 7개 sleep traits와 sleep health score, 935 proteins, MR 및 colocalization까지 수행되었다.

### D. Grip strength → proteomics → healthy aging

2025년 UK Biobank에서 grip strength, healthy aging 및 proteomics를 연결한 연구가 발표되었다.

### E. 단순 proteomic biological age

대규모 UKB 기반 proteomic age 연구가 빠르게 증가하고 있다.

---

## 3.2 본 연구에서 주장할 차별점

본 연구의 novelty는 개별 요소가 아니라 **조합**에서 나온다.

### Novelty 1
단일 질환이 아니라 **multi-domain metabolic health 유지 능력**을 분석한다.

### Novelty 2
cross-sectional metabolic health가 아니라 **20-year longitudinal deterioration/resilience**를 검증한다.

### Novelty 3
EUR pQTL discovery → EAS proteogenomic replication → Korean individual-level longitudinal validation을 연결한다.

### Novelty 4
candidate protein을 단순 biomarker가 아니라 **genetically anchored candidate**로 제한한다.

### Novelty 5
physical activity를 이용해 **precision wellness / modifiable lifestyle interaction**으로 확장한다.

---

# 4. 전체 연구 구조

```text
                    ┌──────────────────────────────┐
                    │      PROTEOMIC DISCOVERY     │
                    │ UKB-PPP / deCODE cis-pQTLs  │
                    └───────────────┬──────────────┘
                                    │
                                    ▼
                 ┌──────────────────────────────────┐
                 │     PROTEOME-WIDE MR SCREEN      │
                 │                                  │
                 │ Adiposity │ Glycemia │ Lipids   │
                 │ BP        │ T2D      │ +/- MetS │
                 └────────────────┬─────────────────┘
                                  │
                                  ▼
                    ┌──────────────────────────┐
                    │     COLOCALIZATION       │
                    │ shared causal variant?   │
                    └────────────┬─────────────┘
                                 │
                                 ▼
                    Candidate proteins
                         ~5–30 proteins
                                 │
                 ┌───────────────┴────────────────┐
                 │                                │
                 ▼                                ▼
       EAST ASIAN REPLICATION             FUNCTIONAL CHECK
     CKB pQTL + KoGES/BBJ/TWB          GTEx / scRNA / tissue
                 │                                │
                 └───────────────┬────────────────┘
                                 ▼
                     KOGES LONGITUDINAL
                        n ≈ 5k–8k*
                                 │
                                 ▼
                Stable health vs deterioration
                  + metabolic burden slope
                  + trajectory subclasses
                                 │
                                 ▼
                  cis-pQTL protein genetic score
                                 │
                                 ▼
                Lifestyle interaction / wellness
                    Physical activity primary
```

\* 실제 분석 N은 변수별 missingness 및 genotype availability audit 이후 확정한다.

---

# 5. 사용 가능한 핵심 데이터 자원

# 5.1 UK Biobank Pharma Proteomics Project (UKB-PPP)

## 역할
European ancestry 중심의 discovery pQTL source.

## 주요 특성

- 전체 proteomics participants: **54,219**
- Olink Explore 3072
- unique proteins: **2,923**
- primary genetic associations: **14,287**
- discovery pQTL analysis는 주로 European ancestry subset에서 수행
- non-European ancestry-specific pQTL 분석도 존재하지만 EAS sample size는 제한적

## 장점

- 대규모
- Olink 기반
- 비교적 robust한 pQTL mapping
- 다수 downstream genetic analyses에서 검증됨

## 주의점

- affinity proteomics이므로 실제 protein concentration과 assay binding signal이 완전히 동일한 개념은 아니다.
- missense variant가 epitope binding에 영향을 줄 수 있다.
- trans-pQTL은 pleiotropy risk가 크다.

## 본 연구 원칙

**Primary MR instrument는 cis-pQTL을 우선 사용한다.**

---

# 5.2 deCODE plasma proteomics

## 주요 특성

- Icelandic participants: **35,559**
- SomaScan
- aptamers: **4,907**
- reported plasma pQTL associations: **18,084 primary associations**

## 역할

1. UKB-PPP signal의 platform-independent replication
2. UKB-PPP에 없는 proteins 확장
3. Olink ↔ SomaScan cross-platform consistency 확인

## 주의점

동일 gene/protein이라도 Olink assay와 SomaScan aptamer가 측정하는 molecular feature가 완전히 동일하다고 가정하면 안 된다.

따라서 cross-platform replication은 다음 기준으로 본다.

- gene/protein mapping 일치
- effect direction 일치
- cis locus 일치
- colocalization 여부
- known epitope-related variant 여부

---

# 5.3 China Kadoorie Biobank East Asian pQTL resource

## 2026 preprint 기준

- Chinese participants: **3,965**
- SomaScan v4.1
- plasma proteins: **7,289**
- pQTLs: **3,212**
- cis-pQTL protein: **1,092**
- EAS phenotypes: 225 unique phenotypes
- sources:
  - BBJ
  - KoGES
  - Taiwan Biobank

## 중요한 특징

해당 연구는 pQTL과 East Asian trait GWAS를 통합해 다음을 수행했다.

- PheWAS
- colocalization
- MR
- cross-biobank comparison

reported high-confidence protein–phenotype pairs 중 다수가 cardiometabolic network에 위치한다.

## Data availability

논문 관련 summary results가 Zenodo에 공개되어 있다.

- `CKB_SomaScan_MR_coloc.gz`
- Supplementary Tables
- analysis code 일부

## 본 연구에서의 위치

**독립적인 primary discovery source가 아니라 EAS supportive / replication evidence로 사용**하는 것을 우선한다.

이유:

1. pQTL N이 UKB/deCODE보다 작음
2. 2026-09 현재 preprint
3. 일부 downstream EAS trait과 본 연구에서 사용할 KoGES GWAS가 overlap할 수 있음

따라서 EAS replication에서 sample overlap을 반드시 기록한다.

---

# 5.4 KoGES public GWAS summary statistics

KoGES GWAS resource:

- participants: **72,298**
- KoreanChip genotype + imputation
- variants: 약 **8.06 million**
- phenotypes: **76**
- genome build: **GRCh37**
- 공개 summary statistics 이용 가능

본 연구에서 우선 사용할 수 있는 trait 예시:

| Trait | Approximate public KoGES GWAS N |
|---|---:|
| BMI | 72,282 |
| SBP | 72,265 |
| HDL-C | 72,297 |
| TG | 72,195 |
| Diabetes | 5,083 cases / 67,127 controls |
| Hypertension | 14,834 cases / 57,395 controls |
| Fasting blood sugar | available |
| DBP | available |
| Waist circumference | available |
| LDL / total cholesterol | available |

### 역할

- EAS outcome replication
- CKB pQTL × Korean outcome MR
- EUR → EAS effect-direction comparison
- candidate protein prioritization

---

# 5.5 KoGES Ansan–Ansung individual-level longitudinal cohort

## 참여자 추적

| Survey | Years | Subjects | Genotyped subjects* |
|---|---|---:|---:|
| Baseline | 2001–2002 | 10,030 | 8,198 |
| 1st | 2003–2004 | 8,603 | 7,184 |
| 2nd | 2005–2006 | 7,515 | 6,358 |
| 3rd | 2007–2008 | 6,688 | 5,693 |
| 4th | 2009–2010 | 6,665 | 5,727 |
| 5th | 2011–2012 | 6,238 | 5,331 |
| 6th | 2013–2014 | 5,906 | 5,025 |
| 7th | 2015–2016 | 6,318 | 5,345 |
| 8th | 2017–2018 | 6,157 | 5,212 |
| 9th | 2019–2020 | 5,854 | 4,949 |
| 10th | 2021–2022 | 5,511 | 4,658 |

\* OPEN KoGES public status table 기준. 실제 연구용 genotype release와 phenotype release의 participant intersection은 별도 audit 필요.

## 장점

- 약 20년 follow-up
- 2년 단위 repeated measures
- anthropometry
- blood pressure
- fasting glucose
- lipid measurements
- smoking/alcohol
- physical activity
- diet
- chronic disease status
- genotype

## 가장 중요한 점

본 연구에서 KoGES는 단순 replication cohort가 아니다.

> **“metabolic resilience”라는 longitudinal phenotype을 실제로 생성하는 핵심 cohort**

이다.

---

# 6. Outcome phenotype 설계

Resilience라는 새로운 phenotype을 사용하므로, phenotype definition이 논문의 가장 중요한 부분이다.

따라서 하나의 임의적인 score에 의존하지 않고 **3-level phenotype strategy**를 사용한다.

---

# 6.1 Primary phenotype: clinically interpretable metabolic transition

한국 성인의 metabolic syndrome criterion을 활용한다.

각 wave에서 다음 5개 component를 계산한다.

1. Abdominal obesity
   - Male: WC ≥ 90 cm
   - Female: WC ≥ 85 cm

2. Elevated BP
   - SBP ≥ 130 mmHg 또는 DBP ≥ 85 mmHg
   - 또는 antihypertensive medication

3. Elevated fasting glucose
   - fasting glucose ≥ 100 mg/dL
   - 또는 diabetes medication

4. Elevated triglycerides
   - TG ≥ 150 mg/dL
   - medication 정보를 사용할 수 있는 경우 dyslipidemia medication 반영

5. Low HDL-C
   - Male: HDL < 40 mg/dL
   - Female: HDL < 50 mg/dL
   - medication 정의는 codebook 확인 후 harmonization

MetS:

\[
MetS_t = I(ComponentCount_t \geq 3)
\]

---

## Baseline healthy definition

Primary:

\[
ComponentCount_{baseline} \leq 1
\]

Sensitivity:

- strict healthy: 0 components
- permissive healthy: <3 components

---

## Stable metabolic resilience

초기 operational definition:

```text
Baseline components <= 1
AND
never develops MetS during eligible follow-up
AND
>= predefined minimum number of valid visits
```

더 엄격한 sensitivity definition:

```text
Baseline components <= 1
AND
components <= 1 in >=70% of observed follow-ups
AND
never develops T2D
```

---

## Metabolic deterioration

Primary deterioration:

```text
Baseline components <= 1
→ incident MetS (>=3 components)
```

Secondary deterioration:

```text
Baseline components <=1
→ incident T2D
OR
persistent >=2 component increase
```

---

## Time-to-event approach

Healthy participants at baseline을 대상으로:

\[
T = time\ to\ first\ incident\ MetS
\]

분석:

- Cox proportional hazards
- interval censoring sensitivity analysis
- competing mortality issue는 data availability에 따라 고려

---

# 6.2 Secondary phenotype: Continuous Metabolic Burden Index (MBI)

Binary metabolic syndrome은 threshold 주변의 정보를 잃는다.

따라서 repeated continuous variables를 이용한 burden score를 만든다.

초기 설계:

\[
MBI_{it} =
\frac{
Z(WC_{it})
+ Z(SBP_{it})
+ Z(FPG_{it})
+ Z(\log TG_{it})
- Z(HDL_{it})
}{5}
\]

성별 차이가 큰 variable은 sex-stratified standardization을 우선 고려한다.

### 후보 variable

- waist circumference
- SBP
- fasting glucose
- log(TG)
- HDL-C inverse

### Optional expansion

- DBP
- HbA1c
- BMI
- HOMA-IR

단, 모든 wave에서 측정되지 않는 variable을 primary score에 포함하면 longitudinal completeness가 크게 감소할 수 있다.

**따라서 Stage 0 missingness audit 이후 최종 variable set을 결정한다.**

---

# 6.3 MBI slope

Linear mixed model:

\[
MBI_{it}
=
\beta_0
+
b_{0i}
+
(\beta_1+b_{1i})Time_{it}
+
\beta_2Age_i
+
\beta_3Sex_i
+\cdots+
\epsilon_{it}
\]

개인별 random slope:

\[
Slope_i = \beta_1 + b_{1i}
\]

해석:

- negative / near-zero slope: metabolic stability
- large positive slope: metabolic deterioration

---

# 6.4 Residualized resilience score

단순 slope는 baseline 상태에 영향을 많이 받는다.

따라서 exploratory analysis에서 expected deterioration을 모델링한다.

\[
ExpectedSlope_i =
f(Age, Sex, BaselineBMI, BaselineMBI, Smoking, ...)
\]

\[
ResilienceResidual_i
=
ObservedSlope_i - ExpectedSlope_i
\]

- 음의 residual: 예상보다 건강하게 유지
- 양의 residual: 예상보다 빠른 악화

주의:

이 변수는 모델 기반 phenotype이므로 **primary endpoint로 두지 않는다.**

---

# 6.5 Exploratory phenotype: trajectory classes

후보 방법:

- latent class mixed model
- group-based trajectory modeling
- growth mixture model
- k-means on standardized longitudinal features는 sensitivity 용도

예상 class:

```text
Class 1: Stable low burden
Class 2: Gradual increase
Class 3: Rapid deterioration
Class 4: Persistent high burden
```

권장:

- 2–5 class 비교
- BIC/AIC
- entropy
- minimum class size
- clinical interpretability

중요:

trajectory class 수를 데이터에 맞춰 무제한 탐색하면 overfitting이 발생한다.

따라서 **Primary analysis는 pre-defined transition / MBI slope**, trajectory class는 exploratory로 둔다.

---

# 7. MHO-specific sub-study

Metabolic resilience 개념을 가장 직관적으로 보여주는 secondary study.

## Korean obesity definition

\[
BMI \geq 25\ kg/m^2
\]

## MHO

Primary exploratory definition:

```text
BMI >=25
AND
MetS components <=1
```

Sensitivity:

- strict MHO: obesity + zero metabolic abnormalities
- conventional MHO: obesity + <3 MetS components

## MUO

```text
BMI >=25
AND
MetS >=3 components
```

## Longitudinal transition

\[
MHO \rightarrow MHO
\]

versus

\[
MHO \rightarrow MUO
\]

질문:

> 어떤 genetically predicted plasma proteins가 비만 상태에서도 대사적으로 건강한 상태를 유지하는 것과 관련되는가?

이 분석은 main analysis와 별도의 secondary manuscript 또는 major sensitivity analysis로 확장할 수 있다.

---

# 8. Stage 1: Proteome-wide causal screening

# 8.1 Exposure

우선순위:

1. UKB-PPP cis-pQTL
2. deCODE cis-pQTL
3. intersection / cross-platform replication

Primary instrument 기준 예시:

```text
cis window: ±1 Mb from coding gene
pQTL P < genome-wide / study-specific threshold
F-statistic > 10
LD clumping appropriate to ancestry
```

trans-pQTL은 primary causal evidence에서 제외하거나 별도 sensitivity로 둔다.

---

# 8.2 Outcome domains

Protein을 한 개의 “MetS GWAS”에만 넣지 않는다.

최소 4개 biological domains를 사용한다.

## Domain A — Adiposity

- BMI
- Waist circumference
- WHR / WHRadjBMI 가능 시

## Domain B — Glycemia

- fasting glucose
- HbA1c
- fasting insulin optional
- T2D

## Domain C — Lipids

- triglycerides
- HDL-C
- LDL-C optional

## Domain D — Blood pressure

- SBP
- DBP
- hypertension optional

---

# 8.3 추천 discovery GWAS resources

## T2D

T2D Global Genetics Initiative / DIAMANTE 계열.

대규모 multi-ancestry study:

- total N >2.5 million
- T2D cases >428k

장점:

- power 매우 높음
- ancestry diversity

주의:

UKB overlap 가능성을 outcome dataset별로 기록한다.

---

## Lipids

Global Lipids Genetics Consortium.

대표 large-scale dataset:

- total N ≈1.65 million
- EAS ≈146k 포함

주요 traits:

- TG
- HDL-C
- LDL-C
- total cholesterol
- non-HDL-C

---

## Glycemic traits

MAGIC.

public summary statistics:

- fasting glucose
- fasting insulin
- HbA1c
- 2h glucose

multi-ancestry / ancestry-specific 결과 사용 가능.

---

## Blood pressure

2024 large European BP GWAS:

- up to **1,028,980**
- SBP
- DBP
- pulse pressure
- GWAS Catalog accession:
  - SBP: GCST90310294
  - DBP: GCST90310295
  - PP: GCST90310296

---

## BMI / anthropometry

GIANT / public anthropometric GWAS 및 필요한 경우 Pan-UKB/other large-scale resources.

주의:

최신 5.1M BMI study의 모든 full GWAS statistics가 동일한 수준으로 공개되어 있는지 실제 download 단계에서 확인해야 한다.

Stage 0에서는 **완전 공개 + reproducible download가 가능한 dataset을 우선**한다.

---

# 9. Candidate protein prioritization rule

단순히 MR P-value 순으로 protein을 고르면 위험하다.

다음 evidence hierarchy를 사용한다.

---

## Tier 1: High-confidence metabolic resilience candidate

다음 조건을 목표로 한다.

```text
1. cis-pQTL MR significant after FDR correction
2. coloc PP.H4 >=0.80
3. >=2 metabolic domains에서 유리한 방향의 evidence
4. 다른 주요 domain에서 강한 반대 방향 evidence 없음
5. Steiger direction consistent
6. EAS에서 directionally concordant evidence 존재
```

---

## Tier 2

```text
MR + coloc strong
but only one domain
OR
multi-domain evidence exists but EAS evidence weak
```

---

## Tier 3

```text
MR significant
but colocalization weak/absent
or platform inconsistency
```

Tier 3는 discovery list에는 남겨도 KoGES individual-level validation 우선순위에서는 낮춘다.

---

# 10. Multi-domain evidence score

formal statistical significance test와 별도로 **prioritization score**를 만든다.

예:

| Domain | Favorable supported | Neutral | Harmful supported |
|---|---:|---:|---:|
| Adiposity | +1 | 0 | -1 |
| Glycemia | +1 | 0 | -1 |
| Lipids | +1 | 0 | -1 |
| BP | +1 | 0 | -1 |
| T2D | +1 | 0 | -1 |

단, correlated outcomes 때문에 단순 합계를 formal P-value처럼 해석하면 안 된다.

따라서 이 score는:

> **candidate ranking / evidence visualization용**

으로만 사용한다.

---

# 11. MR analysis specification

# 11.1 Primary MR

Instrument가 1개인 경우:

\[
\hat{\beta}_{MR}
=
\frac{\beta_{SNP \rightarrow outcome}}
{\beta_{SNP \rightarrow protein}}
\]

Wald ratio.

여러 independent cis instruments가 있는 경우:

- IVW
- correlated instruments 사용 시 LD-aware method 고려

---

# 11.2 Sensitivity analyses

가능한 경우:

- weighted median
- MR-Egger
- MR-PRESSO
- GSMR
- Steiger directionality
- leave-one-out
- heterogeneity Q

하지만 cis-pQTL candidate에서 instrument 수가 적을 수 있으므로 모든 sensitivity method가 항상 가능한 것은 아니다.

---

# 11.3 Major MR risks

### Horizontal pleiotropy
특히 trans-pQTL에서 문제.

### LD confounding
protein pQTL variant와 outcome causal variant가 서로 다른데 LD만 되어 있을 수 있다.

### Reverse causality
disease variant가 protein level을 변화시키는 경우.

### Sample overlap
UKB pQTL과 UKB-derived outcome GWAS가 겹칠 수 있다.

### Winner's curse
discovery pQTL effect가 과대평가될 수 있다.

---

# 12. Colocalization

MR significance만으로는 부족하다.

같은 locus에서 protein과 metabolic trait이 **동일 causal variant를 공유하는지** 평가한다.

---

## Primary method

`coloc`

예상 기준:

\[
PP.H4 > 0.80
\]

강한 evidence:

\[
PP.H4 > 0.90
\]

추가 criterion:

\[
PP.H4 / PP.H3 > 5
\]

CKB 2026 atlas에서도 유사한 조건을 사용했다.

---

## Multiple causal variants

단일 causal variant assumption이 깨지는 locus에서는:

- SuSiE-based coloc
- conditional analysis

를 고려한다.

---

# 13. Stage 2: Cross-ancestry replication

# 13.1 EUR discovery

```text
UKB-PPP/deCODE cis-pQTL
          ↓
EUR or multi-ancestry metabolic GWAS
          ↓
MR + coloc
```

---

# 13.2 EAS replication

```text
CKB cis-pQTL
      ×
KoGES / BBJ / TWB metabolic GWAS
```

가능한 경우 ancestry-matched EAS LD reference 사용.

---

# 13.3 Cross-ancestry comparison items

protein마다 다음을 기록한다.

| Field | Description |
|---|---|
| gene | protein coding gene |
| protein | standardized protein name |
| pQTL SNP EUR | lead EUR cis-pQTL |
| pQTL SNP EAS | lead EAS cis-pQTL |
| EAF EUR | effect allele frequency |
| EAF EAS | effect allele frequency |
| beta protein EUR | pQTL effect |
| beta protein EAS | pQTL effect |
| MR beta EUR | outcome causal estimate |
| MR beta EAS | outcome causal estimate |
| direction concordance | yes/no |
| coloc EUR | PP.H4 |
| coloc EAS | PP.H4 |
| platform | Olink/SomaScan |
| status | replicated / suggestive / discordant |

---

# 13.4 Important interpretation

EAS에서 non-significant라고 해서 EUR finding이 자동으로 반박되는 것은 아니다.

구분해야 한다.

### A. Concordant + significant
strong replication.

### B. Concordant + non-significant
likely underpowered.

### C. Opposite direction
potential heterogeneity.

### D. EUR pQTL absent / rare in EAS
not testable.

---

# 14. Stage 3: KoGES longitudinal individual-level validation

이 stage가 연구의 핵심이다.

# 14.1 대상자 selection

초기 inclusion:

```text
KoGES Ansan/Ansung
baseline genotype available
baseline metabolic variables available
>= predefined follow-up visits
```

Primary recommended minimum:

- baseline + ≥2 follow-up measurements

Sensitivity:

- baseline + ≥4 follow-ups
- ≥10 years follow-up
- complete long-term group

---

# 14.2 예상 N

확정값이 아니다.

가능한 upper bound:

- baseline total: 10,030
- baseline genotype: 8,198
- 10th follow-up: 5,511
- genotype among 10th participants: 4,658

따라서 실제 longitudinal genetic analysis는 대략:

\[
N \approx 5,000-8,000
\]

범위에서 결정될 가능성이 높다.

그러나 이는 **Stage 0 audit 이전 예상치**이다.

실제 N은 다음에 의해 감소한다.

- missing fasting lab
- medication data
- genotype QC
- phenotype exclusion
- visit-count criterion
- baseline healthy restriction

---

# 14.3 Candidate protein genetic score

protein k에 대해 independent cis-pQTL이 여러 개라면:

\[
PGS_{ik}
=
\sum_j
G_{ij}\beta_{jk}^{pQTL}
\]

이때:

- `Gij`: dosage
- `β pQTL`: external pQTL weight
- allele harmonization 필수

한 개의 strong cis-pQTL만 있으면 single-variant genetic proxy로 분석 가능.

---

# 14.4 Primary individual-level models

## Model A — Incident MetS

\[
h(t)
=
h_0(t)
\exp(
\beta_1 ProteinGRS
+
\beta_2 Age
+
\beta_3 Sex
+
PCs
)
\]

Primary genetic model에서는 genetic instrument 특성상 lifestyle covariate adjustment를 과도하게 넣지 않는다.

추가 observational/context model에서:

- smoking
- alcohol
- physical activity
- baseline BMI

등을 추가한다.

---

## Model B — Repeated MBI

\[
MBI_{it}
=
\beta_0+
\beta_1Time+
\beta_2ProteinGRS+
\beta_3ProteinGRS\times Time+
b_{0i}+b_{1i}Time+\epsilon
\]

핵심:

\[
\beta_3
\]

즉 genetically predicted protein level이 metabolic deterioration slope와 관련되는지를 본다.

---

## Model C — trajectory class

Multinomial logistic regression:

\[
TrajectoryClass \sim ProteinGRS + Age + Sex + PCs
\]

---

# 15. Lifestyle interaction / Precision Wellness layer

이 부분이 Digital Health identity를 강화한다.

# 15.1 Primary modifier

**Physical activity**

이유:

1. wellness에서 직접 조절 가능한 행동이다.
2. KoGES에서 관련 questionnaire 변수가 존재한다.
3. exercise와 proteomic aging 간 연결을 보여주는 선행 evidence가 있다.
4. 향후 wearable-based intervention 연구로 확장할 수 있다.

---

# 15.2 Interaction model

\[
MBI_{it}
=
...
+
\beta_4 ProteinGRS_i \times PA_{it}
+
\beta_5 ProteinGRS_i \times PA_{it} \times Time
\]

혹은 simplified model:

\[
MetabolicDeterioration
\sim
ProteinGRS
+
PA
+
ProteinGRS\times PA
\]

---

# 15.3 PA variable의 우선순위

Stage 0 codebook audit에서 확인:

1. MET-min/week 계산 가능 여부
2. moderate/vigorous activity frequency
3. regular exercise indicator
4. sedentary time
5. repeated measurement availability

가장 반복성이 높은 변수를 primary로 선택한다.

---

# 15.4 Secondary lifestyle variables

- diet score
- smoking
- alcohol
- sleep duration
- sleep quality/insomnia
- sedentary behavior

단, 모두 main interaction으로 넣으면 multiple testing burden이 커진다.

따라서:

```text
Primary interaction: Physical activity
Secondary: Sleep
Exploratory: Diet / alcohol / smoking
```

로 제한한다.

---

# 16. Functional annotation

candidate protein을 통계적으로만 끝내지 않는다.

# 16.1 Tissue expression

- GTEx
- Human Protein Atlas
- tissue-enriched expression

우선 tissue:

- liver
- adipose
- skeletal muscle
- pancreas
- vascular/endothelial tissue
- hypothalamus when relevant

---

# 16.2 Single-cell annotation

public scRNA/snRNA datasets를 이용하여 candidate gene expression 위치를 확인한다.

질문:

> candidate protein의 coding gene은 어느 cell type에서 주로 발현되는가?

예:

- hepatocyte
- adipocyte
- myocyte
- beta cell
- endothelial cell
- macrophage
- immune cells

---

# 16.3 Functional priority

다음 evidence를 함께 table화한다.

```text
MR
+ coloc
+ EAS replication
+ KoGES longitudinal
+ tissue expression
+ single-cell expression
+ druggability
```

---

# 17. 통계적 multiple testing 전략

Proteome-wide analysis에서 Bonferroni만 사용하면 매우 엄격할 수 있다.

추천:

## Discovery
FDR q < 0.05

## Strong evidence
MR FDR <0.05 + coloc PP.H4 >0.8

## Cross-domain prioritization
여러 trait에서 반복 확인.

## KoGES validation
candidate protein 수를 사전에 5–30개로 축소.

KoGES에서 다시 3,000–7,000 proteins를 screening하면 N이 부족하고 본 연구의 구조도 흐려진다.

---

# 18. Sample overlap 관리

반드시 `DATASET_REGISTRY.tsv`를 만든다.

예:

| dataset | ancestry | N | phenotype | UKB included | KoGES included | role |
|---|---|---:|---|---|---|---|
| UKB-PPP | EUR-majority | 54,219 | proteins | yes | no | pQTL |
| deCODE | EUR | 35,559 | proteins | no | no | pQTL |
| T2DGGI | multi | >2.5M | T2D | possible | possible | outcome |
| GLGC | multi | 1.65M | lipids | possible | possible | outcome |
| KoGES GWAS | EAS-Korean | 72,298 | 76 traits | no | yes | EAS outcome |
| CKB pQTL | EAS-Chinese | 3,965 | proteins | no | no | EAS pQTL |

같은 cohort가 exposure와 outcome 양쪽에 포함되면 weak instrument bias 방향에 영향을 줄 수 있으므로 기록한다.

---

# 19. Genome build strategy

현재 주요 resources가 GRCh37와 GRCh38이 혼재할 수 있다.

## Master analysis build

**GRCh37를 우선 고려**

이유:

- KoGES public GWAS PheWeb: GRCh37
- 다수 legacy GWAS resources와 호환
- CKB EAS analysis도 GRCh37 기반 summary statistics를 사용

단, 특정 UKB pQTL source가 GRCh38 중심이면 liftover를 명시적으로 수행한다.

필수 columns:

```text
chr
pos
rsid
effect_allele
other_allele
eaf
beta
se
p
N
build
```

---

# 20. Harmonization QC

모든 exposure/outcome pair에서:

```text
1. chromosome normalization
2. build check
3. rsID consistency
4. effect allele alignment
5. strand check
6. palindromic SNP handling
7. allele frequency comparison
8. duplicate SNP removal
9. missing variant proxy lookup
```

Palindromic SNP:

- A/T
- C/G

MAF ≈0.5이면 ambiguity 때문에 제거를 우선한다.

---

# 21. Proteomic assay artifact 관리

Affinity proteomics pQTL 연구의 중요한 문제.

missense variant가 protein abundance가 아니라 **binding affinity**를 변화시킬 수 있다.

따라서 candidate마다:

```text
coding/missense variant?
epitope-binding concern?
Olink replicated?
SomaScan replicated?
cis signal shared across platforms?
```

를 기록한다.

cross-platform replicated cis signal은 높은 우선순위를 준다.

---

# 22. Negative controls / robustness

가능하면 다음을 포함한다.

## Negative control 1
protein과 생물학적으로 무관한 phenotype에서 광범위한 pleiotropy가 있는지 확인.

## Negative control 2
cis vs trans 결과 비교.

## Negative control 3
EUR vs EAS heterogeneity.

## Leave-one-domain-out

예:

```text
candidate score 계산 시 lipids domain 제외
candidate score 계산 시 glycemia 제외
...
```

특정 domain 하나가 전체 signal을 지배하는지 확인한다.

---

# 23. Stage 0 — 실제 분석 전 필수 audit

가장 먼저 할 작업이다.

## 23.1 KoGES codebook audit

찾아야 할 variables:

### ID / visit

- participant ID
- survey wave
- survey date/year
- age
- sex

### Anthropometry

- height
- weight
- BMI
- waist circumference
- hip circumference

### BP

- SBP repeated measurements
- DBP repeated measurements
- antihypertensive medication

### Glycemia

- fasting glucose
- HbA1c
- fasting insulin
- diabetes diagnosis
- diabetes medication

### Lipids

- TG
- HDL-C
- LDL-C
- total cholesterol
- lipid-lowering medication

### Lifestyle

- physical activity
- exercise frequency
- MET if available
- sedentary time
- smoking
- alcohol
- sleep
- diet/FFQ

### Genetics

- sample ID
- genotype ID
- chip version
- genotype build
- PCs
- relatedness
- imputation information

---

# 23.2 Missingness matrix

output:

`KOGES_VARIABLE_AVAILABILITY_BY_WAVE.tsv`

예:

| variable | BL | F1 | F2 | ... | F10 |
|---|---:|---:|---:|---|---:|
| waist | 99% | ... | ... | | |
| SBP | 99% | | | | |
| FPG | | | | | |
| TG | | | | | |
| HDL | | | | | |
| PA | | | | | |
| sleep | | | | | |

---

# 23.3 Participant retention

output:

`KOGES_LONGITUDINAL_ELIGIBILITY.tsv`

columns:

```text
participant_id
n_total_visits
n_valid_metabolic_visits
baseline_healthy
incident_mets
incident_t2d
followup_years
genotype_available
```

---

# 23.4 Feasibility gate

다음 조건을 만족하면 full study 진행.

### GO

- ≥4 metabolic variables가 대부분 wave에 존재
- longitudinal eligible N ≥4,000
- baseline healthy genetic N이 충분
- physical activity가 ≥3개 wave에서 usable
- genotype/pQTL allele harmonization 가능

### MODIFY

- PA repeated data 부족
  → baseline PA stratification으로 변경

- medication variable 불완전
  → lab-only continuous MBI를 primary secondary phenotype으로 강화

### NO-GO / alternative

- genotype-longitudinal intersection이 너무 작음
  → KoGES longitudinal observational validation + EAS summary-MR로 역할 분리

---

# 24. 분석 directory 제안

기존 연구 서버 구조와 분리하기 위해:

```text
/srv/is-analysis/
└── results/
    └── metabolic_resilience/
        ├── stage0_feasibility/
        │   ├── koges_codebook/
        │   ├── availability/
        │   ├── participant_flow/
        │   └── reports/
        │
        ├── stage1_pqtl/
        │   ├── ukb_ppp/
        │   ├── decode/
        │   └── harmonized/
        │
        ├── stage2_gwas/
        │   ├── adiposity/
        │   ├── glycemia/
        │   ├── lipids/
        │   ├── blood_pressure/
        │   └── t2d/
        │
        ├── stage3_mr_coloc/
        │   ├── mr/
        │   ├── coloc/
        │   ├── sensitivity/
        │   └── candidates/
        │
        ├── stage4_eas/
        │   ├── ckb_pqtl/
        │   ├── koges_gwas/
        │   ├── bbj/
        │   └── replication/
        │
        ├── stage5_koges_longitudinal/
        │   ├── phenotype/
        │   ├── genotype/
        │   ├── protein_scores/
        │   ├── models/
        │   └── trajectories/
        │
        ├── stage6_lifestyle/
        │
        ├── stage7_functional/
        │   ├── tissue/
        │   ├── single_cell/
        │   └── druggability/
        │
        ├── figures/
        ├── tables/
        ├── manuscript/
        └── audit/
```

---

# 25. 권장 scripts

```text
scripts/
├── audit_metabolic_resilience_inputs.py
├── build_dataset_registry.py
├── audit_koges_variables.py
├── build_koges_longitudinal_panel.py
├── derive_mets_components.py
├── derive_metabolic_burden_index.py
├── fit_metabolic_trajectories.R
├── harmonize_pqtl_gwas.py
├── run_protein_mr.R
├── run_coloc.R
├── prioritize_metabolic_proteins.py
├── replicate_eas.py
├── build_protein_grs.py
├── run_koges_longitudinal_models.R
├── run_lifestyle_interactions.R
└── build_final_evidence_matrix.py
```

---

# 26. Master candidate table

최종적으로 모든 분석은 하나의 table로 합친다.

`METABOLIC_RESILIENCE_CANDIDATES.tsv`

예상 columns:

```text
protein
gene
uniprot
platform_discovery
cis_pqtl
effect_allele
other_allele
eaf
F_stat

mr_bmi_beta
mr_bmi_p
coloc_bmi_pph4

mr_wc_beta
mr_wc_p
coloc_wc_pph4

mr_fpg_beta
mr_fpg_p
coloc_fpg_pph4

mr_hba1c_beta
mr_hba1c_p
coloc_hba1c_pph4

mr_tg_beta
mr_tg_p
coloc_tg_pph4

mr_hdl_beta
mr_hdl_p
coloc_hdl_pph4

mr_sbp_beta
mr_sbp_p
coloc_sbp_pph4

mr_t2d_beta
mr_t2d_p
coloc_t2d_pph4

domains_supported
discordant_domains

ckb_pqtl_available
eas_direction_concordant
koges_gwas_support
bbj_support

koges_longitudinal_beta
koges_longitudinal_p
protein_grs_x_time

pa_interaction_beta
pa_interaction_p

liver_expression
adipose_expression
muscle_expression
pancreas_expression

single_cell_annotation
druggable
known_drug
assay_artifact_flag

evidence_tier
```

---

# 27. Figure plan

## Figure 1 — Study design

EUR discovery → EAS replication → KoGES longitudinal validation.

## Figure 2 — Proteome × metabolic trait evidence map

Heatmap:

```text
rows = proteins
columns = BMI/WC/FPG/HbA1c/TG/HDL/SBP/T2D
color/value = MR direction/effect
symbol = coloc support
```

## Figure 3 — Cross-ancestry replication

EUR MR beta vs EAS MR beta scatter.

## Figure 4 — KoGES trajectories

MBI trajectory classes 또는 predicted longitudinal curves.

## Figure 5 — Protein GRS longitudinal effect

candidate별 GRS quartile에 따른 metabolic burden trajectory.

## Figure 6 — Lifestyle interaction

PA high vs low에서 protein GRS effect.

## Figure 7 — Functional map

tissue / cell type / pathway network.

---

# 28. Primary tables

## Table 1
Dataset characteristics.

## Table 2
Proteome-wide MR + coloc candidate proteins.

## Table 3
Cross-ancestry replication.

## Table 4
KoGES longitudinal validation.

## Table 5
Lifestyle interaction.

## Supplementary

- all proteins
- all traits
- all MR methods
- sensitivity
- phenotype definitions
- medication adjustments
- variable missingness
- participant flow

---

# 29. Manuscript aims

## Aim 1

Identify plasma proteins with genetically supported effects across multiple cardiometabolic domains.

### Hypothesis
A subset of proteins exhibits concordant causal effects across adiposity, glycemia, lipids, and blood pressure.

---

## Aim 2

Determine whether proteogenomic signals are conserved across European and East Asian populations.

### Hypothesis
Core metabolic proteins are cross-ancestry conserved, while some signals show ancestry-specific effects due to LD and allele-frequency differences.

---

## Aim 3

Test whether genetically predicted levels of prioritized proteins are associated with long-term metabolic resilience in Koreans.

### Hypothesis
Protective protein genetic proxies are associated with lower longitudinal metabolic burden and reduced risk of incident MetS.

---

## Aim 4

Evaluate whether physical activity modifies genetically anchored metabolic resilience.

### Hypothesis
The longitudinal effect of candidate proteins differs according to physical activity level.

---

# 30. Proposed manuscript title options

### Main

**Proteogenomic determinants of metabolic resilience: cross-ancestry causal inference and longitudinal validation in a Korean population**

### Alternative 1

**Genetically anchored plasma proteins underlying long-term metabolic resilience across European and East Asian populations**

### Alternative 2

**From plasma proteomics to healthy metabolic aging: cross-ancestry causal inference and 20-year validation in Koreans**

### Alternative 3

**Molecular determinants of metabolic resilience: proteogenomic discovery and longitudinal validation in the Korean Genome and Epidemiology Study**

---

# 31. 예상 abstract logic

## Background
Traditional cardiometabolic studies focus on disease onset, whereas determinants enabling maintenance of favorable metabolic health during aging remain less well characterized.

## Methods
Plasma cis-pQTLs from large-scale proteogenomic cohorts are integrated with GWAS of adiposity, glycemia, lipid traits, blood pressure and T2D using MR and colocalization. Candidates are assessed in East Asian genetic resources and validated using repeated metabolic measurements in KoGES.

## Results
TBD.

## Conclusion
Genetically anchored circulating proteins may identify molecular mechanisms contributing to long-term metabolic resilience and provide targets for precision prevention.

---

# 32. Major limitations anticipated

## 32.1 Metabolic resilience definition is not standardized

대응:

- primary clinical transition
- continuous MBI
- sensitivity definitions
- trajectory analysis

여러 정의에서 일관성을 확인한다.

---

## 32.2 KoGES에 large-scale plasma proteomics가 없을 수 있음

대응:

실제 protein abundance validation이 아니라:

\[
cis-pQTL \rightarrow genetically\ predicted\ protein
\]

을 검증한다.

논문에서 이를 명확하게 기술한다.

---

## 32.3 KoGES N은 large GWAS보다 작음

대응:

KoGES는 discovery가 아니라 **candidate validation**에 사용한다.

후보를 Stage 1–2에서 5–30 proteins 수준으로 줄인다.

---

## 32.4 Attrition

20년 추적에서 healthy survivor bias가 발생할 수 있다.

대응 후보:

- participant characteristics by follow-up status
- inverse probability weighting
- mixed models using incomplete repeated measures
- complete-case sensitivity

---

## 32.5 Medication

BP, glucose, lipids가 medication에 의해 변한다.

대응:

- clinical component에서는 medication status 반영
- continuous phenotype에서는 medication-adjusted sensitivity analysis
- raw lab-only analysis도 병행

---

## 32.6 pQTL assay artifacts

대응:

- cis-first
- coloc
- missense flag
- cross-platform verification

---

## 32.7 Cross-ancestry power difference

CKB pQTL sample size가 상대적으로 작다.

대응:

EAS non-significance를 곧바로 heterogeneity로 해석하지 않는다.

direction, allele frequency, effect-size CI, pQTL availability를 함께 평가한다.

---

# 33. 석사 논문 범위 조절

모든 분석을 한 번에 끝내려 하면 과도하다.

## Minimum Viable Thesis

```text
1. UKB/deCODE cis-pQTL
2. metabolic trait MR
3. colocalization
4. candidate proteins
5. KoGES GWAS EAS replication
6. KoGES longitudinal phenotype
7. protein GRS → metabolic trajectory
```

여기까지만 성공해도 충분한 논문 구조다.

---

## Full thesis

추가:

```text
CKB pQTL replication
physical activity interaction
MHO sub-study
single-cell annotation
druggability
```

---

## Post-thesis expansion

```text
wearable physical activity / sleep
multi-omics
metabolomics
causal mediation
prospective precision wellness score
intervention cohort
```

---

# 34. CKD 연구와의 관계

CKD pipeline에서 재사용 가능한 모듈:

- pQTL ingestion
- allele harmonization
- instrument QC
- MR
- colocalization
- EAS replication
- candidate ranking
- single-cell annotation
- audit/report framework

새로 필요한 모듈:

- multi-domain metabolic outcome aggregation
- KoGES repeated-measures phenotype construction
- metabolic burden index
- trajectory modeling
- lifestyle interaction

따라서 CKD 연구 자산을 상당 부분 재사용할 수 있다.

---

# 35. 즉시 실행할 Stage 0 작업

## Step 0-1 — 연구 폴더 생성

```bash
mkdir -p /srv/is-analysis/results/metabolic_resilience/{stage0_feasibility,audit}
```

---

## Step 0-2 — dataset registry 생성

파일:

```text
/srv/is-analysis/results/metabolic_resilience/stage0_feasibility/DATASET_REGISTRY.tsv
```

최소 columns:

```text
dataset_id
dataset_name
role
ancestry
platform
phenotype
n_total
n_cases
n_controls
genome_build
download_url
accession
publication
sample_overlap_notes
status
```

---

## Step 0-3 — KoGES data inventory

현재 서버에 이미 확보된 KoGES 자료가 있다면 다음을 출력한다.

```bash
find /srv/is-analysis/data/koges -maxdepth 3 -type f \
  -printf '%p\t%s\n' \
  | sort \
  > /srv/is-analysis/results/metabolic_resilience/stage0_feasibility/KOGES_FILE_INVENTORY.tsv
```

---

## Step 0-4 — KoGES variable audit

목표 output:

```text
KOGES_VARIABLE_CANDIDATES.tsv
KOGES_VARIABLE_AVAILABILITY_BY_WAVE.tsv
KOGES_MEDICATION_VARIABLES.tsv
KOGES_LIFESTYLE_VARIABLES.tsv
```

---

## Step 0-5 — longitudinal sample flow

목표:

```text
N baseline
N genotype
N >=2 metabolic visits
N >=3 metabolic visits
N >=5 metabolic visits
N baseline healthy
N incident MetS
N stable healthy
N incident T2D
```

---

## Step 0-6 — GWAS download audit

확인할 outcome:

```text
BMI
WC
WHR
FPG
HbA1c
TG
HDL
SBP
DBP
T2D
```

각 trait마다:

```text
EUR dataset
EAS dataset
N
build
downloadable?
full summary stats?
sample overlap?
```

를 확정한다.

---

# 36. Stage 0 완료 기준

다음 파일들이 생성되면 Stage 1로 넘어간다.

```text
DATASET_REGISTRY.tsv
KOGES_FILE_INVENTORY.tsv
KOGES_VARIABLE_CANDIDATES.tsv
KOGES_VARIABLE_AVAILABILITY_BY_WAVE.tsv
KOGES_LONGITUDINAL_SAMPLE_FLOW.tsv
GWAS_RESOURCE_MATRIX.tsv
PQTLS_RESOURCE_MATRIX.tsv
STAGE0_GO_NO_GO.md
```

---

# 37. 첫 번째 분석 milestone

**Milestone 1**

> KoGES에서 실제 metabolic resilience phenotype을 만들 수 있는지 확인한다.

성공 조건:

```text
1. WC/SBP/FPG/TG/HDL이 충분한 wave에서 반복 측정
2. baseline + >=2 follow-ups usable N >=4,000
3. genotype intersection 충분
4. incident metabolic deterioration event 충분
5. physical activity usable variable 존재
```

이 단계가 성공하면 이 연구는 실제 수행 가능성이 높다.

---

# 38. 두 번째 분석 milestone

**Milestone 2**

> 공개 pQTL/GWAS만 이용하여 candidate proteins 5–30개를 생성한다.

필수 evidence:

```text
cis-MR
+
coloc
+
>=2 metabolic domains
```

가능하면:

```text
+ deCODE/UKB cross-platform replication
+ EAS direction concordance
```

---

# 39. 세 번째 분석 milestone

**Milestone 3**

> candidate protein GRS가 KoGES metabolic trajectory와 연관되는지 확인한다.

Primary model:

```text
MBI ~ time + protein_GRS + protein_GRS*time + covariates + random effects
```

Key statistic:

```text
protein_GRS × time
```

---

# 40. Decision framework

## Continue as main thesis

다음 중 ≥2개가 강하면 main thesis로 진행.

```text
A. 5개 이상 Tier-1 candidate proteins
B. EAS replication 존재
C. KoGES longitudinal association 존재
D. PA interaction 또는 MHO transition에서 추가 signal
```

## Keep CKD as backup

다음 상황이면 CKD가 더 안전하다.

```text
KoGES metabolic variables severe missingness
genotype-longitudinal N insufficient
candidate protein signals mostly lipid-only
EAS replication impossible
```

---

# 41. 현재 단계의 판단

현재 공개 자원 기준으로는 **GO** 쪽에 가깝다.

근거:

1. KoGES 안산·안성은 약 20년 repeated follow-up이 존재한다.
2. baseline genotype 대상이 8,198명 수준이다.
3. MetS 핵심 component를 구성할 임상 변수가 KoGES에 존재한다.
4. KoGES public GWAS는 72,298명에서 76 traits를 제공한다.
5. UKB-PPP와 deCODE라는 강력한 plasma pQTL source가 있다.
6. 2026년 CKB East Asian pQTL atlas가 추가되어 ancestry validation 환경이 개선되었다.
7. metabolic syndrome 자체의 proteomics 연구는 이미 있으나, **cross-ancestry genetically anchored protein discovery + Korean 20-year resilience validation** 조합은 차별화 여지가 있다.

다만 최종 GO 결정은 반드시 **KoGES raw variable availability / missingness audit 이후** 내려야 한다.

---

# 42. 즉시 다음 분석에서 확인할 항목

다음 대화/작업에서는 아이디어 논의를 멈추고 실제 서버 자료 기준으로 아래를 확인한다.

```text
[1] /srv/is-analysis/data/koges 실제 파일 inventory
[2] KoGES codebook 존재 여부
[3] participant ID key
[4] wave encoding
[5] WC/SBP/DBP/FPG/TG/HDL variable names
[6] medication variable names
[7] physical activity variable names
[8] genotype files and build
[9] available PCs / kinship
[10] actual repeated-measure N
```

이 결과를 이용하여:

```text
METABOLIC_RESILIENCE_KOGES_ANCHOR_PANEL.tsv
```

을 생성하고 Stage 0 audit script를 작성한다.

---

# 43. Recommended anchor phenotype v0.1

Stage 0에서 가장 먼저 구현할 최소 phenotype:

```text
waist
SBP
fasting_glucose
log_triglycerides
HDL
```

score:

\[
MBI =
mean[
z(WC),
z(SBP),
z(FPG),
z(logTG),
-z(HDL)
]
\]

이 5개는 biological interpretability와 longitudinal availability의 균형이 좋다.

추가 변수는 core phenotype을 깨지 않는 선에서 sensitivity로 확장한다.

---

# 44. Primary analysis hierarchy — 최종 권장안

```text
PRIMARY
└── Baseline healthy → incident MetS
    + cis-protein GRS

SECONDARY
└── Repeated continuous MBI
    + protein_GRS × time

SUPPORTIVE
└── EUR pQTL MR + coloc
    → EAS replication

EXPLORATORY
├── MHO → MUO
├── trajectory classes
├── physical activity interaction
├── sleep interaction
└── tissue / single-cell annotation
```

이 순서를 유지하면 결과가 일부 실패하더라도 논문 전체가 무너지지 않는다.

---

# 45. Key references and data resources

## Plasma pQTL

1. **UK Biobank Pharma Proteomics Project**
   - Sun BB et al. *Nature* (2023).
   - 54,219 participants; 2,923 unique proteins.
   - https://www.nature.com/articles/s41586-023-06592-6

2. **deCODE plasma proteome**
   - Ferkingstad E et al. *Nature Genetics* (2021).
   - 35,559 Icelanders; 4,907 aptamers.
   - https://www.nature.com/articles/s41588-021-00978-w

3. **CKB East Asian plasma proteogenomic atlas**
   - Pozarickij A et al. medRxiv (2026), preprint.
   - 3,965 Chinese adults; 7,289 proteins; 1,092 cis-pQTL proteins.
   - https://www.medrxiv.org/content/10.64898/2026.02.05.26345625v1
   - analysis resources: https://zenodo.org/records/21503247

## KoGES

4. **KoGES Cohort Profile Update**
   - International Journal of Epidemiology (2026).
   - https://academic.oup.com/ije/article/55/5/dyag192/8812758

5. **OPEN KoGES**
   - longitudinal cohort status / codebook / data information.
   - https://coda.nih.go.kr/openkoges-v2/

6. **KoGES public GWAS**
   - Nam K et al. *Cell Genomics* (2022).
   - N=72,298; 76 phenotypes.
   - https://pmc.ncbi.nlm.nih.gov/articles/PMC9903843/
   - PheWeb: https://koges.leelabsg.org/
   - summary statistics: https://zenodo.org/record/7042518

## Metabolic syndrome / obesity precedent

7. **Metabolic syndrome and plasma proteome**
   - *Cardiovascular Diabetology* (2021).
   - proteomics + incident/prevalent MetS + MR.
   - https://pubmed.ncbi.nlm.nih.gov/34016094/

8. **Korean MHO/MUO genetic architecture**
   - *Scientific Reports* (2021).
   - N=49,915.
   - https://pmc.ncbi.nlm.nih.gov/articles/PMC7838176/

9. **Korean obesity genetics**
   - GWAS in 93,673 Korean subjects.
   - https://pmc.ncbi.nlm.nih.gov/articles/PMC11359806/

10. **Korean obesity / metabolic syndrome criteria**
    - Korean Society for the Study of Obesity guideline.
    - BMI obesity ≥25 kg/m².
    - abdominal obesity ≥90 cm men / ≥85 cm women.
    - https://pmc.ncbi.nlm.nih.gov/articles/PMC10088549/

## Outcome GWAS

11. **T2D multi-ancestry**
    - >2.5 million participants; 428,452 T2D cases.
    - https://pubmed.ncbi.nlm.nih.gov/37034649/
    - downloads: https://www.diagram-consortium.org/downloads.html

12. **Global Lipids Genetics Consortium**
    - ~1.65 million participants across ancestries.
    - https://pmc.ncbi.nlm.nih.gov/articles/PMC8730582/

13. **MAGIC glycemic traits**
    - public fasting glucose / insulin / HbA1c / 2h glucose resources.
    - https://magicinvestigators.org/downloads/

14. **Large blood pressure GWAS**
    - up to 1,028,980 European participants.
    - SBP GCST90310294
    - DBP GCST90310295
    - PP GCST90310296
    - https://www.nature.com/articles/s41588-024-01714-w

15. **GIANT anthropometric data**
    - https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files

## Recent competing / adjacent research

16. **Sleep proteomics**
    - UKB N=43,709; 935 sleep-associated proteins; MR/coloc.
    - *Molecular Psychiatry* (2026).
    - https://pubmed.ncbi.nlm.nih.gov/42277229/

17. **Grip strength and healthy aging proteomics**
    - UKB N=27,828; proteomics subset N=3,366.
    - *Clinical Nutrition* (2025).
    - https://pubmed.ncbi.nlm.nih.gov/40886562/

18. **Exercise and proteomic aging**
    - UKB N=45,438 + 12-week intervention.
    - *npj Aging* (2026).
    - https://www.nature.com/articles/s41514-025-00318-w

---

# 46. Bottom line

이 프로젝트의 핵심은 새로운 “웰니스 score”를 만드는 것이 아니다.

핵심은:

> **대규모 인간 유전학으로 causal candidate proteins를 먼저 좁힌 뒤, 한국인의 20년 longitudinal data에서 실제로 건강한 대사상태를 오래 유지하는 현상과 연결하는 것**

이다.

따라서 분석의 방향은 다음 한 줄로 고정한다.

```text
Causal protein discovery
→ Cross-ancestry replication
→ Korean longitudinal metabolic resilience
→ Modifiable lifestyle interaction
```

이 구조를 유지하면 **proteomics, genetics, longitudinal epidemiology, precision wellness, Digital Health**를 하나의 석사 연구 안에 일관되게 연결할 수 있다.

---

# 47. Next action

다음 실행 단계는 **Stage 0 KoGES feasibility audit**이다.

가장 먼저 다음을 실제 서버에서 확인한다.

```bash
find /srv/is-analysis/data/koges -maxdepth 3 -type f -printf '%p\t%s\n' | sort
```

그 결과를 바탕으로:

1. KoGES file inventory
2. codebook mapping
3. metabolic anchor variables
4. wave availability
5. longitudinal eligible N
6. genotype intersection

을 계산하고,

```text
STAGE0_GO_NO_GO.md
```

를 생성한다.

이 Stage 0 결과가 본 연구의 실제 착수 여부를 결정한다.
