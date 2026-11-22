# IS 학위 연구 과학적 검증 게이트 감사 — 2026-11-22

> **범위:** MASTER DEGREE 학술 연구 중 허혈성 뇌졸중(IS)만 포함. Novogene/Sales 대상 연구는 제외한다. 이 문서는 GitHub Actions **실행 결과**와 **과학적 결론**을 엄격히 분리한 read-only evidence audit다. 연구 수치·원본·기존 워크플로우를 변경하지 않는다.
>
> **주의:** 기존 **BBJ 중심 IS 연구**와 최근 **CKM→Stroke 교수 제안서 P0–P6**는 **서로 다른 분석 파이프라인**이다. P0–P6의 과학적 0/7을 IS 전체의 0%라고 표현하면 안 된다.

## 1. 이번 감사에서 직접 검증한 원본

| 근거 | 관찰한 실측/로그 | 해석 |
|---|---|---|
| [CKM–Stroke P0–P6 run #45007231962](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45007231962) (2026-11-22, commit `dcb1cfc4b5f74e8ab70a6d079c304c5292c8a4cd`) | 11/11 jobs SUCCESS; aggregate: P0 technical PASS/scientific HOLD, P1 HOLD/HOLD, P2–P6 HOLD/NOT_RUN. **Scientifically validated stages: 0/7.** | 워크플로우 구동 성공, 실제 P2–P6 계산 완료 아님 |
| [TPMI 433.21 N gate #44889438772](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889438772) | TPMI PheWeb N 313,909, PLINK2 OBS_CT 204,099; discrepancy 109,810; `approved_Neff=false` | 표본 정의 및 SAIGE/PLINK2 모델 차이를 해결하지 않은 채 임의 `N_eff` 채택 금지 |
| [SAIGE full N run #44907119885](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44907119885) | 원본 15,588,353 rows: 9,249 cases + 304,660 controls = 313,909; 8,042,797 source-QC PASS variants; scientific HOLD | **SAIGE 소스**에서 확인한 표본 수이며, PLINK2의 204,099와 자동 통합하지 않음 |
| [SAIGE–GIGASTROKE original overlap #44920288867](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44920288867) + [independent conservation audit #44924780834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44924780834) | 5,228,340 allele-matched; 3,745,225 joint INFO≥0.8, MAF≥0.01, AF-Δ≤0.05; technical audit PASS, scientific HOLD | 변이 매칭·보존성만 증명. phenotype 비교 가능성·참가자 비중복·독립 검증·LDSC/MR/coloc 승인은 증명하지 않음 |
| [IS P3 ABF vs SuSiE gate](https://github.com/kimtk94/IS_Analysis_V3/blob/main/docs/is/IS_P3_ABF_VS_SUSIE_GATES_20261017.md) (2026-10-18 snapshot) | 36/393 source-flagged ABF 후보 그룹, 25/36에 제한된 수치 조건 통과 행 존재. **0/36 coloc.abf eligible/executed** (당시) | 원본 source FASTA 플래그는 독립적 allele/FASTA 검증이 아니며, region 전체 통계·sample N·case fraction·priors 필요 |
| [IS P1.5 replication decision](https://github.com/kimtk94/IS_Analysis_V3/blob/main/docs/is/IS_P15_INDEPENDENT_REPLICATION_DATASET_DECISION_20261018.md) (2026-10-18 snapshot) | BBJ-excluded EAS-AIS 독립 GWAS summary statistics 미확인; 확정 독립 재현 유전자 0 | GIGASTROKE에 BBJ 포함 위험. 유사 방향성은 **descriptive concordance**로만 표기 |
| [IS_MASTER.md](https://github.com/kimtk94/IS_Analysis_V3/blob/main/docs/IS_MASTER.md) (updated 2026-10-05 snapshot) | FGF5, ALDH2, SH3PXD2A, COL4A2에 대한 기존 역사적 결과 기재 | 기존 연구의 보고 값. 최신 신규 CKM–Stroke P0–P6의 과학적 PASS나 독립 재실행으로 오인 금지 |

## 2. 방법별 판정 (2026-11-22 GitHub evidence scope)

| 방법/결과층 | 판정 | 근거와 아직 필요한 조건 |
|---|---|---|
| TPMI SAIGE **원본 GWAS QC 및 allele overlap** | **TECHNICAL PASS** | 8,042,797 source-QC PASS, 5,228,340 matched, stringent joint subset 3,745,225; 임상 동등성/independence는 HOLD |
| EAS vs EAS stroke **독립 replication** | **HOLD** | BBJ inclusion/cohort overlap과 phenotype ascertainment를 원문 수준으로 해결해야 함. 독립 재현으로 승격 금지 |
| CKM–Stroke **P0 study manifest** | **TECHNICAL PASS / SCIENTIFIC HOLD** | 14개 phenotype-registry row는 schema 통과하나, EUR exposures의 accession/build/data_file 다수 미완성 |
| CKM–Stroke **P1 원본 수집** | **HOLD** | 11개 core input 미확보, 3개 EAS outcome optional 누락. 영업 데이터로 대체 금지 |
| CKM–Stroke **P2–P5 통계 수행** | **NOT_RUN** | approved tasks의 argv, inputs, outputs 부재. MR/coloc 결과를 생성했다고 보고 금지 |
| CKM–Stroke **P6 replication** | **NOT_RUN** | `data/ckm_stroke/replication.csv` 없음. 존재 여부만으로 독립성 보장되지 않음 |
| **LDSC h² / LDSC-SEG** | **HOLD / 계산 검증 없음** | SAIGE·PLINK2 간 N 정의와 승인된 SNP별 N/N_eff, outcome phenotype, 효과 단위, signed LD reference, HapMap3/common-SNP fit 및 주요 LDSC 통계 확인 전 h²를 실행·보고하지 않음 |
| 기존 BBJ IS **MR / coloc / SuSiE** | **HISTORICAL / 부분 근거; 신규 검증과 구분** | 기존 FGF5 ABF H4≈0.779, SuSiE H4≈0.752라는 문서 기록이 있음. 정확한 원 결과·입력·LD provenance 재조회 없이는 독립 검증 완료 또는 '인과 유전자 확정'으로 격상 불가 |
| 기존 BBJ IS **cell-type context** | **PARTIAL** | tissue/cell expression vs cell-eQTL/ATAC regulatory support 분리. 발현만으로 기능적 인과 주장 불가 |

## 3. 우선순위가 높은 최소 후속 작업 (analysis sprawl 금지)

1. **P0 — 정확한 데이터 모델 확정:** TPMI SAIGE 313,909과 PLINK2 OBS_CT 204,099의 모집단/필터/분석 모델 차이, 표본 수의 column semantics, case/control, GIGASTROKE phenotype 및 overlap roster를 **각 원본 출처 기준으로** 묶어 승인 여부 기록. N_eff 임의 대입 금지.
2. **P1 — LDSC input preflight only:** 검증된 independent effective sample size 또는 문헌/소스가 요구하는 per-SNP sample-size prescription 확보, Z score / allele / variant QC, LD reference ancestry/build, regression SNP set을 체크하고 PASS 전에는 h²/LDSC-SEG 결과 생성 금지. h² estimator에 필요한 입력의 적합성을 별도 검증.
3. **P1 — 기존 IS MR/Coloc 근거 재확인:** BBJ 기존 결과의 input hash, analysis commit, universe attrition, 4 illustrative loci의 regional overlap, allele/MAF, GWAS & QTL N, coloc priors 및 LD diagnostics를 실측으로 재조회하여 gene-wise evidence matrix로 정리. 완결되지 않은 36 P3A 후보 그룹을 무조건 재분석하지 말고 **조건 충족 그룹만** 검토.
4. **P2 — 논문용 제한적 후보 우선순위:** 역사적 유전체 결과와 신규 CKM–Stroke은 통합 추정치를 만들지 말고 계층과 출처를 분리. 교수님 제안서에는 '실측 수행', '부분 QC', '향후 검증'을 분리.

## 4. 재감사 체크

- [x] CKM–Stroke P0–P6 latest runner aggregate 실제 로그 확인
- [x] TPMI SAIGE 원본 N/QC 및 allele overlap GitHub Actions 로그 확인
- [x] 이전 P3/독립 replication 결정·주장 범위 확인
- [ ] LDSC 실행·유전력 수치: **미승인**; 근거와 입력 미충족
- [ ] CKM–Stroke MR/coloc: **NOT_RUN**; approved actual source inputs 미구성
- [ ] 역사적 BBJ IS 대상 4 loci MR/coloc 원 결과 artifact independent readback

**결론:** 이 문서는 audit snapshot이며 새 과학적 결과의 통과를 선언하지 않는다. GitHub Actions 녹색 체크와 과학적 완료는 서로 다른 상태값으로 유지한다.
