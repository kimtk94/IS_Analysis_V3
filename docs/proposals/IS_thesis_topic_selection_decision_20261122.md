# 석사 논문 주제 선정 검토 메모 — IS 기전연구 vs CKM–뇌졸중 아형

**기준일:** 2026-11-22 · **용도:** 지도교수 미팅 전 주제 우선순위 정리  
**판단:** 본 문서는 연구 기획상 판단이며, 과학적 결과의 신규 검증을 주장하지 않는다. 연구별 증거 수준을 구분한다.

## 1. 현재 우선순위

| 항목 | 1순위: BBJ 기반 IS 분자·세포기전 | 2순위: CKM 구성요소별 뇌졸중 아형 이질성 |
|---|---|---|
| 연구질문 | 허혈성 뇌졸중 연관 좌위가 어떤 유전자·세포 맥락에서 기능할 수 있는가 | LDL-C·SBP 등의 유전적 영향이 LAA/CE/SVD 간 실제로 다른가 |
| 이미 기록된 근거 | `docs/IS_MASTER.md`는 4개 핵심 좌위(FGF5·ALDH2·SH3PXD2A·COL4A2)와 기존 ABF/SuSiE 결과를 **기록**함 | GIGASTROKE EUR 아형 공개 원본과 CKM GWAS 구성은 기획됨. 동아시아 TPMI SAIGE와 GIGASTROKE EAS AIS간 엄격한 allele 교집합 **3,745,225 SNP 기술적으로 검증** |
| 이미 실행한 과학적 결과의 수준 | 과거 연구 기록상 FGF5 ABF H4≈0.779, SuSiE H4≈0.752; ALDH2 rs671 신호. **이 수치의 코드·원본 재실행과 LD 적합성은 본 메모에서 재검증하지 않았음** | P0–P6 공식 집계상 과학적 PASS **0/7**; 현재 LDSC h², 아형별 공분산 보정 검정, MR, coloc 결과 없음 |
| 주요 병목 | 동아시아 좌위와 세포·조직별 cis-QTL/ATAC의 donor-aware 재현성, cohort-matched LD 부재, 후보들에 대한 다중 신호 검증 | 아형 결과 GWAS 표본중복·메타분석 `N_eff`·상관 효과 공분산 미확인; 일부 EUR 노출 GWAS 아직 권한/원본 검증 전 |
| 선행연구 중복 | FGF5·ALDH2·pQTL의 EAS/EUR MR은 2025 논문에 이미 등장. **기존 알려진 표적 반복은 신규성이 낮음** | CKM–stroke MR/coloc 및 protein MR subtype 분석도 이미 발표됨. 아형별 정식 CKM 구성요소 효과 대조가 실제 차별화 가능한지 검증 필요 |
| 석사 논문 완성 가능성(정성적 판단) | **상대적으로 높음** — 후보 좌위·이전 분석 기록이 있어 좁은 질문으로 재검증 가능. 단, 실제 원본 재현·기전 검증이 있어야 함 | **현재 불확실성 큼** — 정식 아형 간 차이를 검정하기 위한 자료와 공분산 전제가 해결되지 않으면 중심 가설 자체가 HOLD |
| 추천 역할 | **우선 검토 주제** | **백업 또는 별도 후속 연구** |

**기본 권고:** 석사 논문 **1순위는 BBJ–IS 기전연구**로 두되, 단순 FGF5/ALDH2 재발견이 아니라 **SH3PXD2A 면역·혈관세포 조절 / COL4A2 혈관 mural-cell 조절**을 중심으로 *독립적인 세포 맥락 증거를 확보할 수 있는지* 먼저 판단한다. FGF5와 ALDH2는 잘 알려진 신호에 대한 양성 대조·비교 기전축으로 둔다. 현재 SH3PXD2A/COL4A2의 이전 bulk ABF H4≈0.351/0.331은 결정적인 coloc 증거가 아니다.

**2순위 CKM**은 과학적 전제가 충족될 때 **LDL-C·SBP × LAA/CE/SVD의 정식 효과 차이**에 집중한다. TPMI Phecode 433.21 ↔ GIGASTROKE EAS AIS 3.745M SNP는 데이터 호환성 증거일 뿐 **아형 분석 성공이나 독립 재현 성공의 증거가 아니다.**

## 2. 신규성을 확보하려면 무엇을 보여야 하나

### A안 — BBJ–IS 기전연구 (추천)

- **가설:** 특정 IS GWAS 좌위(특히 `SH3PXD2A`, `COL4A2`)가 혈관 내피·혈관평활근/주위세포와 면역세포에서 서로 다른 조절·발현 맥락을 가질 수 있다.
- **최소 핵심 그림:** (1) 검증된 BBJ EAS 신호 및 독립 locus QC, (2) GTEx/뇌·혈관 QTL의 conditional/coloc 적합성, (3) 독립 sc/snRNA + ATAC에서 donor-aware cell-type 맥락, (4) 좌위별 증거 등급과 반례/음성 결과.
- **차별화 검증:** 유전적 신호→세포 특이성의 연결이 기존 논문에서 이미 보고된 수준을 넘어서는지 문헌 단위로 확인. 유전자 발현만 확인한 경우 '원인 세포'로 단정하지 않음.
- **GO 조건:** 핵심 좌위 1–2개에서 원본 GWAS·해당 조직 QTL과 세포 자료의 **독립 추적·접근성·표본 기원**을 증빙할 수 있음.
- **HOLD/축소 조건:** cohort-matched LD 부재로 fine-mapping/coloc이 불안정하고, 세포 문맥 근거가 단순 expression 수준에 머무르면 결과를 탐색적으로 표시. FGF5/ALDH2 재발견만 남으면 신규성 주장 보류.

### B안 — CKM 구성요소 × stroke subtype

- **가설:** LDL-C·SBP 등 구성요소의 stroke 아형별 유전적 효과가 동일하지 않다.
- **최소 핵심 그림:** (1) 유럽노출·결과 GWAS P0/P1 검증, (2) 표본중복·메타분석 N 정책, (3) 아형별 효과와 **효과 차이의 직접 신뢰구간/검정** (공분산 지원), (4) 검증된 후보만 추가 MR/coloc/세포 문맥.
- **GO 조건:** 원본 EUR LDL-C/SBP 노출 사용승인, 아형별 GWAS의 정확한 공분산 또는 방어 가능한 독립 차이 검정법, LD-reference 적합성, h² 신뢰도.
- **HOLD/축소 조건:** 공유 대조군·N_eff 미해결이면 아형 차이 결론 불가. LAA/SVD LDSC 정밀도 제한시 AIS/CE로 축소할 수 있으나 연구 질문의 중요 부분이 변하므로 교수님 승인 필요.

## 3. 제안된 교수님 미팅 아젠다

1. **방향 선택:** A안(기존 BBJ–IS 분자·세포기전)을 석사 논문으로, B안(CKM 아형)을 예비 후속 프로젝트로 둘지.
2. **구체적인 1차 평가:** A안은 `SH3PXD2A`와 `COL4A2` 중 세포/조직 QTL·ATAC 독립 자료가 더 충실한 쪽을 1차 표적으로 선택할지. FGF5/ALDH2는 외부 선행연구와의 기준점.
3. **시간·범위 원칙:** 사용자 생산 검체·메일 문의는 후순위. 공개 데이터와 GitHub Actions 재현·원본 QC를 먼저 해결.

## 4. 이번 작업에서 직접 확인한 CI 상태

- [CKM P0–P6 workflow #45007231962](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45007231962) **SUCCESS** after adding pinned SciPy for an existing numerical unit test. **391/391 tests PASS**, P0 technical PASS/scientific HOLD, P1 technical HOLD/scientific HOLD, P2–P6 technical HOLD/scientific NOT_RUN; **scientific 0/7**. Do not reinterpret CI success as LDSC/MR progress.
- [TPMI SAIGE 15,588,353 original rows and 8,042,797 QC PASS](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44907119885) completed full-source check.
- [SAIGE/EAS allele join](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44920288867) / [independent audit](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44924780834): **3,745,225** matched nonpalindromic variants under combined INFO≥0.8, MAF≥0.01, |ΔAF|≤0.05; no claim of biological replication.

## 5. 직접 참조한 기록과 선행연구

- Existing study state: [docs/IS_MASTER.md](../IS_MASTER.md) — **last updated 2026-10-05**, so research-stage statuses are historical source statements and need independent latest run verification before final thesis methods.
- CKM data and novelty audit: [docs/ckm_stroke_feasibility_20261116.md](../ckm_stroke_feasibility_20261116.md) and [docs/ckm_stroke_gigastroke_cohort_overlap_neff_20261117.md](../ckm_stroke_gigastroke_cohort_overlap_neff_20261117.md).
- Prior 2025 East Asian/European IS pQTL-MR: https://pubmed.ncbi.nlm.nih.gov/40304040/ (FGF5, ALDH2 etc already reported).
- 2026 proteogenomic stroke integration: https://www.nature.com/articles/s43856-026-01734-z
- 2026 CKM metabolic stroke MR: https://pubmed.ncbi.nlm.nih.gov/42290487/
- 2022 GIGASTROKE original: https://www.nature.com/articles/s41586-022-05165-3

**Decision remains conditional:** This is a ranked **study-plan recommendation**, not a claimed statistical power calculation, independent validation, or final dissertation authorization. The ranking should be revisited if current independently validated molecular-QTL and cell-type evidence of A is lacking, or the EUR covariance/N-eff gate for B is resolved.
