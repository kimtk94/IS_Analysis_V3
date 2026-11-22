# 석사 연구주제 제안 요약 (교수님 검토용)

**연구 가제**  
**심혈관–신장–대사 위험요인의 유전적 영향은 허혈성 뇌졸중 아형별로 다른가?**  
*Genetic heterogeneity of CKM-related liabilities across ischemic stroke subtypes*

## 연구 필요성 및 차별화

LDL-C, 혈압, 신장기능, 당뇨 등 CKM 관련 요인은 뇌졸중 위험과 밀접하지만, 이들 **개별 요인의 유전적 영향이 LAA(대혈관죽상경화), CE(심장색전), SVD(소혈관폐색)에서 실제로 다른지**는 비교 분석 설계와 공통 대조군에 대한 처리가 중요합니다. 동아시아·유럽 단백질 MR, CKM 대사체 MR·colocalization, 뇌졸중 아형 단백질 MR은 이미 발표돼 있어 **같은 방법의 단순 반복을 신규성으로 주장하지 않습니다.**

**핵심 가설:** CKM 구성요소별 유전적 효과는 뇌졸중 아형 간 이질성을 보이며, 일부는 분자 QTL과 혈관·면역세포 수준의 기전 근거로 연결될 수 있습니다.

## 연구 설계 (승인 조건부)

| 단계 | 수행 내용 |
|---|---|
| 1. 아형 GWAS 비교 | GIGASTROKE EUR AIS·LAA·CE·SVD 원본의 ancestry/build/allele/표본수 검증 → 단일형질 LDSC 유전력 진단 → 가능한 경우 유전상관 |
| 2. 정식 아형 대조 | **LDL-C, SBP**를 주노출로 시작해 LAA–CE–SVD 간 효과 차이를 공통 대조군 공분산을 고려해 검정. eGFR/UACR은 보조, T2D/BMI 등은 확장 |
| 3. 기전 검증 | 선정된 후보에 한해 cis-pQTL/eQTL MR, colocalization, GTEx·뇌/혈관 단일세포 데이터의 생물학적 맥락 점검 |
| 4. 동아시아 보조 분석 | GIGASTROKE EAS AIS와 TPMI Phecode 433.21을 **표현형 정의가 다르다는 전제 아래** 별도 민감도/방향성 분석 자료로 검토 |

## 현재 확보한 분석 근거 (실측)

TPMI SAIGE 원본 **15,588,353개 변이** 중 제공자 `QC=PASS`는 **8,042,797개**였습니다. GIGASTROKE EAS AIS와 GRCh38→GRCh37 변환 및 비회문 대립유전자 매칭을 수행해 **5,228,340개 SNP**가 정렬됐으며, `INFO≥0.8`, `MAF≥0.01`, `AF 차이≤0.05`를 동시에 적용하면 **3,745,225개**가 남았습니다. 전체 계산 및 독립 집계 재감사는 GitHub Actions에서 **PASS**했습니다.

- [실제 원본 분석 #44920288867](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44920288867)
- [독립 재감사 #44924780834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44924780834)

**해석 제한:** 이 수치는 GWAS summary statistics 간 좌표·allele 기반 기술적 교집합입니다. 질환 인과관계, 독립 재현, 유전력, MR·colocalization 결과는 아직 없습니다. TPMI `433.21`과 any ischemic stroke는 동일한 진단이 아니며, 개인별 표본 중복도 입증되지 않았습니다.

## 교수님 논의가 필요한 결정

**① 주제 선택:** CKM 구성요소별 뇌졸중 아형 이질성을 중심으로 할지, 기존 BBJ 기반 FGF5/ALDH2 등 유전자 기전연구를 중심으로 유지할지.  
**② 분석 범위:** LDL-C·SBP 중심의 작은 정식 검정으로 시작하고, eGFR/UACR 및 단백질·세포 분석은 후속 단계로 둘지.  
**③ 진입 조건:** 유럽 아형 GWAS의 표본중복·유효 표본 수·LD 참조 및 공분산 문제를 해결한 뒤 LDSC/MR을 실행할지.

**현재 판정: 연구제안·데이터 접근성 검토 GO / 인과·아형별 통계 결과 HOLD.**

상세 방법·선행연구 대조표: [연구제안서 본문](IS_CKM_STROKE_professor_proposal_20261122.md).
