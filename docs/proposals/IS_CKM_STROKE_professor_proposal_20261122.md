# 교수님 연구제안서 초안 — CKM 구성요소와 허혈성 뇌졸중 아형별 유전적 이질성

> 문서 상태: **PROPOSAL_READY / DATA_FEASIBILITY_PASS / SCIENTIFIC_RESULTS_HOLD**  
> 작성 기준일: 2026-11-22 · 연구 브랜치: `research/ckm-stroke-staged-actions-20261111`  
> 용도: 석사 연구주제 논의용 1차 제안서. **이미 수행된 GWAS 기술 QC와 추후 수행할 과학적 검정을 명확히 구분한다.**  
> 범위 주의: 이 문서는 `docs/IS_MASTER.md`에 기록된 BBJ 기반 FGF5·ALDH2·SH3PXD2A·COL4A2 연구를 대체하지 않는 **CKM→stroke subtype 후보 연구설계**다.

## 1. 가제

**심혈관–신장–대사 위험요인의 유전적 영향과 허혈성 뇌졸중 아형(LAA·CE·SVD)의 이질성: 유럽·동아시아 GWAS 및 분자 QTL 통합 연구**

English: **Genetic heterogeneity of cardiometabolic and kidney-related liabilities across ischemic stroke subtypes: integrative GWAS and molecular-QTL analyses**

## 2. 연구배경과 기존 문헌 대비 공백

CKM (cardiovascular–kidney–metabolic)은 단일 유전형질이나 특정 단계가 아닌 임상적 위험요인 묶음이다. 뇌졸중 자체는 CKM 임상 분류와 얽힐 수 있으므로 본 연구는 **CKM stage 4 → stroke 인과효과**를 직접 추정한다고 주장하지 않는다. 대신 측정 가능하고 정의가 명확한 LDL-C, 수축기혈압(SBP), eGFR/UACR, T2D/BMI 등 **구성요소의 유전적 연관성**을 평가한다.

이미 발표된 연구는 (a) 동아시아·유럽 단백체 기반 허혈성 뇌졸중 MR, (b) 단백체 MR+coloc과 LAA/CE/SVD 비교, (c) CKM 대사체 시그니처와 뇌졸중 MR+coloc, (d) 대규모 CKM 유전구조 연구까지 포괄한다. 따라서 **“EAS/EUR + pQTL MR + coloc” 자체**를 신규성으로 주장할 수 없다.

**검증할 차별화 질문:** CKM의 개별 구성요소가 공유하는 유전적 효과와 아형별로 차이를 보이는 효과가 **LAA, CE, SVD 간 정식 통계적 대조에서도 유지되는가?** 그 차이가 혈관벽·면역·대사 분자기전과 합치하는가? 이 질문도 완전한 세계 최초 주장으로 취급하지 않으며, 최신 경쟁 문헌과 중복 분석 여부에 대한 표를 유지한다.

### 반드시 비교할 선행연구

| 문헌 | 이미 수행한 내용 | 본 제안과의 구별점/경고 |
|---|---|---|
| Mishra et al., Nature 2022 (GIGASTROKE), doi:10.1038/s41586-022-05165-3 | 다중 조상집단 뇌졸중 GWAS, 아형과 다중 오믹스 | 아형 차이 자체 및 ancestry 차이는 이미 일부 보고됨 |
| *Proteome-Wide Genetic Study in East Asians and Europeans…*, PubMed 40304040 (2025) | EAS/EUR cis-pQTL MR로 FGF5 등 후보 | 기존 후보 15개 재사용만으로 신규성 부족 |
| Lei et al., JCBFM 2026, doi:10.1177/0271678X261462718 | CKM 대사체–뇌졸중 연관성, MR·coloc | CKM과 stroke 연결 자체는 선행 |
| Song et al., Transl Stroke Res 2026, doi:10.1007/s12975-026-01509-7, PMID 42911303 | protein MR + stroke subtype + coloc | 단백질별 subtype MR/coloc만으로 신규성 부족 |
| *Integrative proteogenomic analyses…*, Communications Medicine 2026, PMID 42399650 | 광범위 pQTL·subtype·인과 후보 표적 | 분자기전의 탐색은 대조·재현 근거 필요 |

참조 URL: https://pubmed.ncbi.nlm.nih.gov/40304040/ · https://pubmed.ncbi.nlm.nih.gov/42911303/ · https://pubmed.ncbi.nlm.nih.gov/42399650/ · https://doi.org/10.1177/0271678X261462718

## 3. 연구가설

- **H1 (주가설):** CKM 구성요소의 유전적 효과는 허혈성 뇌졸중 아형 간에 동일하지 않다. 사전 지정한 `LAA vs CE`, `LAA vs SVD`, `CE vs SVD` 대조에서 **공유 대조군/표본 중복에 따른 공분산을 고려한 정식 검정**을 수행한다.
- **H2 (보조):** 아형 특이성을 보인 유전적 신호 중 일부는 cis-pQTL/eQTL과 동일한 원인 변이 후보를 공유할 수 있다. 그러나 MR과 colocalization은 **각각 다른 인과 전제**를 갖는다.
- **H3 (탐색):** 동아시아 AIS 관련 자료에서 일부 방향성이 유지될 수 있다. 이는 검증할 가설이지 이미 확인된 효과 재현 결과가 아니다. **TPMI Phecode 433.21은 LAA·CE·SVD 정의에 해당하지 않는다.**

## 4. 구체적 목적 및 통계 분석 설계

### Aim 1 — 아형·구성요소별 유전적 신호의 비교 가능성 확보

**주요 결과:** GIGASTROKE EUR AIS `GCST90104540`, LAA `GCST90104542`, CE `GCST90104541`, SVD `GCST90104543`. 모든 파일은 원본 GRCh37, effect allele, log(OR), SE, P 및 연구별 표본 기여를 확인한다. **CE·LAA accession 반전 금지.**

**노출(우선순위 제한):** 주분석 `LDL-C` 및 `SBP` 두 구성요소; 이차 분석 `eGFR`, `UACR`; 확장 분석 `T2D`, `BMI`, `HDL-C`, `TG`. 서로 다른 단위(연속형·log/정규화·이분형)와 효과 방향을 명시하고, 검증되지 않은 공개 EUR 요약통계는 아직 입력으로 승인하지 않는다.

**사전 조건:** 원본 파일·조상집단·좌표·효과 대립유전자·중복/회문 SNP·INFO·P-value/SE QC, 노출–결과 참여 집단 중복 위험, 적절한 EUR LD 참조집단, outcome 메타분석에서 **실제 사용 가능한 (N_{mathrm{eff}}) 근거**를 확인한다. 원자료의 total N이나 사례/대조군 조화평균을 임의 LDSC N으로 고정하지 않는다.

**단계별 평가:** 승인된 입력부터 단일 형질 LDSC `h²`/SE, mean χ², intercept 및 regression SNP count를 먼저 보고, 불안정한 LAA/SVD의 `rg` 추정을 중단하거나 탐색적으로 제한한다. 충분한 유전력 및 참조 SNP 정합성을 만족하는 경우에만 다른 구성요소와 `rg`를 수행한다.

### Aim 2 — 아형 간 **정식 차이** 및 조건부 유전적 관련성

LD-pruned instruments를 노출 GWAS에서 고르고, 도구강도·독립성·다효성 가정을 평가한 다음, 조상집단과 효과방향을 맞춘 분석을 수행한다. 단순히 `LAA 유의, CE 비유의`를 **두 아형 효과가 다르다는 증거**로 쓰지 않는다. 아형 간 `β_LAA−β_CE` 등의 대조는 가능한 **공통 대조군·상관 GWAS의 공분산을 반영**한다. 필요한 공분산의 신뢰 가능한 추정이 없다면 subtype 차이에 대한 정식 결론을 **HOLD**한다.

MR은 도구 수에 따른 IVW/weighted median/MR-Egger 등 적절한 방법 및 leave-one-out/heterogeneity/directionality 검증을 조건부로 적용하며, 해석에 충분한 instrument가 없는 경우 **인과 추정으로 확대하지 않는다**. LDL-C와 SBP의 상관 및 다효성은 관찰/조건부 추정에서 나누어 다루고, 다변량 MR을 하려면 조건부 F-statistic과 적합한 도구집합을 사전에 확보한다. 다중비교 FDR 정책은 비교 대상 수를 확정한 후 등록한다.

### Aim 3 — 분자 QTL·세포 유형에서 메커니즘 우선순위 평가

Aim 2에서 정식 통계 검증을 통과한 후보에 한해 UKB-PPP·deCODE·GTEx 등 cis-pQTL/eQTL과 coloc/필요시 multi-signal 접근을 수행한다. eQTL/pQTL의 조직/세포 출처, genomic build, allele, LD 및 overlap을 구분한다. rs671/ALDH2처럼 coding-driven 가능성이 큰 후보는 bulk-eQTL 결손을 곧바로 생물학적 부정으로 해석하지 않는다.

이전에 검토한 `FGF5`, `ALDH2`, `SH3PXD2A`, `COL4A2`는 **우선 검토 후보**일 뿐, 이번 CKM-아형 제안에서 재현된 최종 표적이라고 적지 않는다. snRNA/ATAC는 후보세포의 **표현·조절 가능성**을 설명하는 후속 해석층이며, 단일세포에서 MR 효과를 직접 검증하는 것으로 설명하지 않는다.

## 5. 이번에 실제 확인한 동아시아 데이터 타당성 — 검증 결과만 기재

| 실측 항목 | 검증 결과 |
|---|---:|
| TPMI SAIGE Phecode 433.21 원본 전체 변이 | 15,588,353 |
| SAIGE `QC=PASS` | 8,042,797 |
| SAIGE 회귀분석 환자·대조군 수 | **9,249 cases / 304,660 controls** (전 변이 일관) |
| GIGASTROKE EAS AIS discovery | **19,032 cases / 237,242 controls** |
| 실제 GRCh38→GRCh37 좌표 + 비회문 allele 일치 | 5,228,340 |
| 추가 `INFO ≥0.8`, `MAF ≥0.01`, `|ΔAF|≤0.05` 동시 통과 | **3,745,225** |
| 엄격한 교집합에 대한 독립 계산 보존 재감사 | **PASS** |
| EAS 독립 인과 효과·대립유전자별 복제·MR·coloc | **NOT RUN / HOLD** |

근거: [whole-genome Action #44920288867](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44920288867) **attempt 2 SUCCESS**; [independent conservation audit #44924780834](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44924780834) **attempt 2 SUCCESS**. SAIGE 효과는 **Allele2 기준**, TPMI PheWeb 좌표는 **GRCh38**, GIGASTROKE EAS 원본은 **GRCh37**. 원본 샘플 독립성은 참여 코호트 명단만으로 참가자 단위 0명 중복임을 증명할 수 없다.

**주의:** TPMI `433.21` (*Cerebral artery occlusion, with cerebral infarction*)은 GIGASTROKE *any ischemic stroke*와 정의가 완전히 동일하지 않고 TOAST 아형도 아니다. 따라서 **동아시아 SNP 공통분모의 기술적 검증**에 사용하고, 별도 phenotype sensitivity 데이터로만 검토한다. TPMI의 다른 `PLINK2 OBS_CT=204099`는 SAIGE N_case/N_ctrl로 교체할 수 없다.

## 6. 의사결정 기준 및 예상 결과물

| Gate | 진입 조건 | 승인 시 결과물 | 현 상태 |
|---|---|---|---|
| P0 데이터 접근·원본 기술 QC | 출처/파일해시/build/효과대립유전자 및 allele 매칭 | 데이터 목록, reproducibility evidence | **PASS (위 EAS 기술검증 범위에 한함)** |
| P1 아형 outcome / exposure 사용승인 | 수정본 accession·전장 GWAS 원본·조상집단·표본중복·LD·실제 N 정책 | 분석 준비표 | **PARTIAL / HOLD** |
| P2 단일형질 LDSC | P1 PASS + 실제 regression SNP / h²·intercept 진단 | h² QC 및 불안정 trait 표시 | **NOT RUN** |
| P3 아형 차이 검정 | 공분산 지원, 효과 비교가능성, 검정력 | formal contrast + effect/CI/FDR | **NOT RUN** |
| P4 MR/coloc | IV strength, cis source, LD reference, pleiotropy 등 | instrument table, MR/coloc evidence | **NOT RUN** |
| P5 분자·세포 문맥 | P4에서 제한된 후보, 독립 데이터 | candidate mechanism figure | **NOT RUN** |

**제안서 단계 GO:** 데이터 접근 및 계산 재현성의 근거가 있고 구체적 차이 검증 가설이 있으므로 교수님과 주제 적합성 검토가 가능하다.

**논문 결과 주장 HOLD:** 현재 결과는 원본 GWAS의 allele/좌표 교집합이지, `subtype heterogeneity`, pQTL 인과관계, LDSC h², cell-type mechanism을 검정한 결과가 아니다.

**주제 중단·축소 기준:** (a) 노출 GWAS·stroke GWAS의 표본중복/메타분석 N 전제가 해소되지 않는 경우, (b) LAA/SVD 유전력 신호가 불안정해 아형 차이 검정의 정밀도가 확보되지 않는 경우, (c) 이미 출판된 CKM-subtype 분석과 차별화되는 통계 대비가 성립하지 않는 경우. 이때 **AIS/LAA/CE/SVD 전체로 무리하게 확장하지 않고** 하나의 구성요소·아형 또는 기존 BBJ IS 기전연구로 축소한다.

## 7. 교수님께 확인할 핵심 결정 3가지

1. 이번 논문의 **중심**을 `CKM 구성요소×stroke subtype 정식 이질성 검정`으로 할지, 기존 BBJ `FGF5/ALDH2` 기전연구를 주연구로 유지할지.
2. 통계적 1차 노출을 **LDL-C 및 SBP**로 제한하고 신장·당뇨 관련 요인은 보조로 둘지.
3. EUR subtype 결과를 주분석으로 두되 EAS AIS 및 TPMI를 **표현형 불일치가 있는 보조 탐색층**으로 다루는 것에 동의하는지.

### 참고: 파일·자동화 범위

이 문서는 **연구기획 및 검증 상태의 스냅샷**이다. 데이터/원본 GWAS/사례별 결과값을 Git에 넣거나 MR을 이미 실행했다고 작성하지 않는다. 기존 P2–P5의 승인되지 않은 입력 계획은 자동 실행하지 않는다. 표본·검체 관련 별도 연락이 필요한 사안은 본 제안 단계의 필수 진입조건에서 제외한다.
