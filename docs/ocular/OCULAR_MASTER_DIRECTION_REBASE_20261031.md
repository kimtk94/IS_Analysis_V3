# OCULAR MASTER — 연구 방향 재설정 (2026-10-31)

## 2026-11-01 신규성 및 논문 원고 검증 — 이전 자료의 출판 수준 재분류

[선행연구별 Claim Overlap / 신규성 / Erratum 감사](THESIS_NOVELTY_CLAIM_OVERLAP_AND_ERRATUM_AUDIT_20261101.md)와 [영문 초록·Introduction·Methods·Results·Discussion 전체 초안](THESIS_MANUSCRIPT_METHODS_RESULTS_DISCUSSION_DRAFT_20261101.md)을 작성했다. 새로운 원시 데이터 분석은 수행하지 않았으며 **기존 검증된 Figure/edgeR/교차연구 산출물에만 의존**한다.

**출처 관련 핵심 수정:**
1. Singh et al. (PMID 41890121, DOI 10.64898/2026.03.12.711121)은 검색·PMC 검증 기준 **2026년 3월 bioRxiv 프리프린트**이다. 현재 단계에서 동료심사 저널 논문으로 설명하면 안 된다.
2. Mehrotra et al. 100nM 연구는 **2026년 2월 IOVS 정식 출판**(PMID 41636427, DOI 10.1167/iovs.67.2.13)이다. 그보다 앞선 bioRxiv 버전은 PMID 41509484. **2026년 10월 1일 정정 공고**(PMID 42895266, DOI 10.1167/iovs.67.12.34)가 실제 있으나 정정 원문/수정 보충표 영향을 확인하지 못해 **외부 결과 최종 출판본 검증 HOLD**다. 정정 내용 확인 전에는 이전 그림 수치를 확정된 정정 후 보충표라고 주장하지 않는다.
3. GSE313736는 공개된 GEO 데이터세트임을 확인했지만 **63↑·17↓=80개 DEG의 저자 원본 유전자 리스트와 대응하는 독립 동료심사 본문은 미확인**이다. 우리의 83/88 집합을 저자의 정확한 80개로 단정하거나 중복집계하지 않는다.
4. DEX TM 전사체·장기 노출 효과는 2022년 논문, TM ROS/TGFβ/NOX4 기전은 2023년 연구, ECM·시간별 변화는 Mehrotra 2026 논문, TM chromatin/GWAS는 Singh 2026 프리프린트에서 이미 분석했다. **새 기전 발견 주장은 STOP**, 실제 검증된 **독립 실험 간 전사 방향 재현성·동일 donor 분류 강건성·재현 가능성 감사**만 학위논문 주요 기여로 남긴다.

**원고 진행 판정:** Methods/Results/Discussion 초안 작성 GO, 최종 원고 제출 **HOLD** (2026-10-01 Erratum 원문 영향 검증, 원 저자 80개 리스트·방법 대응, 지도교수와 학위논문 수준 신규성 합의). 전체 표·그림 설명은 저널 출판 수준과 GEO 자료 상태를 혼동하지 말 것.

## Figure·논문 Results 1차 패키지 — 실제 GitHub Actions 검증 완료 (2026-10-31)

기존 검증 산출물만 재사용한 [최종 Figure CI #41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011) **SUCCESS**. 원본 SHA256 고정 확인 후 3개 PNG와 3개 SVG, 16개 후보 Supplementary Table S1, 출처 manifest, 해석용 Figure legends 생성. [artifact #13045363678](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011/artifacts/13045363678) 직접 다운로드 및 압축 CRC 확인, **그림 3개를 열어 본문과 하단 주석의 겹침이 없도록 수정**한 최종본을 기준으로 사용. 초기 도표 run #41208620110은 CI 자체는 성공했으나 하단 설명문이 겹쳐 최종본에서 제외.

- **Fig 1:** GSE313736 500 nM 48h / Mehrotra 100 nM 48h / GSE301517 500 nM 7d. MYOC↑, TGFBI↓, MMP2↓ 세 연구 일치, ACTA2는 7d에서 반대 방향. 7d CPM+1 로그 비율과 48h edgeR 로그 효과크기는 **동일 추정치가 아니므로 방향을 중심으로 그렸다**.
- **Fig 2:** 내부 57명 BeamA confidence .75/.50, 동일 55명 BeamA/scPred TM_Fibro. 핵심 3개 방향·유전자별 FDR 안정성. **동일 원자료 재분석이지 독립 복제가 아님**.
- **Fig 3:** Hallmark EMT 194개 중 방향 일치: 48h 두 연구 141개, 500nM 48h↔500nM 7d 156개, 100nM 48h↔500nM 7d 145개, 세 실험 모두 124개. 탐색적·상관된 유전자 기반 비율이지 경로 검정 p값이 아니다.

[검증한 논문 Figure와 Results 해석 문서](THESIS_FIGURES_AND_MANUSCRIPT_RESULTS_VERIFIED_20261031.md)를 본 연구의 공식 작성 기준으로 사용한다. 다음은 figure 생성 반복이 아니라 (1) 원 논문 80개 유전자 원본 목록 확보 (2) 기존 논문의 핵심 발견 대비 novelty table (3) 결과 Methods·Discussion 통합이다. **아직 wet-lab ECM 기능, 임상 스테로이드 ∆IOP, 약물 인과 효과는 입증되지 않았다.**

## Gate B 공통 donor 기준 검증 완료 — 55명 BeamA/scPred raw UMI 이중 재분석

2026-10-31 [공통 55명 두 주석 원시 UMI/edgeR run #41206289766](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41206289766)에서 **18개 원본 한 번 집계 및 두 edgeR paired model 적합이 실제 성공**했다. 다만 해당 run의 전체 결론은 **FAIL**: 마지막 비교 스크립트가 scPred에서 발현 필터로 제거된 `ERFE`를 반드시 양쪽에 있어야 하는 유전자로 간주했기 때문이다. 유전자별 결과를 다운로드하여 ERFE는 *untestable*로 표시하는 수정 후, **[집계된 결과만 다시 비교하는 CI #41207583296](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296) SUCCESS**, [artifact #13045350639](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296/artifacts/13045350639) 실제 TSV/JSON 점검 완료. [최종 검증 문서](GATE_B_MATCHED55_BEAMA_TM_FIBRO_ACTUAL_VERIFIED_20261031.md)에 전체 근거와 실패·수정 기록 보존.

- 실제 유효한 **동일 55 biological donor pairs**를 BeamA score≥0.75와 scPred TM_Fibro 각각 분석 (55 vs 55로 donor가 다른 비교 아님).
- BeamA **16,051개**, scPred TM_Fibro **14,003개** 유전자 검정; 공통 **13,975개**. 두 효과 log2FC 상관 **r=0.94990493**, 방향 **12,634/13,975**, 두 모델 모두 FDR<0.05 **10,255개**.
- **MYOC (+4.006/+3.525), TGFBI (−1.868/−1.941), MMP2 (−0.591/−0.459)** 모두 같은 방향 및 두 모델 FDR<0.05로 computational **GO**. ACTA2는 두 주석에서 +/+로 안정적이나 기존 7일 독립 TM 연구에선 −라서 *3-study universal candidate*가 아니다.
- 16개 후보 가운데 **15개 검정 가능하고 15개 방향·FDR 모두 일치**; ERFE는 scPred TM_Fibro에서 `filterByExpr`로 제외, FDR이 없으므로 비유의 또는 반대라고 해석하지 않는다.
- 55명 2-arm×2-label cell 합계와 source UMI QC 확인, 원시 donor count 결과는 GitHub runner 내에서만 사용 후 제거하고 artifact에는 공개 aggregate만 보존.

**논문의 새로운 강점:** 3개 별도 연구에서 MYOC/TGFBI/MMP2 발현 방향을 확인했고, 내부 동일 원본 57명 confidence 재적합 및 55명 고정 alternate-label 재분석까지 수행했다. 다만 여러 cell-label 분석은 **동일 donor·세포를 공유**하므로 독립적 복제 연구 2개를 추가한 것이 아니다. 처리 후 cell classification 및 DEX/10x library 준비 완전 교란은 여전히 남는다. 이후 작업은 추가 대형 workflow 확장보다 source provenance Figure·3-study comparison·독립 기능/임상 증거 계획 정리에 우선 집중한다.

## Gate B 최신 — full raw UMI cutoff 검증 완료 / 대체 주석 확장 결정

[완료된 0.50 재집계 + 18개 10x gzip QC + 57 donor edgeR run #41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259) **SUCCESS**, 원시파일 SHA/CRC와 실제 분석결과 [보고서](GATE_B_BEAMA_050_075_ACTUAL_RAW_EDGER_VERIFIED_20261031.md) 확인. 같은 57명, 16,050 공통 유전자 효과 상관 **r=0.99911597**, Δ|log2FC| 중앙값 **0.0112462**, MYOC/TGFBI/MMP2 **3/3 방향/FDR 안정**, ACTA2도 동일 연구 내 기준 간 안정이나 독립 7일 TM 연구와 방향 불일치. 고효과 내부 유전자 |log2FC|≥3은 0.75 기준 **83**, 0.50 기준 **88**, **81** 교집합. 기존 83의 2개 ACP5/HLF는 threshold 경계 이탈, 신규 7개 중 5개는 이전 0.75에서 검정되지 않았다. **원 저자 63 증가·17 감소 80개의 원본 명단 비교가 아니라 우리 내부 분석의 cutoff 민감도 결과**이다.

[대체 세포주석 metadata run #41201790738](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41201790738)도 SUCCESS: scPred **TM_Fibro 조건별 ≥30세포 paired donor 57명**, 9 pools, 이 중 **55명** 원래 BeamA57과 공통. JCT any 56명, JCT confidence≥0.75는 **8명**에 불과. [실제 주석 준비상태 보고서](GATE_B_ALTERNATE_CELL_STATE_METADATA_20261031.md)를 참고. **다음 우선 순위는 단일 패스 scPred TM_Fibro/BeamA 공동 55 donor raw-count 분석 타당성 및 구현**; 이전 BeamA .50/.75 원시 매트릭스 작업을 기계적으로 반복하지 않는다. 해당 대체 주석 유전자 발현/FDR은 아직 분석하지 않았으므로 HOLD.

중요한 한계: 효과 수치가 0.50/0.75에서 일치해도 DEX와 다른 10x 처리 라이브러리의 완전한 교란, 처리 후 레이블 선택과 독립 임상 안압 결과 부재는 해결되지 않았다. **논문 주장: 3개 독립 연구의 전사 반응 후보 우선순위 + 내부 cutoff 민감도**, 기전/인과/임상 예측은 불가. 원 저자 고효과 80개 개별 원본 유전자 명단 역시 확보 전까지 HOLD.

## Gate A 업데이트 — 제3 DEX–TM 코호트 실제 확인 완료 (2026-10-31)

이전 마스터 문서에서 최우선 P0로 지정했던 GSE301517의 원본·메타데이터 검증은 이제 **완료**되었다. [최종 3연구 비교 Actions #41186736618](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186736618) PASS, 공개 artifact #13037335888 검증. 자세한 결과는 [GSE301517_THIRD_COHORT_7DAY_ACTUAL_VERIFIED_20261031.md](GSE301517_THIRD_COHORT_7DAY_ACTUAL_VERIFIED_20261031.md)를 기준으로 사용한다.

중요한 새로운 사실: 원 연구 GEO 처리는 **500 nM DEX, 7일**이고 기존 GSE313736 **500 nM, 48시간**, Mehrotra **100 nM, 48시간**이므로 셋을 시간/농도의 동일 조건 재현성으로 취급하지 않는다. 공식 CPM 7,905,333 bytes SHA256 `717f4c9861de4d6b169439ef56a781f50d93ec918c456a9d0e61873e4a6e5f09`; 21 GSM / 85,298 Ensembl rows / 14 Control-Dex 발현 샘플 / TM62 4+4, TM70 2+2, TM72 1+1 / 별도 Resource 7건. 3개 TM strain에서 각 처리군 반복을 먼저 평균 내므로 독립 donor를 14명으로 과대계산하지 않았다. 57-donor 원 연구와 15,876 Ensembl ID 매칭.

**세 연구에 걸쳐 방향이 모두 일치한 핵심은 MYOC↑, TGFBI↓, MMP2↓ 총 3개**다. ACTA2는 두 48h 연구에서 ↑, 새 7일 연구에서 ↓로 **핵심에서 조건 의존성 후보로 내려야 한다**. EMT shared 194 중 156개가 기존 500 nM 48h와 일치, **124/194개가 세 연구 모두 같은 방향**. P9 핵심 4개를 모두 일반화된 안정 핵심이라고 설명한 기존 결론은 이 후속 결과로 **수정한다**. MYOC/TGFBI/MMP2의 7일 CPM +1 log2 비율은 원자료 기반 edgeR log2FC와 같은 모수로 비교하거나 임상 의미를 단정하면 안 된다.

**다음 실제 최우선은 Gate B:** BeamA post-exposure selection (score .5 vs .75) / cell-state annotation/배치 민감도 및 원 논문 63↑/17↓ 정확 유전자 목록 비교. Gate A를 반복해서 '확인 예정'이라고 쓰지 않는다. GSE301516 ATAC/POAG GWAS 연동은 이 기반을 정리한 후 조건부 검토. 새로운 P번호를 계속 늘리지 않는다.

## Gate B 실행 갱신 — 실제 메타데이터 PASS / 원시 count 재적합 미완료 (2026-10-31)

[공식 GEO 메타데이터 Gate B run #41198576725](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198576725) PASS. BeamA confidence 0.75와 0.50을 비교하면, **각 조건에 최소 30세포인 paired donor는 57명으로 동일**하다. score0.5에서 포함되는 메타데이터 세포는 Vehicle **48,903**(0.75에서는 40,115), DEX **55,652**(0.75에서는 46,454)이고, donor 추가는 **0명**이다. 이 수치는 10x 원본과 실제 barcode join 전 메타데이터 개수이므로 원시 count 기반 유전자 효과 재현성이라고 해석하지 않는다.

원 저자 GEO 요약은 63개 증가·17개 감소라고 보고했으나 GEO Series 공식 보충자료 목록에서 그 80개 Ensembl 유전자 명단은 제공되지 않아 우리 고효과 83개(67 증가·16 감소)와 정확한 목록 비교는 **HOLD**다.

원시 score0.5 재집계는 공식 GEO 다운로드가 간헐적으로 잘못된 gzip·중도 종료 응답을 반환하여 [첫 실행 #41198827259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198827259)와 [재실행 #41198957708](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198957708)이 실패했다. 전송 재시도 및 대형 행렬 전체 CRC를 추가해 [검증 재실행 #41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259)를 시작했다. **실제 최종 결과 확인 전에는 edgeR score0.5의 효과크기·FDR에 어떤 결론도 부여하지 않는다.** 상세 문서는 [GATE_B_BEAMA_SELECTION_AND_AUTHOR_DEG_AUDIT_20261031.md](GATE_B_BEAMA_SELECTION_AND_AUTHOR_DEG_AUDIT_20261031.md).

## 현재 확인된 상태

GitHub [PR #253](https://github.com/kimtk94/IS_Analysis_V3/pull/253)은 검토 시점 OPEN, mergeable, 미병합. 기준 커밋 `1a7d2a61e`, 변경 98파일, 159커밋. 최신 [경량 CI #41185786190](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185786190) 성공, [P9 최종 3개 검증 #41185708990](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185708990) 성공. 구현과 CI의 성공은 생물학적 검증 완료가 아님.

| 부분 | 실제 성과 | 해석 |
|---|---|---|
| 원자료 | GSE313736 18/18 10x 원본 무결성 확인 | GO |
| paired 분석 | 57 donor, 114 pseudobulk, 9 pool | GO |
| edgeR | 16,050 검정, 13,717 FDR 유의 (85.46%) | 통계 실행 PASS, 기술적 교란 주의 |
| 독립 비교 | Mehrotra 2026 100 nM TM 보충표: 12,251 ID 매칭; 824개 외부 유의 유전자 중 659개 같은 방향·양쪽 유의 | 외부 연구의 발현 방향 재현성 GO |
| EMT | 194개 중 141개 부호 일치; 49개 양쪽 유의, 그중 46개 일치 | 탐색적 GO, 기능/인과 HOLD |
| P9 패널 | MYOC/TGFBI/MMP2/ACTA2 등 16개 분류 | 실험 설계 GO, wet-lab 미실시 |
| 미생물 | 16S 및 shotgun 접근성 검증 | 환자군/대조군 매핑, 오염 대조 HOLD |
| 임상 | 임상 steroid ΔIOP 결과 없음 | 예측 모델 STOP |

### 근본적 설계 한계

DEX와 vehicle은 각 pool 내 서로 다른 10x 라이브러리에서 처리되었고, 분석하는 BeamA 고신뢰도 세포는 노출 **후** 분류되어 선택되었다. Donor blocking, 9-pool leave-out, 외부 100 nM 연구와의 일치성은 유용하지만 약물 인과효과를 완전히 분리하지 못한다. 85% 이상의 대규모 DEG 비율은 반드시 기술적·생물학적으로 설명해야 한다.

원본 보고의 74명/316,799 singlets와 우리 GEO GSM 70명, 최종 metadata 59명, 엄격 BeamA 57명의 차이와 세포 선택 편향을 투명하게 Figure로 보고한다. EMT 유전자 전사량 감소가 곧 ECM 축적 증가·유출저항 증가를 의미하는 것은 아니다.

## 선행연구와 신규성 제한

- GSE313736 원 연구가 이미 74 donor 기반 sc-eQTL 및 DEX 발현 변화를 보고하므로, 원 저자의 63 증가/17 감소 고효과 유전자와 **유전자 ID 단위 비교**가 필요하다.
- Mehrotra et al., IOVS 2026 (PMID 41636427): 9개 QC 통과 TM strain, 100nM 1/6/48h, 48h 857 high-confidence DEG와 ECM/adhesion 경로를 이미 보고.
- Singh et al. 2026 (PMID 41890121): DEX–TM chromatin·GWAS 통합으로 26 IOP + 52 POAG 후보 유전자 보고. 단순 ATAC/GWAS overlap은 신규성이 아니다.
- 2026 native human glaucoma TM atlas (PMID 42396508): 건강·녹내장 각각 7명, 17개 cell population, ECM/염증/스트레스 서술. 세포별 원시 데이터 공개 다운로드 여부는 별도 검증 전 HOLD.
- 실제 steroid 이후 ΔIOP GWAS (PMID 42388255)는 별도 임상 phenotype. POAG나 baseline IOP와 같은 결과로 취급하면 안 되며 접근 승인 데이터와 혼동하지 않는다.

## 수정된 단일 연구 질문

**인간 섬유주의 스테로이드 노출 연관 반응 중 세포상태 및 실험 코호트에 관계없이 반복되는 전사 프로그램은 무엇이며, 이 중 녹내장/안압 유전학 및 실제 TM 조직에 연결될 수 있는 후보는 무엇인가?**

영문 가제: *Cross-study and cell-state-aware prioritization of glucocorticoid-associated trabecular meshwork programs with glaucoma relevance.*

현재는 `association` 및 `candidate prioritization` 연구. 새로운 원인 유전자, 용량별 반전, 단백질 기능, 환자별 위험 예측이 검증된 연구가 아님.

## 우선순위 재편 — 세 가지 핵심 Gate

### Gate A. 제3의 독립 TM DEX 연구 확인 — 최우선

[GSE301517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517)은 21개 RNA-seq GSM 중 제목상 Resource 7, Control 7, Dex 7개이며 TM62/TM70/TM72의 3개 서로 다른 세포주에 Control/Dex 기록이 있다. **14개 대조/처리 파일은 14명 독립 donor가 아니다.** 공개 처리 결과 `GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz`는 약 7.5MB. 원본 제목, donor/technical 반복 관계, dose/time, 컬럼 ID, 파일 완전성 및 유전자 매핑을 실제로 검증한 뒤 진입한다.

CPM은 정규화된 수치이므로 raw integer count로 간주하여 edgeR/DESeq2에 넣지 않는다. 원시 count가 확보되지 않으면 후보군 발현 방향 검증만 수행하고 통계적으로 독립 재분석했다고 쓰지 않는다.

동일 SuperSeries의 [GSE301516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301516)에는 ATAC 데이터 12개 실행과 공개 peak 파일이 있어 2차 chromatin 검증 대상이지만, 기존 Singh 논문의 기발표 유전자를 신규로 주장하지 않는다.

**GO 기준:** 정확한 생물학적 복제 단위가 ≥3 strain으로 입증되고, 처리/대조 3쌍 및 정량 파일의 유전자와 컬럼 매핑 성공. 적절한 가설과 독립 비교 결과를 작성할 수 있어야 한다. **HOLD:** strain/replicate/time/dose 미확인 또는 처리 파일이 확인 불가능한 경우.

### Gate B. 기존 57 paired donor 분석의 생물학적 타당성

BeamA confidence 0.5/0.75, alternative cell label, donor당 적격 cell 수, 각 조건의 donor N 변화, pool 간 편향 및 이상치를 비교한다. 원 연구의 63↑/17↓ 명단과 우리 83개 고효과 DEG를 정확한 gene ID로 대조한다. MDS/PCA와 동일 발현량 범위에서 일치성을 검토한다. 완전한 treatment/library 교란을 통계적으로 제거했다고 쓰지 않는다.

**GO 기준:** 사전 설정한 핵심 패널·경로의 부호 안정성, 공통 유전자 universe, donor N, 정확한 QC 및 원 논문 대비 차이와 제한을 보고. **HOLD:** 결과가 post-treatment classifier cutoff 또는 특정 pool에 민감하거나 복제 단위가 불명확한 경우.

### Gate C. 논문 생물학·유전학적 연결 — 조건부

확인된 독립 RNA 결과에 연결되는 candidate만 disease-related ATAC/IOP/POAG/public native TM tissue와 교차검증한다. eQTL/MR/coloc은 SNP effect, allele harmonization, matched LD, full summary statistics와 충분한 도구 변수가 확인될 때만 진행. **특히 steroid-induced ΔIOP와 일반 POAG/IOP GWAS를 동일 outcome으로 혼동하지 않는다.**

실제 0/100/500nM 공여자별 paired 48h 실험, RT-qPCR MYOC/TGFBI/MMP2/ACTA2, MYOC protein, MMP2 activity 및 ECM 기능 검증은 협력 실험이 가능할 때만 별도 의사결정한다. 표본 수는 pilot donor variance 기반으로 산정하고 무작위 배치/용매매칭/오염 QC 기록. **실험이 없다면 생물학적 기전이 아닌 계산적 재현성과 후보 우선순위 논문으로 범위를 낮춘다.**

## 실행·문서 정리 기준

1. **더 이상의 개별 Pxx 워크플로 난립 금지**. 기존 P0–P9 실행 ID·산출물은 그대로 보존하며 본 문서를 단일 기준점으로 사용.
2. 159커밋/98파일 PR은 지나치게 거대하므로 병합 전에 실행 증거 인덱스·소스 해시·필수 재현 진입점·이전 실패와 수정 내역·GitHub 경량 CI 확인. 필요하면 향후 변경을 작은 독립 PR로 분리. 이 문서는 병합을 실행하지 않는다.
3. 외부 독립 DEX-TM 자료를 우선하고 미생물 앰플리콘/shotgun은 **별도 녹내장 프로젝트로 보류**. 오염대조/환자군 매핑 없이 미생물 기전과 스테로이드 안압 상승을 연결하지 않는다.
4. 연구 그림: (i) 데이터 계보/선택 편향 (ii) QC+모형 (iii) 3개 독립 실험의 동일 유전자 효과 방향 (iv) cell-state threshold sensitivity (v) ATAC/IOP 조건부 통합 (vi) 근거 수준별 GO/HOLD 표.
5. **논문의 완료 조건**: 제3 코호트 GO, 원논문 완전 DEG 매핑, 선택/교란 한계의 정량화, 결과 해석을 뒷받침하는 그림, 선행 2026 문헌 대비 비중복 기여, 외부 원자료/실험 접근 여부에 맞는 주장 수준. 임상 예측과 novel causal mechanism은 현재 STOP/HOLD.

*전략 수립 시 출처를 확인한 참고 논문 및 데이터: GSE313736; GSE301517; GSE301516; PMID 41636427, 41890121, 42396508, 42388255; dbGaP phs004037. GitHub CI 통과는 생물학적 결론 PASS와 다르다.*
