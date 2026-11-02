# GSE313736/TM DEX 석사학위논문 보완 우선순위 (2026-11-02)

**기준:** 이미 계산한 통계와 산출물은 GO. 학위논문 제출과 심사 통과는 미확정. 필수 6개, 선택 3개. 원저자 DEG80 유전자 ID 요청 및 연구진 연락은 비필수 HOLD.

## 실제 검증 완료

- 57 paired donor edgeR, 16,050 tested, 13,717 BH-FDR<0.05, 내부 기준 83 high-effect.
- 동일 57 donor BeamA confidence threshold comparison; matched 55 donor BeamA/scPred; 발현량 및 효과크기별 방향 민감도.
- Mehrotra 100nM 48h, GSE301517 500nM 7d 세 TM strain; MYOC ↑, TGFBI ↓, MMP2 ↓의 방향 일치 및 EMT aggregate 124/194 versus strict all-three-strain 90/194.
- 9 original source pool x 2 original 10x libraries, 원본 library+treatment 19 columns/rank18; pooled BeamA delta +2.406339 pp = +0.868638 pp within-pool +1.537701 pp weights. 교란 원인별 인과 추정은 불가능.
- 메인 Figure 1–3, 실제 시각 검수를 거쳐 수정한 Figure S1, GitHub CI와 Google Drive 보관 완료.

## 필수 미완료 6개: 구체적인 완료 조건

| ID | 우선 | 보완 항목 | 완료 산출물 | 주의 및 판정 |
|---|---|---|---|---|
| P0-A | 최상 | 지도교수의 방법론 중심 석사논문 주제 승인 및 소속 디지털헬스 과정의 학위심사 자격·등록/학점/예심/심사 일정 확인 | 승인 여부 및 과정별 적용 규정 | 연구 CI PASS로 승인/수료/학위신청 자격을 추정하지 않는다 |
| P0-B | 최상 | 74 original source persons, 59 author-QC persons, 자체 paired 57명, fixed same-annotation 55명: 표본과 QC 기준의 서로 다른 분모 명시 | 표본별 출처·포함/제외 기준 Table 1, Figure 0 donor/cell/pool inclusion flow, source GSM 시간충돌 표시 | 독립 ID 연결 없이 59명에서 정확히 2명이 빠졌다고 단정 금지 |
| P0-C | 최상 | 9개 pool x 2개 처리군 x 18 source 10x library의 관계를 시각적으로 설명 | 별도 schematic Figure, paired donor pseudobulk 단위 대 9 pool unit 비교, rank19/18 증명 캡션 | DEX 독립 계수 식별 불가. 원본 약물 효과가 0이라는 뜻은 아님 |
| P0-D | 최상 | 13,717/16,050 DEG와 83 내부 high-effect 후보의 FDR/효과크기 절단값 및 사후 선택성 분리 | 유전자 검정 universe, FDR, abs(log2FC) 경계, 유의성·탐색·방향 검증 차이 명시한 Table 2와 원고 검토 | GSE301517 7일 CPM ratio를 edgeR log2FC나 독립 gene-wise FDR로 취급 금지 |
| P1-E | 높음 | 선행연구 2022·2023·2026 IOVS Mehrotra·2026 Singh preprint 대비 새로운 기여 수준 정리 | 논문별 prior observation과 우리 analytical novelty를 1:1 매칭한 Table 3 + Intro/Discussion 반영 | MYOC, ECM/EMT 발견 자체를 신규 발견으로 주장 금지 |
| P1-F | 높음 | 현재 Markdown manuscript를 대학원 정식 논문 DOCX/PDF로 완전 조립·교정 | 국문/영문 제목·초록, 목차, 방법, 결과, 그림 1–3, **수정 Figure S1**, 표, 실제 DOI/참고문헌, 코드/자료 출처, 부록, PDF 페이지 단위 검수 | 심사용 최종본 내용 및 학과 양식·서명/전자/인쇄 제출 절차 확인 전 제출 준비 완료 주장 금지 |

## 선택 연구 보강 3개 — 학위논문 완성을 위한 필수조건 아님

| ID | 목적 | 가장 적은 계산으로 할 일 | 제한 |
|---|---|---|---|
| O1 | 83 high-effect 기준 사후성 민감도 보강 | 기존 edgeR 통계 아티팩트에서 abs(log2FC)≥1·2·3과 FDR<0.05을 사전에 고정해 후보 수·핵심 3유전자 방향 표 작성 | 절단값 이동으로 약물 인과 검증을 했다고 주장 금지 |
| O2 | 독립 Mehrotra 연구의 실질적 교차 재현성 향상 | 2026 정정공고에 기재된 SRA PRJNA1243143의 1h/6h/48h, TM/SCE, strain/arm accession-level 원본 metadata crosswalk; 결과가 맞을 때에만 실제 raw-count 재분석 고려 | 기존 논문 FDR와 7일 CPM sign을 동일 추정량으로 합치지 말 것 |
| O3 | 향후 기전/임상 결론 강화 | 별도 독립 코호트 또는 라이브러리/처리군 교차 균형된 실험과 단백질·IOP 평가 설계 | 석사 computational reproducibility 논문에서는 미래 연구 제안으로 충분 |

## 중요도 구분

- **논문 작성 즉시 P0:** P0-B (표본 선정/흐름도), P0-C (교란·실험단위 구조도), P0-D (통계적 주장/유전자 기준 명확화), P1-E (선행연구 차별성).
- **제출 전 P0:** P0-A (학과/지도교수), P1-F (완성형 논문/교정).
- **연구 자체 진행에는 nonblocking:** 원저자 80개 유전자 ID, 원저자 연락, 임상 steroid-induced ∆IOP, 신규 wetlab, 제3 연구 새 유전자별 p/FDR. 이를 논문에서 미확인이라고 공개한다.

## 학사 일정 주의 (2026-11-02 기준)

2026-2 일반대학원 *일부 학과* 공개 공지에서 논문지도/예심/심사 신청 마감은 2026-10-26이며 이미 지났다. 이 공지가 디지털헬스 과정에도 동일 적용된다는 확증은 없다. 당학기 심사 희망 시 반드시 소속 과정의 GLS 등록/신청 상태와 사무실 공지 확인. 학위 심사 eligibility는 분석 완료와 별개.

학교 공개 참고 (과별 차이 주의): https://soc.skku.edu/soc/notice.do?articleNo=226550&mode=view

## GitHub Actions 권장 순서

1. 병렬 A: sample membership/source-unit flow 및 9pool/18library 설계 그림 (표본 crosswalk 없는 부분은 unknown 화살표).
2. 병렬 B: cutoffs/selection/FDR 수치표 및 선행연구 기여 대조표를 문서·CI로 자동 생성.
3. 마지막: 학과 논문 템플릿 기반 DOCX/PDF 조립, 보조그림은 corrected Fig S1만 사용, 최종 PDF 시각/참조 검수.
4. 학사 단계: 지도교수 연구 방향 승인 및 학사 자격 확인.

관련 검증본: [현행 논문 초안](THESIS_MANUSCRIPT_METHODS_RESULTS_DISCUSSION_DRAFT_20261101.md), [최종 과학적 판단](MASTER_THESIS_FINAL_SCIENTIFIC_GO_NO_GO_20261102.md), [그림 및 출처 인덱스](THESIS_PUBLICATION_READINESS_EVIDENCE_INDEX_20261102.md).
