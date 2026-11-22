# IS P3A 실제 36개 그룹 / 20개 유전자 일괄 복원·검증 — 2026-11-22

**연구 범위:** MASTER DEGREE 허혈성 뇌졸중. Sales 목적 연구 제외.

**검증 결과:** 과거 [학술 연구 브랜치의 원본 커밋 0bf646f](https://github.com/kimtk94/IS_Analysis_V3/tree/0bf646fd6f29322d17da878253f37def24c01c74/docs/data/is)에 보존된 393개 gene×QTL×region 원본 행에서 `possible_single_signal_ABF_followup=True` **36개**를 정확한 자연키 `region_id × QTL_dataset × gene_id`로 복원했다. 36개는 서로 다른 **20개 Gene ID**이며 모두 `IS_XDATA_G0015` 영역에 속한다. 검증된 QTL별 수는 `QTD000131=10`, `QTD000171=11`, `QTD000434=15`다.

당시 서버 전용 상세 `science/archived_fasta_abf_audit.json`의 **25개 source-stat-positive / 11개 source-stat-missing 개별 배정은 아직 복원되지 않았다.** 이 숫자와 독립적으로, 실제 보존된 **98개 공식 BH-FDR eGene 표에서 36개 중 25개**가 발견되고 11개는 발견되지 않았다. 두 분할이 숫자가 같더라도 **동일한 그룹 구성이라고 가정하면 안 된다.**

## 20개 유전자 — 모든 36개 검증 그룹 포함 (인과성 순위 아님)

원본 QTL BH q는 해당 molecular-eQTL 검정의 증거일 뿐 stroke-coloc/MR/독립 재현 확률이 아니다. 표는 유전자명 알파벳순 정렬이다.

| Gene | QTL 분석 그룹 수 | eGene BH q <0.05 그룹 | 가장 작은 원본 QTL BH q | QTL dataset suffix |
|---|---:|---:|---:|---|
| AS3MT | 1 | 1 | 1.09e-29 | 434 |
| C10orf32 | 3 | 3 | 8.69e-59 | 131, 171, 434 |
| CALHM2 | 3 | 3 | 2.26e-3 | 131, 171, 434 |
| CCDC147 | 1 | 1 | 2.74e-2 | 434 |
| CNNM2 | 1 | 1 | 1.48e-5 | 434 |
| CYP17A1 | 1 | 0 | — | 171 |
| GSTO1 | 3 | 2 | 5.85e-8 | 131, 171, 434 |
| GSTO2 | 3 | 3 | 2.54e-57 | 131, 171, 434 |
| INA | 1 | 1 | 1.57e-3 | 434 |
| ITPRIP | 2 | 1 | 1.21e-2 | 171, 434 |
| NEURL | 2 | 1 | 1.03e-4 | 131, 434 |
| NT5C2 | 2 | 0 | — | 131, 171 |
| PDCD11 | 1 | 1 | 1.75e-2 | 131 |
| RP11-127L20.3 | 2 | 2 | 1.03e-11 | 171, 434 |
| RP11-416N2.4 | 2 | 1 | 1.73e-3 | 131, 171 |
| RP11-724N1.1 | 2 | 1 | 2.48e-18 | 171, 434 |
| SFR1 | 1 | 0 | — | 131 |
| SLK | 2 | 2 | 3.85e-5 | 131, 434 |
| USMG5 | 1 | 0 | — | 434 |
| WDR96 | 2 | 1 | 1.04e-10 | 171, 434 |

- 36개 그룹의 원본 식별자와 공식 eGene 매칭값: [복원된 실제 group manifest](evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json)
- 원본 393×98×393개 검증 테이블 **전체를 CI에서 다시 읽고 Git blob SHA1와 행별 값을 교차 대조**: [Actions #45019559971](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45019559971), 12/12 tests PASS.
- 과거 원본 GWAS 두 릴리스의 publisher MD5/SHA256/bytes 및 eQTL source checks: [Actions #45018300939](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45018300939), technical PASS; scientific HOLD.
- 이미 검증된 원본 GWAS artifact만 재사용해 **36개 실제 group ID를 통합한 경량 evidence audit**: [Actions #45019782997](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45019782997), 13/13 tests PASS.

## 임상·통계적 해석 제한

1. 36개 모두 원본 데이터 유래 ID지만, source-stat 25/11의 개별 배정은 서버 원본 미확보로 **UNVERIFIED**. 그룹별 QTL/GWAS beta/SE/MAF/효과 allele, trait case fraction·표본 수, 완전한 cis 통계량의 **동일 원본 출처**를 아직 대조해야 한다.
2. 36개 전부 molecular-QTL donor ancestry 및 genotype-matched **signed LD 미승인**, 수집된 eGene BH05 수치만으로 ABF/SuSiE를 실행할 수 없다.
3. EAS GIGASTROKE 발견 코호트에는 BBJ가 포함될 수 있어, BBJ 대비 독립 검증이라는 주장은 승인되지 않는다.
4. 현재 신규 P3A 연구의 **공식 Coloc 실행 승인 0/36, 신규 ABF·SuSiE 실행 0/36**. "36개 부정"이 아니라 **36개 아직 검증 필요**.
5. 이전의 FGF5·CALHM2 등 **15개 BBJ 레거시 탐색적 ABF 그룹은 이번 36개 전부가 속한 chr10 영역과 다른 입력 모델**이므로 합쳐서 H4 순위를 비교하면 안 된다.

## 후속 연구 계획

**배치 처리 단위:** 36개 고정된 자연키 × 원본 GWAS·QTL source 검증 → 그룹별 `HOLD/ELIGIBLE` → `ELIGIBLE`에 한해 정확한 원본 cis window에서 `coloc.abf`와 사전확률 민감도를 GitHub Actions 매트릭스로 계산 → LD 기준을 충족한 그룹만 `coloc.susie` → 전체 36개 결과와 20개 유전자 요약을 통합한다. 분석 대상에서 제외되는 그룹도 사유와 함께 남긴다.

**다음에 필요한 외부 원본:** 과거 서버의 36개 그룹별 source-stat-positive membership, 원본 GWAS effect-allele signed beta/SE·N/case fraction의 archived row 동등성 증명, QTL donor N·영역별 완전한 nominal cis 통계 및 QTL-matched LD. 위 요건 중 하나라도 충족되지 않는 그룹은 추가 p-value/H4를 임의 생성하지 않는다.
