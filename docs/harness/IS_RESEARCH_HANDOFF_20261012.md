# IS Research — Harness 인수인계 및 완료 계약
**기준:** 2026-10-12 KST | **Project:** MASTER DEGREE / Ischemic Stroke | **상태:** 공개자료 기반 탐색적 분석 완료, 정식 인과 검증 미완료

## 1. 목표 및 현재 위치
일본인·동아시아 허혈성 뇌졸중(IS) GWAS(BBJ, GIGASTROKE)와 뇌 조직 eQTL(EBI QTL Catalogue GTEx v8 Brain Cerebellar Hemisphere QTD000161)을 통합해 재현 가능한 석사 논문용 후보 유전자·좌위 분석 파이프라인을 완성한다. 연구 구조는 GWAS → 정밀 대립유전자 매칭 → LD/conditional GWAS → fine-mapping → multi-signal colocalization/MR → 조직·세포 기능 주석 → 논문 Figure/Methods다. **정식 통계 방법이 실행되지 않았다면 인과 후보를 확정하지 않는다.**

**현재 최우선 문제:** chr10의 CALHM2 eQTL 신호가 BBJ IS lead와 동일 신호에 속하는지, GIGASTROKE의 낮은 LD lead와 독립적인지 검증. 네 유전자(FGF5/C4orf22/CALHM2/NEURL)는 집중검증 하위집합일 뿐이며 **전체 IS discovery 2,225개 유전자 범위 유지**가 필수다.

## 2. 환경 및 서버 경로
- 서버: `taehoonserver`, 서버 루트 `/srv/is-analysis/`.
- 원본 저장소: `/srv/is-analysis/IS_Analysis_V3`.
- **현재 작업 트리:** `/srv/is-analysis/worktrees/is-broad-discovery-20261009`; GitHub `kimtk94/IS_Analysis_V3`, 브랜치 `research/is-broad-discovery-20261009`; 기준 HEAD **`9c1a641`**.
- BBJ canonical GRCh37 전체 GWAS(13,435,531행): `/srv/is-analysis/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz`.
- GIGASTROKE EAS AIS canonical GRCh37(6,725,216행): `/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz`.
- 두 GWAS × EBI 분자 QTL 통합: `/srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv` (20,689 gene×SNP).
- BBJ/EBI 검증 입력: `/srv/is-analysis/results/is/audits/is_priority4_ebi_bbj_abf_pilot_20261011_v2/IS_PRIORITY4_EBI_BBJ_ABF_QC_VALIDATED_INPUT.tsv` (11,050 gene×SNP).
- EAS504 VCF 및 PLINK RAW: `/srv/is-analysis/results/is/audits/is_priority4_eur503_eas504_conditional_20261011_v1/`.
- EUR503 PLINK RAW: `/srv/is-analysis/results/is/audits/is_priority4_eur503_fullcis_signed_ld_qc_20261011_v1/`.
- JPT104 식별용 공식 1000G sample metadata: `/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel` (JPT 104명은 EAS504에 포함; 그 외 EAS400명).
- JPT/EAS/EUR candidate 부트스트랩: `/srv/is-analysis/results/is/audits/is_priority4_jpt_bootstrap_20261011_v1/`.
- BBJ COJO 사전 QC: `/srv/is-analysis/results/is/audits/is_priority4_bbj_eas_cojo_readiness_20261011_v1/`, `is_priority4_bbj_eas_cojo_reference_20261011_v1/`, `is_priority4_bbj_eas_cojo_final_20261011_v1/` (정식 COJO 결과 아님).
- 제한 없는 full GWAS lead: `/srv/is-analysis/results/is/audits/is_priority4_unrestricted_gwas_leads_20261011_v1/`.
- 두 GWAS 교차 p·대립유전자 검증: `/srv/is-analysis/results/is/audits/is_priority4_unrestricted_cross_alleles_20261011_v1/`.
- chr4 palindrome LD: `/srv/is-analysis/results/is/audits/is_chr4_palindrome_lead_ld_20261011_v1/`.
- chr4/chr10 GWAS–QTL multi-track Figure: `/srv/is-analysis/results/is/audits/is_priority4_locus_multitrack_20261011_v1/` 및 Git `docs/figures/IS_BBJ_IS_L001_BBJ_GIGA_EBI_QTL_MULTITRACK_20261011.svg`, `IS_BBJ_IS_L002_BBJ_GIGA_EBI_QTL_MULTITRACK_20261011.svg`.
- **최신 chr10 두 lead LD/분자 QTL 감사:** `/srv/is-analysis/results/is/audits/is_chr10_two_lead_ld_tags_20261011_v1/` — `IS_CHR10_EAS_JPT_TWO_LEAD_PER_SNP_LD.tsv`, `IS_CHR10_EAS_JPT_TWO_LEAD_LD_QC.tsv`, `IS_CHR10_EAS_JPT_TWO_LEAD_LD_MANIFEST.json`.
- 상세 방법·해석은 작업 트리 `docs/IS_PRIORITY4_*.md`, `docs/IS_CHR4_*.md`, `docs/IS_CHR10_*.md` 참조.

## 3. 검증된 연구 결과 및 해석 수준

### chr4 FGF5/C4orf22
BBJ 전체 좌위 lead GRCh37 `4:81182554:T:C`, p=1.2119e-8. GIGASTROKE EAS AIS 전체 좌위 lead `4:81181072:T:A`, p=2.583e-8. **두 위치 모두 양 GWAS에서 p<5e-8**. 두 lead 위치 간 공개 참조 LD: **EAS504 r²=.942015, JPT104 r²=.974855**, 같은 고LD 지역 신호를 태깅할 가능성이 큼. 다만 GIGASTROKE lead는 A/T palindrome이므로 가닥 provenance는 **UNRESOLVED**. 단일 인과 변이 또는 FGF5 인과 유전자 확정 불가. 가정적 eQTL effective N=175, p12=1e-5에서 단일 신호 ABF PP.H4 약 **FGF5 .913 / C4orf22 .421**; 공식 eQTL N은 검증되지 않음.

### chr10 CALHM2/NEURL
BBJ lead `10:105459834:T:C`, p=8.762e-10. GIGASTROKE lead 위치 `10:105599770` (원본 ID `G:T`, 참조 `T:G`), p=5.584e-11. 상호 교차 p: GIGASTROKE에서 BBJ lead **7.792e-7**, BBJ에서 GIGASTROKE lead **9.246e-8**; 효과 방향은 harmonization 후 일치하지만 동일 SNP의 양쪽 genome-wide 유의성 동시 달성은 아님. 두 lead 간 **EAS504 r²=.007372, JPT104 r²≈.02375**; 별도 신호 가능성을 시사하지만 formal conditional 독립성은 아직 없음.

**최신 EAS/JPT lead tagging:** chr10 EAS504에서 r²≥.8 SNP 수 BBJ lead **6**, GIGASTROKE lead **11**; JPT104에서는 각각 **8, 11**. CALHM2 EBI nominal QTL p<1e-4 SNP 가운데 r²≥.8인 변이는 BBJ lead 태그 **3개**, GIGASTROKE 태그 **0개** (두 집단 모두). 이 집합은 **LD tag**이지 SuSiE credible set이 아님. 단일 신호 exploratory ABF 가정 N175/p12=1e-5에서 CALHM2 PP.H4≈.905, NEURL≈.061; 공식 multi-signal coloc 아님.

EUR503에서 나타난 heuristic 잔여 후보 `10:105694487:G:T`는 EAS/JPT에서 같은 cut-off를 통과하지 못함: alpha=.05 approximate z **EUR5.498, EAS4.011, JPT3.508**. 참조 LD 기반 approximate residual z는 정식 조건부 GWAS p-value가 아님. 2,000회 reference-genotype bootstrap은 참조패널 불확실성만 나타냄.

### COJO/인과 추론 블로커
공개 EAS504·JPT104는 GCTA가 제시한 **>4,000명 수준의 참조표본 권장 규모에 크게 미달**하며 1000G LD는 COJO용으로 권장되지 않음. 서버에서 GCTA 실행 바이너리도 당시 PATH에 없었음. BBJ raw per-SNP N=174,686을 소스에서 보존했지만 binary-trait case/control N 및 overlap 확인 필요. BBJ 7,472명 WGS LD panel(NBDC `JGAS000738/JGAD000873`)은 **controlled access**로 접근 권한 없음. EBI/GTEx molecular donor LD와 QTL effective N도 미검증. 따라서 정식 conditional GWAS, SuSiE fine-mapping, multi-signal coloc, 확정적 MR/인과 표현 **BLOCKED**.

## 4. 기완료 작업·마지막 검증
- BBJ/GIGASTROKE full regional GWAS 전수 스캔, 표준 GRCh37 및 reciprocal allele/p lookup 완료.
- EAS504/JPT104/EUR503 ALT dosage LD, palindrome 위치 및 bootstrap 민감도 완료.
- chr4/chr10 BBJ·GIGASTROKE·QTD000161 multi-track SVG 생성 및 GitHub 커밋 `286bd0f`.
- chr10 두 lead LD+CALHM2 QTL 감사 커밋 **`9c1a641`**, 코드 `scripts/is/audit_is_chr10_two_lead_ld_tags.py`, 문서 `docs/IS_CHR10_TWO_LEAD_EAS_JPT_LD_TAGS_20261011.md`, 테스트 `tests/test_is_chr10_two_lead_ld_tags.py`.
- 기준 시점 repo root에서 `python3 -m unittest discover -s tests -p 'test_is_*.py' -q` 실행 **329개 통과**. 변경 후 반드시 재검증.
- 역사적 **646개 coloc**, 새로운 **120개 ABF 민감도 분석**, **2,225 후보 유전자** 원본은 보존.
- **작업 소유권 주의:** `scripts/is/audit_is_col4a_eas_eur_ais_cross_ancestry_transfer.py` 및 `server/is_audit_col4a_eur_ais_v1.sh`는 타 작업의 untracked 파일이다. 임의로 수정/추가/삭제하지 않는다.

## 5. 남은 작업 — Harness 우선순위

**P0 — 공개자료만으로 즉시 진행:** CALHM2 nominal QTL p<1e-4에 해당하는 **정확한 3개 SNP**의 GRCh37/38 좌표·REF/ALT·eQTL beta/SE/p·BBJ/GIGASTROKE beta/p·EAS/JPT lead별 r²을 개별 표로 생성. 3개가 동일 LD cluster인지 확인해 pseudo-replication 배제. NEURL 대조, FGF5 양성 비교 포함. 신규 일자별 감사 폴더+manifest+Figure+테스트+보고서 Git commit.

**P1 — 통계적 확증의 자료 접근성:** 공식 BBJ/Japanese 규모 큰 ancestry-matched LD 및 독립 conditional summary GWAS 공개 여부 확인. EBI/GTEx molecular QTL donor LD, effective N 및 tested cis universe 검증. controlled 자료 필요시 `BLOCKED_REFERENCE_DATA` 명시하고 진행 가능한 다른 공개 분석 수행. 승인 없이 restricted BBJ genotype 접근 시도 금지.

**P2 — Broad discovery:** 4개 유전자에 한정하지 말고 2,225 유전자/전체 lead locus 후보를 출판된 IS 연구 수준으로 평가. genomic position, molecular pQTL/eQTL, colocalization 증거 등급, ancestry, 세포·조직 annotation 및 선행 논문 비교를 locus별 정리.

**P3 — 논문/웹:** Figure, data QC, GWAS-QTL/LD 신호, 가정적 PP.H4 민감도, prior/N 조건, negative results, limitations, source literature를 방법·결과·논의에 반영. 연구 웹은 배포 수단과 도메인을 확인한 후 업데이트 및 배포 검증.

## 6. 실행 안전 계약
1. 기존 canonicals 및 결과 디렉터리는 **READ ONLY**; 역사적 646 coloc·120 ABF·2,225 discovery 후보 **절대 덮어쓰지 않기**.
2. 임의의 p12/effective N·EAS/EUR/JPT LD를 정식 GWAS·QTL cohort LD로 설명하지 않는다.
3. 각 SNP는 좌표 빌드, REF/ALT, effect allele, palindrome, signed beta, p/SE, EAF 및 샘플 overlap 검증. GWAS nominal signal / LD tag / credible set / causal colocalization을 분리해 표시.
4. 모든 신규 결과는 고유한 날짜·버전 suffix의 결과 폴더, provenance JSON, TSV, 테스트, 문서에 저장.
5. `cd /srv/is-analysis/worktrees/is-broad-discovery-20261009`에서 전체 `test_is_*.py` 테스트를 실행하고 성공 내역만 보고.
6. `git status --short` 및 remote divergence 확인; **본인이 생성한 파일만 git add**, fetch/rebase 후 push, force push 금지. Unrelated untracked 파일 보존.
7. GitHub나 연구 웹 배포는 실제 커밋 SHA/배포 상태 확인 전 성공 주장 금지.

## 7. 최종 완료 조건 (Definition of Done)

**A. PUBLIC_DATA_AUDIT_COMPLETE**: 2,225-gene 범위를 유지한 locus ledger, chr4/chr10 full regional GWAS와 분자 QTL 개별 변이 audit, CALHM2 3-SNP LD cluster/두 lead 민감도, 원본 출처·allele/strand QC, 버전고정 Figure+TSV+Methods+negative controls, 테스트/CI 전부 통과, GitHub 커밋 및 웹 사용 시 실제 배포 검증. 접근 불가능한 자료로 인한 공식 분석은 명시적으로 **BLOCKED**, 가짜 결론 없음.

**B. FORMAL_CAUSAL_VALIDATION_COMPLETE**: 충분히 크고 적절한 일본인/BBJ LD genotype, 완전한 GWAS summary 및 binary-trait N·overlap 검증, study-compatible GTEx QTL N/LD, 조건부 GWAS·SuSiE credible set 수렴 및 다중 signal coloc/prior sensitivity, 필요 시 MR pleiotropy/heterogeneity와 기능/조직 replication 및 지도교수 검토. 충족 전 **논문용 확정 인과 유전자**라고 보고하지 않는다.

**현재 결론:** `PUBLIC_DATA_EXPLORATORY_COMPLETE / CALHM2_3SNP_AUDIT_NEXT / FORMAL_CAUSAL_CLAIM_BLOCKED`.

## 8. 첫 Harness 실행 순서
`git status`/HEAD 확인 → input/manifest 존재·hash 점검 → chr10 CALHM2 3개 QTL SNP 개별 조회 → EAS/JPT 두 lead r² 교차 계산 → 새 TSV+Figure+보고서 생성 → `python3 -m unittest discover -s tests -p 'test_is_*.py' -q` → 본인 파일만 커밋·푸시 → 완료사항/블로커/다음 단계 보고. **확인되지 않은 결과는 완료로 표시하지 않는다.**
