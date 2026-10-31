# OCULAR MASTER — 연구 방향 재설정 (2026-10-31)

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
