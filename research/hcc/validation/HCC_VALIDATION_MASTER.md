# HCC Multi-omics — 연구 검증·논문 유의성 MASTER v1.0

**작성:** 2026-10-30 KST · **용도:** 분석 재개 시 가장 먼저 확인할 연구 운영 문서  
**상태:** PROTOCOL, 아직 논문의 연구가설·통계결과를 확증하지 않음  
**관리 문서:** research/hcc/validation/checkpoints.json (기계 판독 기준), research/hcc/validation/DECISION_LOG.md (의사결정 이력)  
**기존 실행 작업:** IS_Analysis_V3 PR #251, Google Drive HCC_MULTIOMICS_FEASIBILITY  
**주의:** 기존 HCC Radiogenomics 수업의 제출 양식은 AI 작성 제한이 명시된 바 있으므로 이 문서는 **내부 기술 검증 기준서**이지 제출 원고가 아니다.

## 1. 최종 연구 질문과 주장 범위

**우선 가설 RQ1:** HCC의 간세포 대사·분화·콜레스테롤/지질·증식 신호에 대해, RNA와 Protein의 경로 수준 일치/불일치가 병인(ALD, HCV, HBV, NAFLD 등)에 따라 달라지는가? 그 차이는 **독립 코호트에서 재현되는가**?

**후속 가설 RQ2:** Bulk RNA–Protein 경로 차이가 종양 세포의 상태, 면역·기질 세포 구성, 조직 위치와 어떤 관계가 있는가? 단일세포에서는 환자 단위 pseudobulk 또는 환자 단위 혼합모형으로 검증한다.

**조건부 가설 RQ3:** 외부 재발 정보가 적절한 경우, 사전 확정된 생물학적 pathway 점수가 병기·연령·AFP·치료 등을 포함한 **기존 임상모델 대비** 독립 검증에서 재발 예측을 개선하는가?

**데이터 접근 실패 시 RQ4:** 환자별 RNA–Protein 매칭이 불가능하면, 독립 연구 간 **pathway-level triangulation/replication**으로 전환한다. 이 경우 “matched proteotranscriptomic integration”이라는 문구 사용 금지.

연구의 핵심 기여는 (A) 병인 간 신호의 **이식성/비이식성**, (B) 세포 기원의 독립적인 정량 근거, (C) 충분한 자료가 존재할 때만 예후 모델의 임상적 부가가치이다. 단순 새 유전자 리스트·새 아형·UMAP 이미지만으로 주요 신규 발견을 주장하지 않는다.

## 2. 선행 논문 → 재현할 결과와 반드시 피할 중복

| ID | 핵심 논문/근거 | 논문 결과/방법에서 확인할 체크포인트 | 기존 결과를 넘어서는 우리 검증 |
|---|---|---|---|
| L01 | Jiang et al. Nature 2019, DOI 10.1038/s41586-019-0987-8 | HBV 관련 조기 HCC 110 paired tumor/non-tumor, 3개 단백체 아형, 불량 예후 S-III, 콜레스테롤 항상성 및 SOAT1, 기능적 실험 | SOAT1·콜레스테롤 pathway는 **양성 대조군**, 독립 코호트·병인으로의 전이 평가 |
| L02 | Gao et al. Cell 2019, DOI 10.1016/j.cell.2019.08.052 | HBV-HCC 159명, 종양/인접 조직, proteomic 3아형, 대사 변화, PYCR2/ADH1A 예후 관련, CTNNB1–ALDOA phosphorylation | Gao 코호트는 독립 HBV 맥락의 재현성 검사; 논문에 나온 유전자 자체를 신규 발견으로 소개 금지 |
| L03 | Ng et al. Nat Commun 2022, DOI 10.1038/s41467-022-29960-8 | 114명 122 종양 RNA/WES, 그중 51 종양 global proteome/phosphoproteome; cell cycle, metabolic repression, Wnt/AKT/mTOR/Notch, CTNNB1/TP53 연관 인산화; 일부 전사체–단백체 불일치 | 병인별 effect modification을 **동일 환자 단위**로 추정하고 신호의 외부 이식성 및 세포기원 해석. 원 논문이 이미 한 분석을 novelty로 주장 금지 |
| L04 | Lu et al. GEO GSE149614 | HCC scRNA 환자 10명·샘플 21개·71,915 cell metadata; 종양/인접/PVTT/전이 site 구분 | 환자 단위 pseudobulk + 조직부위 통제. 세포 하나하나를 독립 환자로 취급하지 않음 |
| L05 | Liu et al. Cell 2018 TCGA-CDR, DOI 10.1016/j.cell.2018.02.052 | OS/DSS/DFI/PFI 표준화 및 암종별 endpoint 권고 | RFS/TTR/OS를 혼용하지 않음. 시간영점·사건 정의·치료 적격성 검증 |
| L06 | TRIPOD+AI BMJ 2024, DOI 10.1136/bmj-2023-078378; PROBAST+AI BMJ 2025, DOI 10.1136/bmj-2024-082505 | 예측모델의 모집단, predictor/outcome, 누출, 외부 검증, calibration·discrimination·bias 평가 | 내부 CV는 외부 검증이 아님. 임상모델 대비 추가 가치가 없으면 임상 유용성 주장 금지 |
| L07 | J Proteome Res 2025 DOI 10.1021/acs.jproteome.5c00828 | HBV-HCC proteomic/transcriptomic subtype과 예후 연구 추가 존재 | 투고 전 2025~현재 novelty search 재수행, 단순 HBV multi-omics subtype 재발견 회피 |

**고임팩트 연구와의 범위 정합성:** 논문의 핵심적인 paired design, QC, multi-omic concordance, pathway, external evidence, clinical confounding을 커버한다. 새로운 기전의 실험적 인과성은 사용자 연구에서 별도 실험 없이는 주장하지 않는다. DOI 및 GSE 링크는 논문 인용·근거용이지 데이터 이용권 보증이 아니다.

## 3. 검증 게이트 — 실행 증거와 실패 시 연구 전환

| Gate | 검증 시점 | 반드시 확인할 수치/증거 | PASS 판정 근거 | 실패 시 처리 |
|---|---|---|---|---|
| V00 | 시작 | 원 논문 DOI, data accession, 이용권, 공개/제한 자료 구분 | 라이선스·접근 경로·출처 문서화 | 공개 CI에서 EGA/민감 자료 취급 금지 |
| V01 | 수집 | 데이터 파일 URL·SHA256, format, assay/platform, genome/gene ID | 실제 파일 스키마 + 다운로드 로그 | 재시도 또는 다른 공개 코호트 선택 |
| V02 | **가장 중요** | RNA 발현 매트릭스와 단백질 정량 매트릭스에 **동일 조직 biopsy ID가 존재**, 각각 distinct patient/biopsy N | 데이터값이 있는 두 매트릭스의 교집합 산출물 | 매칭 불가 시 RQ4로 전환 |
| V03 | 분석 준비 | RNA/Protein QC, batch, 결측, feature mapping, tumor/adjacent 상태, 중복 | QC 보고서 + 샘플 제외 사유 및 시각화 | 교란/품질 문제 시 비교 중단 |
| V04 | 병인 비교 | ALD/HCV/HBV/NAFLD/MASH 단독·복합 병인, 환자별 n, missingness | 환자 단위 codebook, 생물학적 복수 biopsy 클러스터 관리 | n 부족하면 기술·탐색 연구로 축소 |
| V05 | 교란 검증 | BCLC, grading, cirrhosis, tumor purity, sex/age, treatment, batch/site | 균형 표 + 조정/감도 분석 | 조정 불가능하면 ‘병인 효과’ 주장 금지 |
| V06 | 양성 대조군 | Jiang SOAT1/cholesterol; Gao PYCR2/ADH1A, ALDOA; Ng metabolism / cell-cycle | 측정 가능한 feature 기준으로 방향·효과·CI를 제시, 차이 설명 | assay mapping/QC 원인 분석 |
| V07 | RNA–Protein | 같은 환자/조직의 1:1 gene–protein 매핑, 상관/불일치/표본 수/CI | 사전 정의 pathway, 부트스트랩, FDR 명시 | 매핑 약하면 경로 수준 결과로 제한 |
| V08 | 병인×경로 | 직접적인 group interaction coefficient, CI, q-value | 집단 비교 차이의 통계적 검정 + sensitivity | 한 집단만 p<0.05여도 상호작용 주장 금지 |
| V09 | 외부 재현 | 다른 환자·연구·assay에서 방향·효과크기·이질성 확인 | 재현 로그 + heterogeneity | ‘study-specific’로 표기 |
| V10 | 세포기원 | scRNA patient/sample/barcode QC, site 통제, pseudobulk | 최소 환자별 효과/CI, 환자 영향도 확인 | UMAP/Dotplot은 기술적 참고자료로만 |
| V11 | 공간분석 선택 | 치료상태·조직 부위와 spot QC, 독립환자 | 적절한 comparable spatial tissue | 불가능하면 과감히 생략 |
| V12 | 임상 endpoint | recurrence event/status codebook, time zero, follow-up, censoring, patient n/events | patient-level event table + 적격성 설명 | RQ3 제거, OS/DFI/PFI로 임의 치환 금지 |
| V13 | 임상모델 | 고정 score, baseline 임상모델, 진정한 외부 코호트 | ΔC-index+CI, calibration, time-AUC/Brier/decision utility | HR 또는 KM p만으로 ‘임상적 유용성’ 불가 |
| V14 | 강건성 | alternate normalization, study/batch/site, negative controls, label permutation | 1차 효과 방향/크기 안정성 | 탐색으로 격하/한계 인정 |
| V15 | 재현성·보고 | source hashes, code SHA, config, seed, environment, figures+tables, null results | CI + bio evidence ledger + 누출 검증 | 논문 확정 보류 |

**게이트 상태:** NOT_STARTED / PARTIAL / PASS / BLOCKED / FAIL. PASS는 원본으로부터 직접 계산한 **필수 증거**를 연결한 경우에만 설정한다. 논문 발표 수치, 저장소 파일명, CI 녹색 표시만으로 승격 금지. 데이터 접근 불능은 BLOCKED, 충분한 정밀도 아래 예측 실패는 FAIL, 작은 n에서 귀무가설 기각 실패는 대개 INCONCLUSIVE이며 FAIL과 다르다.

## 4. 통계적 판단 규칙 — 분석 실행 전에 잠금

**사전 가설:** 간세포 분화·대사(HNF4A/ALB/APOA), 지질·콜레스테롤 대사, 증식/cell-cycle, 면역/myeloid, 혈관/기질을 우선 후보로 둔다. 문헌에 제시된 경로는 **재현 대조군**과 신규 **병인 상호작용 탐색**을 구분한다. 경로 정의/가족 수, 결측 기준, 유전자 매핑은 결과를 보기 전에 freeze.

- 고차원 탐색: BH-FDR q<0.05를 잠정 기본값으로 하고 effect size·95% CI·가용 환자 수를 함께 보고. 작은 사전 확정 확인 가설은 다중검정을 반영한 별도 가족으로 관리. 임계값은 통계적 유효성 보증이 아니며, 실제 n/power 검토 후 고정해야 한다.
- **병인별 유의성 차이 ≠ 유의한 병인 상호작용**. 실제 interaction term과 CI 평가.
- 다중 종양을 가진 환자 및 단일세포의 세포 반복은 patient cluster와 pseudobulk로 통제. 무작위 split을 biopsy/cell 레벨로 하면 leakage.
- RNA–Protein discordance는 MS missingness, 단백질 검출성, batch/tissue purity, feature identity를 먼저 검사.
- 외부 코호트에서 **방향이 일치**하면서 크기·CI·이질성이 충분한지 확인. 외부 검증 데이터를 feature selection/threshold 튜닝에 사용하지 않는다.
- 재발 예후: Cox PH 진단, competing event 검토, censoring/immortal time, 사건 수와 모델 복잡도, patient-level nested CV와 별도 외부 샘플. 작은 event N에서는 확증적인 예측모델을 추진하지 않음.

## 5. 검증 결과별 논문의 유의성·주장 가능 범위

| 경로 | 실제 검증 결과 | 논문에서 허용할 주장 | 연구 가치 / Next action |
|---|---|---|---|
| **A** | 동일 환자 멀티오믹스 + 병인 interaction + 독립 proteome replication + 환자단위 cell origin | 병인별 차등 이식성, 세포 환경과의 연관 | **가장 높은 우선순위**의 관찰연구. 실험적 인과성은 주장 불가 |
| **B** | paired + 여러 코호트에서 pathway 보존, interaction null + **좁은 CI** | 병인 간 보편적인 대사 신호와 효과 한계 | 검정력이 충분할 때 가치 있는 negative/boundary finding |
| **C** | 내부 q<0.05, 외부 효과 역방향 또는 0 근처 | 내부 코호트 특이적 연관 | 일반화된 바이오마커 아님. batch/stage/assay 조사 |
| **D** | 비유의하지만 CI 넓고 n 작음 | 판단 불가 / underpowered | ‘차이 없음’ 주장 금지, 비교 축소 |
| **E** | RNA/Protein 정량값 조인 불가 | 독립 코호트 pathway triangulation | paired integration 문구 삭제, 연구 주제 변경 |
| **F** | 세포별 패턴이 한 환자/부위에 의존 | 기술적 cell context | 세포기원 검증 통과로 표시 금지 |
| **G** | 재발 정의 또는 독립 사건 코호트 부재 | biological reproducibility | Clinical prediction 부록/가설로만 |
| **H** | Cox HR 유의하지만 ΔC-index 미확인/외부 calibration 불량 | outcome 연관 존재 가능 | 임상 추가 가치 미입증 |
| **I** | 통계는 우수하나 Jiang/Gao/Ng 또는 최신 문헌과 동일 발견 | independent replication | novelty 낮음; 새로운 병인·코호트·기전 추가 없으면 상위 저널 주장 금지 |
| **J** | 이식성/세포 맥락/독립 재발예측 incremental validation 전부 충분 | 생물학적 재현성과 외부 임상 부가가치 | 높은 수준의 translational 관찰연구 후보, **전향적 유용성은 별도** |

논문의 "유의성"은 통계 q 값과 연구 독창성, 생물학적 설득력, 다른 환자군에 대한 재현성, 실제 임상적 부가가치를 **서로 독립적인 축**으로 평가한다. 저널 Impact Factor/게재는 보장할 수 없다.

## 6. 현재 실사 사실(확인과 미확인 구분)

| 자료 | 직접 확인한 항목 | 아직 확인되지 않은 것 |
|---|---|---|
| Ng 2022 | 114 patient / 122 tumor biopsies; 같은 biopsy에서 RNA/Proteome/Phosphoproteome 측정 Y인 51 tumor = 49 patients | 두 실제 **정량 matrix**의 gene×sample join. EGA restricted RNA / public processed matrix 여부 |
| Ng 병인(49명) | ALD only 13, HCV only 11, NAFLD only 7, HBV only 4, HCV+ALD 7, ALD+HBV 3, ALD+NAFLD 2, no liver disease 2 | BCLC/therapy/batch confounding 및 interaction power |
| Gao | PDC000198/199 각 171 repository cases, 331 aliquots. 2019 Cell 논문은 paired patients 159명 | 171/159 불일치 원인 및 PDC Protein Assembly ↔ BioSino OEP000321 RNA patient mapping |
| Gao 후보 파일 | Zhou_Liver_Cancer_Proteome.tmt11.tsv 및 .sample.txt; Phosphoproteome.phosphosite.tmt11.tsv | 실제 matrix values 및 동일 환자 ID |
| Ng PRIDE | PXD025705 198개 / PXD025836 195개 파일 목록 확인. Spectronaut SWATH/Progenesis Result XLSX 발견 | full quantitative matrix format, biopsy ID mapping |
| GSE14520 | 488 sample rows; recurrence status/time 242 rows, 원본 값 1=136/0=106 | 실제 환자 수, status code meaning, RFS/event eligibility |
| GSE149614 | 71,915 cells / 21 samples / 10 patients | patient-level pseudobulk 효과와 independent replication |

이 수치는 **2026-10-30 실사 증거**이며 연구 유의성 결과가 아니다. 상세: https://github.com/kimtk94/IS_Analysis_V3/pull/251 및 해당 PR의 GitHub Actions artifact. 안전 임계값(최소 paired 30명, 병인군 각 10명, recurrence event 30/20 등)은 초기 **engineering screening heuristic**이지 power calculation 대체 불가.

## 7. 다음 작업 인수인계

1. 최초 작업 시 이 문서 + checkpoints.json + 기존 feasibility/datasets.json을 읽고 기록된 마지막 evidence timestamp 확인.
2. 해당 단계에 속하는 첫 **PARTIAL/BLOCKED** 중 실행 가능한 public-data checkpoint 하나 이상 처리.
3. 코드·원본 버전·환자 수·분석 결과·실패 이유를 증거 파일로 생성하고 CI 검증.
4. **PASS 승격은 증거 산출물 존재와 내용 확인 후에만**, 검증 결과를 DECISION_LOG와 Google Drive에 반영.
5. 통계 결과가 유의/비유의인지뿐 아니라 **어느 주장이 가능해졌고 어느 주장은 폐기되는지** 반드시 명시.
6. 매칭 없이 RNA/Protein 통합, 환자 검증 없이 세포 기반 유의성, 외부 검증 없이 예측 성능 주장은 명시적으로 금지.

**NEXT P0:** PDC Protein Assembly 실제 정량 TSV 헤더와 샘플 매핑; Ng SWATH/Progenesis 워크북 정량 matrix 존재 여부; legal processed RNA mapping.  
**NEXT P1:** GSE14520 환자 단위 재발 정의 검증, GSE149614 pseudobulk 설계.  
**Go/No-go:** paired 정량 매트릭스 매칭 실패시 RQ4로 전환, 재발 cohort 실패시 RQ3 삭제.

## Sources
- Jiang 2019: https://doi.org/10.1038/s41586-019-0987-8
- Gao 2019: https://doi.org/10.1016/j.cell.2019.08.052
- Ng 2022: https://doi.org/10.1038/s41467-022-29960-8
- HCC single-cell: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614
- TCGA CDR: https://doi.org/10.1016/j.cell.2018.02.052
- TRIPOD+AI: https://doi.org/10.1136/bmj-2023-078378
- PROBAST+AI: https://doi.org/10.1136/bmj-2024-082505
- 2025 novelty sentinel: https://doi.org/10.1021/acs.jproteome.5c00828

## 8. Drive research handoff (2026-10-30)

**Google Drive living copy (internal reference):** https://docs.google.com/document/d/1Ipk7kS6ukGYqwgwhSWHICs42EBqpu3eRFd82mBOkaAE/edit

This Drive document is a readable summary; the canonical evidence status and fail-safe claim eligibility logic are in the GitHub checkpoint registry and CI. Both governance PR #252 and molecular feasibility PR #251 remain draft until reviewed and merged. If the copies differ, do not silently overwrite either; reconcile with a dated DECISION_LOG entry.
