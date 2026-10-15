---
title: "MASTER DEGREE — 원홍희 교수 선행연구 기반 IS·CKD 연구 방향 및 High-impact workflow 종료 기준"
date: "2026-10-15"
status: "RESEARCH_DIRECTION_V1 / AUDIT_REQUIRED"
projects: ["ischemic_stroke", "chronic_kidney_disease"]
owner: "MASTER_DEGREE"
scope_policy: "high-impact-paper-relevant-workflow only"
evidence_policy: "documented != scientifically validated"
source_repository: "kimtk94/IS_Analysis_V3"
---

# 원홍희 교수 선행연구 기반 IS·CKD 연구 방향 (2026-10-15)

> **핵심 결정:** 새로운 분석법을 계속 덧붙이지 않는다. *Nature*, *Nature Medicine*, *Nature Genetics*급 **질문 관련 workflow**와 원홍희(Hong-Hee Won) 교수의 **유전학→molecular QTL→MR/SMR→colocalization→조직·세포 근거→후보 우선순위화**라는 연구 문법을 충족하면 분석을 중단하고 논문 작성으로 전환한다. 저널 IF는 문헌 선별 힌트일 뿐 분석 타당성의 직접 척도가 아니다.

## 0. 문서 적용 범위와 출처/한계

- 본 문서는 **연구 방향과 중단 기준**을 확정하는 설계 문서이며, 서버 원자료와 전체 결과를 직접 재실행하여 검증한 문서가 아니다.
- 기존 원홍희 교수 자료: Google Drive [147편 verified working catalog](https://docs.google.com/spreadsheets/d/1eL6uYOjrDhOAKYhQvuuVciHW40r3GiqADnDGSVf9P6k/edit), [핵심 25편](https://docs.google.com/spreadsheets/d/1v-nKAIrbEwJ7JVI436uh2naVb0tXbx_n51d6b5aEIQ0/edit), [Top 8 상세 리뷰](https://docs.google.com/document/d/1KvOz9LzlhSOV-O69CCVGUNeVOyAnR9UdH5Q101PlDXs/edit), [교수 연구 종합 리뷰](https://docs.google.com/document/d/1n6YS4cgwW1HU20qyOT8nkpIYSbGyFa6nhSJGKhMoPCY/edit). 147편은 **공식 완전 출판목록이 아닌 기존 working catalog**다.
- 기존 상위 저널 기준 문서: [HIGH_IMPACT_WORKFLOW_BASELINE_20261015.md](./HIGH_IMPACT_WORKFLOW_BASELINE_20261015.md).
- 프로젝트 현황 참조: [IS_MASTER.md](./IS_MASTER.md), [ckd_stage1_mr.md](./ckd_stage1_mr.md), [ckd_stage2c_susie.md](./ckd_stage2c_susie.md), [ckd_stage3a_kidney.md](./ckd_stage3a_kidney.md).
- **연구주제 행정적 최종 확정은 별개**: 이전 논의에서는 CKD를 논문 후보로 결정했으나, 2026-10-05 GitHub IS_MASTER에는 IS가 thesis_main으로 표기됨. 어느 것이 최종 학위논문인지 지도교수와 결론을 내리기 전에는 두 연구의 진척을 하나의 검증된 학위논문으로 혼용하지 않는다.
- **날짜 주의:** GitHub main의 일부 문서는 2026-10-05 또는 그 이전 스냅샷이고, 이후 IS 분석은 별도 브랜치/서버에서 진행됐을 수 있다. 아래 상태는 **문서 기반 provisional**이며 서버 실행 산출물 감사로 업데이트해야 한다.

## 1. 연구 방향 결론

### IS (허혈성 뇌졸중) — 추천 질문

**영문 가제:** *Ancestry-aware integrative genetic and molecular prioritization of ischemic-stroke susceptibility loci with vascular cell-type context*

- 주요 질문: **동아시아(BBJ)의 ischemic-stroke loci에서**, 기존 서구 중심 GWAS로 충분히 설명되지 않는 후보 유전자/단백질 기전을 **대조 가능한 molecular QTL과 뇌혈관 세포 근거**로 우선순위화할 수 있는가?
- 흐름: BBJ/EAS stroke GWAS → 전 유의 locus 적격성 QC → EUR/EAS 방향/LD 비교 → credible-set/조건부 분석(가능 범위) → brain/artery eQTL·plasma pQTL → MR/SMR/coloc → brain/vascular sc/snRNA 기능적 해석 → 최종 증거표.
- **FGF5, ALDH2, SH3PXD2A, COL4A2**는 핵심 기전 예시일 뿐 초기 발견 후보군을 4개로 제한하지 않는다. CALHM2, NEURL1, C4orf22, INA, COL4A1 등은 탐색군으로 유지하고 사전 정의된 선정 규칙을 적용한다.
- **ALDH2 rs671**은 coding/protein/metabolic branch: bulk eQTL 부재를 실패로 처리하지 않는다. rs671 통계·allele·conditional/LD 근거는 별도 QC한다.
- 세포 분석의 목적은 **어느 vascular/immune/mural cell type에서 후보의 기능적 개연성이 있는지**를 설명하는 것. expression-only는 causal regulatory link와 구분한다.

### CKD (만성콩팥병) — 추천 질문

**영문 가제:** *Proteogenomic prioritization of chronic kidney disease susceptibility proteins integrating plasma pQTL, causal inference, kidney tissue evidence and East Asian replication*

- 주요 질문: 혈장 cis-pQTL 유래 단백질–신기능 연관성 중 **MR와 regional colocalization**에서 지지되고 **신장 tissue QTL/세포 위치**와 일관적인 후보는 무엇인가? 동아시아 데이터에서 얼마나 재현되는가?
- 흐름: UKB-PPP cis-pQTL → eGFRcrea/CKD/BUN 등 harmonization → MR/Wald·다중검정 → ABF coloc → 필요시 SuSiE/LD 민감도 → kidney tissue pQTL/eQTL → kidney cell atlas → EAS outcome 검증 → 후보 증거 등급.
- 기존 9개 seed: **ACP1, CPVL, F12, GSTA1, GSTA3, HLA-E, INHBC, SDCCAG8, UMOD**. Seed 선정 시점 및 전체 eligible protein screen을 보존한다.
- INHBC/UMOD/SDCCAG8/HLA-E는 2024 Si et al.에서 선행 근거가 존재하므로 신규 발견이라는 표현 금지; 추가 다층근거 및 집단별 재현성에 기여를 둔다.
- KoGES controlled-access 개인 genotype/phenotype이 없으면 **종단 eGFR 및 incident CKD 직접 재현을 주장하지 않는다**. 공개 데이터 기반 논문은 별도 완성 가능하며 KoGES는 데이터 확보 시 추가 검증으로 표기한다.

## 2. 원홍희 교수 논문: 직접 적용할 핵심 6편

| 우선 | 원홍희 교수 참여 논문 | 서지/DOI | 원 논문의 실제 workflow | IS/CKD에서 가져올 것 | 그대로 따라 하지 않을 것 |
|---|---|---|---|---|---|
| P0 | **Kim et al. (2023)**, *Prioritization of therapeutic targets for dyslipidemia using integrative multi-omics and multi-trait analysis* | *Cell Reports Medicine*, DOI [10.1016/j.xcrm.2023.101112](https://doi.org/10.1016/j.xcrm.2023.101112), PMID 37582372 | Lipid PRS→PheWAS→MR/MVMR, 병렬로 eQTL SMR/HEIDI→pQTL/network→target scoring→phenome safety; 30 targets | 다양한 근거를 **gene-wise evidence matrix**로 보존하고 우선순위화하는 논문 서술/figure 구조 | 모든 프로젝트에 1,129 Phecodes PheWAS, PRS, MVMR, network score 7종을 일괄 의무화하지 않음 |
| P0 | **Lee et al. (2025)**, *Identification of potential therapeutic targets for problematic alcohol use using multi-omics data* | *eBioMedicine*, DOI [10.1016/j.ebiom.2025.105939](https://doi.org/10.1016/j.ebiom.2025.105939), PMID 40966854 | Brain/blood eQTL·pQTL→TWMR/PWMR→HEIDI/pleiotropy check→Bayesian coloc→resource concordance→druggability/PheWMR | **MR 단독 hit이 아닌 sensitivity+coloc+독립 근거로 타깃 선정**; ancestry 별 LD를 분리하는 원칙 | 후보 근거 없는 전체 phenome-wide MR, 모든 hit에 druggability scoring |
| P0 | **Kim et al. (2026)**, *Genetic architecture of plasma pTau217 and related biomarkers in Alzheimer's disease via GWAS* | *Alzheimer's & Dementia*, DOI [10.1002/alz.71275](https://doi.org/10.1002/alz.71275), PMID 41804841 | 한국인 K-ROAD biomarker GWAS (n=1,972)→brain QTL→SMR/coloc→Korean snRNA cell context | 동아시아 cohort에서 **locus→regulatory gene→cell type** 근거 사슬 구성, 한 후보 figure의 모범 | 뇌졸중 또는 CKD에 무관한 AD-specific biomarker 실험 |
| P0 | **Jang et al. (2026)**, *A meta-analysis of single-nucleus expression QTL linking genetic risk to brain disorders* | *Nature Genetics*, DOI [10.1038/s41588-026-02541-x](https://doi.org/10.1038/s41588-026-02541-x), 2026-04-07 publisher correction 존재 | 4 cohorts, 유럽인 983명·핵 580만개→cell type/subtype eQTL→GWAS coloc/MR/fine-mapping; SingleBrain resource 공개 | IS의 **microglia/neuron/glia cell-eQTL annotation** 및 bulk eQTL null 해석; 원시 재분석 대신 **공개 QTL summary** 사용 | 580만 핵 eQTL meta-analysis 자체 재수행, 새 모델 개발; 혈관세포 표현형에 SingleBrain을 무리하게 대입하지 않음 |
| P1 | **Cho et al. (2025)**, *Functional insight into East Asian-specific genetic risk loci for Alzheimer's disease* | *Alzheimer's & Dementia*, DOI [10.1002/alz.14553](https://doi.org/10.1002/alz.14553), PMID 39991798 | EAS AD 7 GWAS/26 loci→cell type eQTL·peak-to-gene links→functional gene inference | IS의 EAS locus annotation 시 **variant→regulatory peak→gene→cell** 연결이 가능하면 우선 보고 | 모든 locus의 peak-to-gene link를 데이터 없이 강제 생산하지 않음 |
| P1 | **Won 참여 (2024)**, *Multivariate genomic analysis of 5 million people elucidates shared MetS components* | *Nature Genetics*, DOI [10.1038/s41588-024-01933-1](https://doi.org/10.1038/s41588-024-01933-1) | MetS multiple traits→Genomic SEM multivariate GWAS→functional annotation→MR/PheWAS/PRS & cross-ancestry | 복수 phenotype 정합성과 **EAS replication의 별도 해석**; 대사 회복탄력성 연구 설계 참고 | IS·CKD에 genomic SEM/PRS 추가는 필수 아님. 현재 분석 질문을 바꿀 때만 선택 |

> 원홍희 교수 논문은 '필수분석 목록'이 아니라 **연구실의 증거 적층·통계 보고·Figure 설계 양식**을 배우는 기준이다. 저널 IF가 높은 다른 대표 연구는 GWAS/조직 분야의 **분석 검증 기준**으로 함께 사용한다.

### 교수 논문에서 직접 확인한 방법론/통계 기준 (논문마다 다른 기준을 통합해 오용 금지)

| 방법 | 원 논문 사례 | 우리 적용 규칙 |
|---|---|---|
| cis-QTL instrument | PAU 2025: eQTL cis ±1 Mb, GWAS 수준 P<5×10⁻⁸ | cis 창, allele, P, strength/F statistics를 분석별 사전 정의. 기존 CKD cis pQTL 정의와 혼용 금지 |
| MR | Dyslipidemia 2023: IVW, weighted median, MR-Egger/MR-PRESSO 민감도 분석; PAU 2025: instrument 조건에 맞춰 TWMR/PWMR 수행 | k=1 Wald ratio, 적절한 다중 IV는 IVW/가중중앙값, MR-Egger는 충분한 IV에서만 해석. 검정력 없는 sensitivity test를 억지로 PASS 처리하지 않음 |
| SMR·HEIDI | Dyslipidemia 2023: 40,302 gene/tissue tests에 대해 SMR Bonferroni ≈1.24×10⁻⁶를 보고. PAU 2025: HEIDI P<0.01 heterogeneity 신호 | **HEIDI의 귀무가설과 결정방향 반드시 확인**. HEIDI 비유의성=공통 인과변이 증명 아님. 원 논문의 test-family 및 다중검정 분모 그대로 복사 금지 |
| Bayesian coloc | PAU 2025: PP.H4>0.75를 strong colocalization 기준으로 사용 | 기존 프로젝트의 주 플래그 PP.H4≥0.8과 비교; 0.75~0.8은 참고 구간으로 별도 표기. PP.H3, prior sensitivity, SNP coverage 및 credible set을 병행 |
| SuSiE | Nature/현 파이프라인: multiple-signal loci에서 coloc.susie 사용 | 좋은 LD reference + 양 trait의 credible set 수렴 확인 시 사용. UMOD 같은 미수렴 loci는 **UNRESOLVED**로 남기고 무한 재시도 금지 |
| ancestry validation | PAU 2025: 비유럽인 GWAS에 ancestry-matched LD 경로; pTau217 2026: 한국인 GWAS와 brain QTL·Korean snRNA 결합 | BBJ EAS/JPT와 EUR LD 구분, allele frequencies·effect 방향 비교. EAS outcome support와 EAS QTL colocalization은 **서로 다른 검증** |
| multiple tests | 각 논문에서 Bonferroni/FDR 분모가 연구 설계와 다름 | discovery family/targeted validation family를 분리해 분모 보고. nominal P만으로 최종 target 선언 금지 |
| cell context | SingleBrain 2026: eQTL/coloc; East Asian AD 2025: eQTL + peak-to-gene | (1) expression/localization (2) QTL (3) allelic regulation을 다른 층으로 명시. **표현=인과작용 아님** |

## 3. 상위 저널 benchmark와 교수 논문 교차 매핑

| 문제/단계 | High-impact benchmark | 원홍희 교수 reference | 최소 충족 산출물 |
|---|---|---|---|
| IS EAS/다인종 GWAS | Mishra et al., *Nature* 2022 ([DOI](https://doi.org/10.1038/s41586-022-05165-3)) | Cho 2025 EAS locus functionalization, Kim 2026 Korean biomarker GWAS | 적격 locus **전체 목록**, EAS/EUR comparable-variant 표, allele harmonization 및 locus plots |
| IS molecular target | Mishra 2022 functional gene prioritization | Lee 2025 PAU eQTL/pQTL–MR–coloc | MR/SMR/ABF/SuSiE 단계별 evidence/failure 표, gene별 주장 등급 |
| IS vascular cell context | Wälchli et al., *Nature* 2024 ([DOI](https://doi.org/10.1038/s41586-024-07493-y)) | Jang 2026 SingleBrain, Cho 2025 peak-to-gene | 적절한 혈관/면역세포 atlas 기반 expression/QTL/ATAC 링크 및 한계 |
| CKD tissue proteogenomics | Hirohama et al., *Nature Medicine* 2025 ([DOI](https://doi.org/10.1038/s41591-025-03872-8)) | Kim 2023 multilevel drug-target priority, Lee 2025 MR/coloc | plasma MR/coloc, kidney tissue pQTL/eQTL, 9 genes evidence integration |
| CKD kidney cell resolution | Sheng et al., *Nature Genetics* 2021 ([DOI](https://doi.org/10.1038/s41588-021-00909-9)) | Jang 2026와 Cho 2025의 세포 수준 검증 문법 | kidney cell-type localization과 가능한 cell-type QTL; expression-only 별도 표시 |
| Cross-ancestry | Mishra 2022 / CKD 공개 EAS GWAS | 2024 MetS cross-ancestry, 2025 PAU ancestry-matched LD | ancestry-specific evidence & LD provenance 명시 |

### 피해야 할 과잉 개발

1. 초기 후보를 무한 확장하거나 신규 분석 엔진·모델을 계속 추가하지 않는다. **eligible universe QC**가 끝나면 그 결과 안에서 통계적으로 우선순위화한다.
2. PRS, Genomic SEM, 대규모 PheWAS, 약물재창출·off-target 안전성 스캔은 연구의 **주요 주장**에 필수일 때만 추가한다. 교수님 논문에서 수행했다고 자동으로 의무화하지 않는다.
3. SingleBrain/brain vascular atlas/신장 atlas는 **공개 데이터 재사용**이 기본. 수백만 세포 재분석이나 신규 wet-lab 실험은 필수 단계가 아니다.
4. SuSiE가 수렴하지 않는다고 LD 불일치 상태에서 반복횟수를 무한정 늘리지 않는다.
5. 개인 genotype 부재 상태에서 KoGES 직접 종단 복제를 완료했다고 쓰지 않는다.

## 4. IS 단계별 실행·중단 기준

| Gate | 질문 | 완료 판정에 필요한 최소 증거 | 현 상태(문서 기준) |
|---|---|---|---|
| IS-01 GWAS scope | 모든 적격 BBJ/EAS 신호에서 출발했나? | 고정된 discovery 입력 버전, 포함/제외 규칙, 전체 후보 table, allele/build QC, checksum | PARTIAL: 4 core gene 이외 broad screen 검증 필요 |
| IS-02 LD/fine-map | ancestry mismatch를 적절히 처리했나? | reference LD 출처, 비교 가능 SNP, credible set 또는 unresolved/sensitivity 기록, EAS/EUR effect 비교 | PARTIAL: BBJ cohort-specific LD 부재 명시 |
| IS-03 QTL/causal | molecular gene evidence가 일치하나? | eQTL/pQTL 출처, MR/SMR, coloc PP.H3/PP.H4/prior, 다중신호 QC, 실패 reason | PARTIAL: 문서상 구현 일부 확인; 실행물 확인 필요 |
| IS-04 mechanism | gene과 관련 세포 맥락을 설명하나? | FGF5·ALDH2·SH3PXD2A·COL4A2 + 독립 선정 후보의 세포 localization / QTL / ATAC evidence level | PARTIAL: Phase10 목표, 미확인 결과를 COMPLETE로 표시 금지 |
| IS-05 manuscript | 재현 가능한 논문 그림인가? | Figure 1 flow; Figure 2 ancestry/locus; Figure 3 molecular evidence; Figure 4 cell-type/prioritization. Supplementary QC+full candidate table | NOT AUDITED |

**IS 정지 조건:** IS-01~05에 대해 산출물이 모두 PASS, 혹은 외부 LD 데이터 불가 때문에 UNRESOLVED인 부분은 claim을 한정하여 명확히 설명하고 핵심 결론을 유지할 수 있을 때. **cortex sn-eQTL ≠ cerebrovascular endothelial eQTL**이며 데이터가 맞지 않으면 annotation only.

## 5. CKD 단계별 실행·중단 기준

| Gate | 질문 | 완료 판정에 필요한 최소 증거 | 현 상태(문서 기준) |
|---|---|---|---|
| CKD-01 exposure/outcome | pQTL·GWAS 입력 정의와 allele이 유효한가? | UKB-PPP ST16 provenance, cis instrument/QC, F statistics, eGFR/CKD/BUN phenotype harmonization | DOCUMENTED, runtime PASS 미확인 |
| CKD-02 causal | MR·coloc 해석이 재현 가능한가? | Wald/IVW/WM/MR-Egger 조건별 선택, multiple-testing, coloc ABF, SuSiE 수렴/credible sets | PARTIAL: 특정 locus(UMOD 등) unresolved |
| CKD-03 kidney tissue | 혈장 단백질이 신장 조직에서 추가 근거가 있나? | Hirohama 2025 supplementary + kidney eQTL 일치·불일치/coverage table | DOCUMENTED, 후보별 출력물 최종감사 필요 |
| CKD-04 cell and ancestry | cell specificity와 EAS signal에 일관성이 있나? | HPA/kidney atlas 정확한 cell labels, EAS outcome direction, EUR/EAS LD provenance | PARTIAL |
| CKD-05 replication claim | longitudinal Korean confirmation을 주장하나? | KoGES controlled phenotype/genotype와 종단 모델 결과 **또는 NOT AVAILABLE 별도 명시** | BLOCKED by access, 논문 전체는 별도 추진 가능 |
| CKD-06 manuscript | figure + candidate novelty audit 완료? | Figure 1 overview; Figure 2 MR screen; Figure 3 coloc examples; Figure 4 tissue/cell & EAS evidence; 별도 Si 2024 overlap table | NOT AUDITED |

**CKD 정지 조건:** CKD-01~04,06은 통계·출처·그림까지 검토 가능한 상태; CKD-05가 불가할 때 'KoGES validation' 주장을 삭제하고 **public-data-only scope**로 종료. 이미 보고된 UMOD/SDCCAG8/INHBC/HLA-E를 신규 유전자라고 하지 않음.

## 6. 최종 Evidence tier / 해석 규칙

- **Tier A — strongest supported:** 정확한 allele/LD/QC + 다중검정 통과한 molecular association + 사전 지정된 MR/SMR 방향성과 통계 + 지역 coloc/credible-set/민감도 검토 + 적절한 tissue/cell 또는 독립집단 근거. 다만 **causal gene proven**이라고 단정하지 않음.
- **Tier B — convergent but incomplete:** GWAS/QTL+MR 또는 coloc이 충분하지만 LD mismatch, tissue specificity, 독립 validation 등의 제한이 있음.
- **Tier C — exploratory:** expression-only, 단일 명목 P, 약한 coloc, alternative candidate, unresolved credible set. 구조적으로 제거하지 않고 Supplementary에 남김.
- **Unresolved/blocked:** 미수렴, 자료 부족, SNP 매칭 실패는 결과 실패가 아닌 **측정·식별 한계**로 관리. 모든 후보를 Tier A에 올리는 것이 목표가 아님.
- 근거가 모두 유의하더라도 종간 조직·phenotype 불일치와 horizontal pleiotropy를 무시하지 않는다. PP.H4가 높다는 이유로 MR의 exclusion restriction 또는 기능 기전이 자동 입증되는 것은 아니다.

## 7. 지도교수 검토를 위한 연구방향 요약 (질문 포함)

**연구 핵심 제안:** "공개 다인종 GWAS와 혈장/조직 QTL을 활용하여 IS 또는 CKD의 유전적으로 지지되는 후보 유전자·단백질을 우선순위화하고, MR·colocalization 및 적절한 세포 자료로 근거를 보강한다. 검정력·LD·데이터 접근 제약을 논문에 명시하며 상위저널의 질문 관련 분석 workflow를 충족하는 선에서 추가 개발을 종료한다."

**검토 요청 포인트:**
1. 학위논문 **1순위 질환은 IS인가 CKD인가**, 두 연구를 분리할지?
2. 원홍희 교수 연구실 방식 중 **원고에 반드시 포함해야 할 SMR/HEIDI, pQTL MR, regional coloc, validation 조합**이 무엇인가?
3. **동아시아 replication만으로 충분한가**, 아니면 controlled-access 한국인 longitudinal phenotype 분석을 학위 요건으로 필수 지정할 것인가?
4. ISSUES: EAS cohort-specific LD가 없는 IS locus, CKD SuSiE 미수렴 locus의 **해석·보고 수준**이 적절한가?
5. 유전적 지지 타깃 발굴 연구에서 **PheWAS/druggability를 본문에 넣을지**, Supplementary/후속연구로 미룰지? 기본은 후자.

## 8. 즉시 진행할 검증 작업 — 추가 분석이 아니라 갭 클로징

1. GitHub 최신 **research/is-broad-discovery**, **research/is-evidence-audit**, **feat/is-harness**, **research/ckd-thesis** 브랜치와 서버 artifact를 기준으로 IS/CKD 표의 PASS/PARTIAL/BLOCKED를 업데이트한다. 메인 브랜치의 오래된 docs만으로 결론 내리지 않는다.
2. IS 전 locus candidate universe·allele audit, EAS/EUR/LD provenance를 **한 표**로 통합한다.
3. CKD 9개 gene마다 MR, ABF, SuSiE convergence, tissue QTL, HPA, EAS replication, prior-literature overlap을 **한 표**로 통합한다.
4. 선택된 **학위논문 1개**에 대해 4개 main figures를 제작하고 원문 및 Supplementary 결과를 검증한다. 실제 자료 없는 민감도 검사나 wet-lab 단계는 신설하지 않는다.
5. Figures와 Methods의 용어 사용(association/support/shared signal/causal)을 일치시키고 analysis freeze를 선언한다.

## 9. 주요 문헌 및 링크 (검증 가능한 원문·DOI)

### 원홍희 교수 연구

1. Kim MS et al. (2023). *Cell Reports Medicine*. **10.1016/j.xcrm.2023.101112**. [원문](https://pmc.ncbi.nlm.nih.gov/articles/PMC10518515/) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/37582372/). Figure 1: two-phase discovery + prioritization; Figure 3: molecular target scoring.
2. Lee DJ et al. (2025). *eBioMedicine*. **10.1016/j.ebiom.2025.105939**. [원문](https://pmc.ncbi.nlm.nih.gov/articles/PMC12478272/) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/40966854/). MR+HEIDI/pleiotropy+coloc PP.H4>0.75; cross-resource replication.
3. Kim JP et al. (2026). *Alzheimer's & Dementia*. **10.1002/alz.71275**. [PubMed](https://pubmed.ncbi.nlm.nih.gov/41804841/). Korean plasma-biomarker GWAS→brain molecular evidence.
4. Jang B et al. (2026). *Nature Genetics*. **10.1038/s41588-026-02541-x**. [원문](https://www.nature.com/articles/s41588-026-02541-x) | [SingleBrain eQTL data](https://doi.org/10.5281/zenodo.14908182). Correction: **10.1038/s41588-026-02581-3**.
5. Cho M et al. (2025). *Alzheimer's & Dementia*. **10.1002/alz.14553**. [원문](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848548/).
6. (2024) *Nature Genetics* multivariate MetS. **10.1038/s41588-024-01933-1**. [원문](https://www.nature.com/articles/s41588-024-01933-1).

### 질환별 high-impact 비교 문헌

7. Mishra A et al. (2022). *Nature*. **10.1038/s41586-022-05165-3**. [원문](https://pmc.ncbi.nlm.nih.gov/articles/PMC9524349/).
8. Wälchli T et al. (2024). *Nature*. **10.1038/s41586-024-07493-y**. [원문](https://www.nature.com/articles/s41586-024-07493-y).
9. Hirohama D et al. (2025). *Nature Medicine*. **10.1038/s41591-025-03872-8**. [원문](https://pmc.ncbi.nlm.nih.gov/articles/PMC12435990/). Kidney pQTL/eQTL, 337 samples; supplementary 1–30.
10. Sheng X et al. (2021). *Nature Genetics*. **10.1038/s41588-021-00909-9**. [원문](https://pmc.ncbi.nlm.nih.gov/articles/PMC9338440/).
11. Si J et al. (2024). *Genome Medicine*. **10.1186/s13073-024-01356-x**. [원문](https://link.springer.com/article/10.1186/s13073-024-01356-x). CKD prior-target overlap.

---

**Document acceptance criteria:** 연구별 질문, 분석에 기여하는 원홍희 교수 논문, 고IF peer benchmark, 통계 해석 방향, 현재 제약, 중단/완료 조건 및 Figure가 모두 한 문서에 연결되어 있어야 한다. **코드 존재 또는 실행 성공 ≠ 과학적 검증**.
