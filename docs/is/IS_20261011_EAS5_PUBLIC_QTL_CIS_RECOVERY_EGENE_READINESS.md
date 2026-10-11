# IS 연구 — 공개 QTL cis 검정 범위 완전성 복구 및 eGene 근거 분리 (2026-10-11)

**분석 브랜치:** `research/is-broad-discovery-20261009`
**연구 판정:** QTL 입력 검정 범위 및 공식 permutation 근거 **VERIFIED**(아래 한계 내), IS 공동위치화·인과 유전자 **NOT ESTABLISHED**.
**보존한 기존 결과:** 2,225개 고유 위치 기반 후보, 동아시아 영역 869개 유전자, 기존 BBJ–GTEx coloc 646개 분석(43개 유전자), BBJ와 GIGASTROKE GWAS, 15개 첫 공개 QTL 추출본. 모든 원자료 수정 금지.

## 연구 배경: 왜 새로운 coloc을 즉시 실행하지 않았는가?

앞선 5개 IS 잠정 영역에서 eQTL Catalogue 공용 nominal QTL 15개 파일을 확보하고 BBJ 일본인 IS 및 GIGASTROKE EAS AIS의 GRCh37 대립유전자와 매칭했다. 기존 ±500kb 주변 QTL 추출본에는 **원래 cis 검정 대상의 일부만** 포함되어 있었다. 따라서 GWAS와 일치한 SNP 수가 많다는 이유로 바로 coloc.abf 또는 coloc.susie를 실행하면 gene-to-SNP 검정 범위를 체계적으로 잘라내는 문제가 발생한다.

독립된 eQTL Catalogue `*.permuted.tsv.gz`의 유전자별 `n_variants`를 읽어 각 QTL 데이터셋의 원래 cis-tested universe와 대조했다. 분류의 단위는 **240개 위치 기반 유전자 × 3개 분자 QTL 데이터셋 = 720개 조합**이다. 해당 연구에서는 `permuted`의 공식 `variant` lead ID가 확장 nominal 입력 파일에 포함되는지도 검사했다. 다만 동일 변이 *수* 및 대표 SNP 포함은 공식 전체 개별 변이 목록과의 완전한 동일성 입증보다 약한 검증이다.

## 1. 공용 소스 및 기술적 입력 확보

| 공용 molecular-QTL 연구 | accession | EBI permutation 연구 전체 유전자 기록 수 | 감사 범위 |
|---|---|---:|---|
| GTEx v8 Artery Aorta 재처리 | QTD000131 | 21,888 | EAS 5개 영역 |
| GTEx v8 Brain Cortex 재처리 | QTD000171 | 24,675 | EAS 5개 영역 |
| ROSMAP DLPFC | QTD000434 | 23,076 | EAS 5개 영역 |

공식 public endpoint는 `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{QTS}/{QTD}/{QTD}.permuted.tsv.gz`이다. GTEx 2개는 `QTS000015`, ROSMAP은 `QTS000025`. 다운로드된 3개 파일은 원본 SHA-256을 기록했고, genomic source provenance는 보존했다. 이는 **공개 재처리 분자-QTL 결과**이므로 기존 GTEx V8 FastQTL의 동일한 donor 분석이라는 증거는 아니다.

대상 영역의 새 nominal QTL raw 자료는 공식 EBI 인덱스 `tabix`로 GRCh38 확장 genomic window에서 확보했다. 각 대상 유전자에 대해 이전 15개 파일의 모든 record의 **중복까지 고려한 multiset 포함**을 검증한 후 SNP allele set으로 중복 제거하여 독립된 `n_variants`와 비교했다. 공식 permutation의 lead variant genomic ID가 포함되지 않거나 과도한 변이 수가 발견되면 실패 처리하도록 했다.

## 2. 720개 검정 조합의 전체 cis 입력 복구

| 잠정 IS EAS 영역 | 위치 기반 유전자 | 조합(유전자×3 QTL 연구) | 확장 전 변이 수 일치 | 확장 후 변이 수 일치 | 공식 permutation 유전자 기록 없음 |
|---|---:|---:|---:|---:|---:|
| IS_XDATA_G0015 | 71 | 213 | 10 | **121** | 92 |
| IS_XDATA_G0023 | 14 | 42 | 0 | **22** | 20 |
| IS_XDATA_G0007 | 10 | 30 | 0 | **11** | 19 |
| IS_XDATA_G0025 | 92 | 276 | 0 | **167** | 109 |
| IS_XDATA_G0018 | 53 | 159 | 0 | **72** | 87 |
| **총합** | **240** | **720** | **10** | **393** | **327** |

기존에는 **383개**가 부분 cis 추출본이었다. 확장 후에는 **393/393개의 공식 검정 대상**에서 `n_variants`가 일치하고 공식 permutation lead genomic variant도 새 nominal 파일에 존재한다. 결과적으로 **383개 조합의 기술적 추출 공백을 복구**했다.

G0018의 `AC129929.5`(ENSG00000238184)는 이전 3Mb 조회의 마지막 좌표 `GRCh38 chr11:3377434` 오른쪽에 남아있던 **QTD000131 5개, QTD000171 4개, QTD000434 4개**의 유전자-SNP 기록으로 인해 전체 cis 검정 수가 일치하지 않았다. 공용 EBI 구간 `11:3377435-3387434`를 추가 조회하여 기존 입력에 순서대로 붙인 **새 파일**만 생성했다. 원본 3MB file은 건드리지 않았으며 3개 dataset 모두 G0018의 최종 72/72를 통과했다.

**327개는 음성 eQTL 유전자로 간주하지 않는다.** 해당 공개 데이터셋의 permutation 유전자 검정 기록을 확인하지 못한 조합이므로 원인(발현량, 필터, phenotype/QC, assay scope)은 추가 조사 대상이다.

## 3. 전체 연구 유전자 대상의 FDR 검증 — IS 인과성과 분리

공식 permutation 파일의 `p_beta`를 각 원본 분자-QTL 연구 **전체 유전자 검정 집합**(2~2.5만여 개)에 대해 Benjamini–Hochberg 보정했다. IS 후보만 선택한 뒤 재보정하여 q값을 유리하게 만들지 않았다.

| 통계 단위 | 결과 |
|---|---:|
| 공식 cis 검정 기록 있는 EAS5 유전자×QTL dataset | **393** |
| 전체 연구 단위 BH q≤0.05인 cis-eQTL 유전자×dataset | **98** |
| 세 데이터셋 중 하나 이상에서 FDR 통과한 고유 위치 기반 유전자 | **66** |
| 위 66개 중 기존 646회 GTEx–BBJ ABF 분석에 포함되지 않았던 유전자 | **56** |

| 영역 | QTL source eGene(FDR q≤0.05) | 기존 coloc 미포함 eGene | GWAS 근거 |
|---|---:|---:|---|
| G0015 | 22 | **16** | T1, BBJ·EAS AIS 유전체 유의 신호가 같은 넓은 검색 구간에 존재 |
| G0023 | 4 | 0 | T1 |
| G0007 | 0 | 0 | T1; 세 QTL source에서는 genomewide eGene FDR 통과 없음 |
| G0025 | 27 | **27** | T4, AIS 연관 신호는 탐색적 |
| G0018 | 13 | **13** | T4, AIS 연관 신호는 탐색적 |
| **합계** | **66** | **56** | T1과 T4 분리 |

새로운 T1의 eGene 16개는 모두 G0015 영역의 유전자이다. 대표 유전자 `C10orf32`, `GSTO2`, `AS3MT`, `TMEM180`, `GSTO1`, `CNNM2` 등은 **재처리된 분자 QTL 연구에서 유전자 발현 연관을 확인한 것**이지 IS GWAS와 동일 인과변이를 공유한다거나 임상 뇌졸중 기전을 의미하지 않는다. T4의 40개 신규 분자-QTL eGene에 대해서도 같은 주의가 필요하다. 어느 유전자도 새 인과 후보로 승격하지 않았다.

원래 393 gene–tissue 좁은 구간의 GWAS matched nonpalindromic SNP에 대해 `QTL P≤1e−4` 및 `EAS AIS P<5e−8`를 **같은 물리적 SNP에서** 요구했을 때 충족 조합은 **0/393**이었다. 이 조건은 정식 coloc 검정이 아니며, 이 0을 공동 인과변이 부재의 증거로 주장할 수도 없다. 확장 QTL 자료의 모든 SNP에 대해 GWAS REF/ALT 재매칭을 다시 수행한 것은 **아직 아님**.

## 4. AIS 표본 구성 감사

로컬 EBI YAML `GCST90104545`의 `sample_size=256,274`와 후속 출판물 보고 수치의 합을 대조했다: **19,032 AIS cases**, **237,242 controls**, 환자 비율 **0.0742642640**. 대조로 AS는 27,413 cases +237,242 controls =264,655. 참고: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11354504/> (AIS), <https://kp4cd.org/node/1138> (AS), <https://www.nature.com/articles/s41586-022-05165-3> (GIGASTROKE 출판본).

**분명한 한계:** AIS 환자비율은 로컬 합계와 후속 논문 표본 구성의 교차 확인이지, 원 논문의 개별 cohort별 참여자와 변이별 유효 표본 수/겹침을 재검증한 결과가 아니다. GIGASTROKE EAS는 BBJ 참여자 일부가 포함될 수 있으므로 독립적인 BBJ 재현이라고 표현하지 않는다.

## 5. 연구 상태 및 다음 coloc 실행 전 차단 기준

**VERIFIED:** 15개 원래 QTL source input 및 추가 15개 extended-cis genomic slice 무결성, EBI 공식 cis-tested `n_variants`, 공식 permutation 대표 genomic SNP의 nominal 입력 포함, 기존 raw 행 multiset 보존, 393 cis 검정 조합의 분자 QTL 기술적 입력 범위, 각 연구별 원 유전자 전체를 사용한 BH FDR.

**EXPLORATORY:** 66개 source-eGene과 새 56개 eGene, T1의 신규 source 16개. 뇌졸중 질환-인과성과 분리해야 하는 자료 가용성 및 분자 eQTL 증거.

**BLOCKED/NOT VERIFIED:** 확장 후 추가된 QTL SNP의 GRCh38→37 FASTA 및 GWAS 대립유전자 매칭 재실행, 원 QTL donor genotype 구성과 ancestry·LD, GIGASTROKE original case-control source-by-cohort table 및 BBJ 표본 비독립성, SNP-wise complete colocalization posterior/SMR/multisignal SuSiE, 허혈성 뇌졸중 환자조직 발현 효과, 최종 인과 유전자 지정.

**기존 논문 후보 집합 2,225개와 646개 coloc posterior는 변경하지 않았다.** 원본 QTL 압축 파일, 리포트 전용 derivative, GitHub 소스 및 자동화 거래·CSS 서비스의 경계도 그대로 유지했다.

## 6. 재현성·코드·원자료 위치

공식 permutation 3개 원본:
`/srv/is-analysis/data/is/qtl/eqtl_catalogue/eas5_permuted_20261011_v1`

추가 raw QTL:
`/srv/is-analysis/data/is/qtl/eqtl_catalogue/eas5_G0015_extended_cis_20261011_v1`,
`.../eas5_T1_extended_cis_20261011_v1`,
`.../eas5_T4_extended_cis_20261011_v1`,
`.../eas5_G0018_edge_cis_20261011_v1`.

실행 감사 테이블:
- `/srv/is-analysis/results/is/audits/is_eas5_full_cis_tested_universe_20261011_v1` (원본 720 조합 cis completeness)
- `/srv/is-analysis/results/is/audits/is_eas5_extended_cis_consistent_20261011_v3` (5개 영역 × 3 소스 per-gene completeness)
- `/srv/is-analysis/results/is/audits/is_eas5_full_cis_master_readiness_20261011_v1` (전후 720 통합 상태와 Figure)
- `/srv/is-analysis/results/is/audits/is_eas5_studywide_BH_egenes_20261011_v1` (240개 후보 QTL eGene BH)
- `/srv/is-analysis/results/is/audits/is_eas5_240_egenes_gwas_triage_20261011_v1` (GWAS T1/T4 분리 QTL 근거표)
- `/srv/is-analysis/results/is/audits/is_gigastroke_eas_ais_case_fraction_20261011_v1` (표본비 검증)
- [논문용 cis completeness before–after Figure](../figures/is/IS_EAS5_720_FULL_CIS_RECOVERY_FIGURE.svg)

코드:
- `scripts/is/audit_is_eas5_qtl_full_cis_permuted_universe.py`
- `scripts/is/extend_is_eas5_public_qtl_complete_cis.py`
- `scripts/is/merge_is_eas5_full_cis_readiness.py`
- `scripts/is/audit_is_eas5_permuted_studywide_bh_egenes.py`
- `scripts/is/build_is_eas5_240_qtl_gwas_source_triage.py`
- `scripts/is/audit_is_gigastroke_eas_ais_case_fraction.py`
- `scripts/is/render_is_eas5_full_cis_recovery.R`

공개하고 있는 유전자·조직 TSV는 **원본의 일치 검증 결과만** 보관한 작은 분석표이며, 대용량 raw EBI QTL 및 private/controlled genotype 자료는 GitHub에 업로드하지 않는다.
