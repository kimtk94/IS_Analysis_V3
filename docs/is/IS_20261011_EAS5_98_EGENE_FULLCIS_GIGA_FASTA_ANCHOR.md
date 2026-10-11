# IS — 66개 공개 eGene의 완전한 cis 범위와 GIGASTROKE FASTA 기반 대립유전자 감사

**실행일:** 2026-10-11 (KST)
**연구 범위:** 5개 잠정 IS 동아시아 영역 중 3개 공개 QTL 연구에서 연구별 BH-FDR q≤0.05를 통과한 66개 eGene, 98개 유전자×조직 비교
**판정:** 입력 provenance 및 SNP 교집합은 **VERIFIED**, 동시 유의성은 **EXPLORATORY**, 신규 공동위치화와 IS 인과성은 **NOT ESTABLISHED**.

## 1. 왜 재분석했는가?

기존 15개 공개 nominal QTL 파일은 각 영역 주변 약 2Mb만 잘라낸 자료였다. 공식 eQTL Catalogue 유전자별 `permuted.tsv.gz`의 `n_variants`와 대표 SNP ID를 대조한 결과, 전체 720개 조합 가운데 permutation 검정 기록을 가진 393개 조합 전부에서 **전체 cis SNP 개수와 permutation 대표 SNP가 확인되도록 확장**했다. 미보고 327개는 음성 eQTL로 취급하지 않았다. 출처: `docs/is/IS_20261011_EAS5_PUBLIC_QTL_CIS_RECOVERY_EGENE_READINESS.md`.

그 후 공식 permutation `p_beta`를 **각 원래 QTL 연구의 모든 유전자**에 대해 Benjamini–Hochberg 보정하여 98 gene–QTL tissue assays, 66개의 eGene을 확보했다. 기존 646개의 coloc posterior와 최초 2,225개 IS 후보 유전자 목록은 변경하지 않았다.

### source eGene 수와 선행 연구 비교

GTEx v8 Aorta/Cortex 및 ROSMAP DLPFC의 eGene 검정은 **분자 형질과 cis 유전형의 연관**이지, 뇌졸중 환자의 차등발현이나 GWAS 공유 인과변이를 뜻하지 않는다. 본 GTEx v8 연구는 838명 중 **715명(85.3%)이 European American**이라고 보고했으며, 본 387명 Aorta·205명 Cortex 재처리 자료의 세부 유전조상 조합과 LD는 독립 검증하지 않았다. 참고: GTEx Consortium, *Science* 2020, DOI 10.1126/science.aaz1776; <https://pmc.ncbi.nlm.nih.gov/articles/PMC7737656/>.

eQTL Catalogue 공식 정의상 GRCh38 `ALT`는 beta의 **effect allele**이며 `AC/AN`은 대체대립유전자 빈도이다. 또한 dbSNP 중복 rsID에 따른 동일 genomic variant의 여러 행은 **coloc 이전에 중복 제거**해야 한다. 참고: <https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md> 및 <https://www.ebi.ac.uk/eqtl/Data_access/>.

## 2. 확장된 cis 영역의 입력 QC

- 각 유전자–조직의 전체 원 QTL 변이 개수 대조, 공식 permutation 대표 SNP의 nominal 파일 존재, 기존 원본 파일의 전 행 multiset 포함 검증.
- 98개 FDR+ 조합에서 **490,373개**의 유전자–QTL 데이터셋–SNP 매핑 기록이 고유 매핑 및 REF, ALT, beta/SE/P, AC/AN 수치 검증을 통과한 뒤 후속 GWAS 검색 대상이 됐다. 이는 독립 490,373개 SNP라는 뜻이 아니다.
- GRCh38→GRCh37 liftover는 `hg38ToHg19` 역방향 체인 및 독립 GRCh37 `human_g1k_v37.fasta`를 사용했다.
- G0018에서 서로 다른 GRCh38 SNP가 같은 GRCh37 좌표 및 REF/ALT로 매핑되는 **11개 gene–dataset–SNP 충돌 키**가 확인됐다. 오류를 임의 합치거나 beta를 선택하지 않고, 두 원본 SNP를 모두 교집합에서 격리했다.
- BBJ GWAS 원본 canonical의 REF/ALT 및 효과대립유전자 관례는 따로 검증했다. GIGASTROKE EAS AIS는 다음 항목에서 별도 처리했다.

## 3. GIGASTROKE AIS의 역사적 variant ID와 실제 REF/ALT

원본 `GCST90104545_buildGRCh37.tsv.gz`는 **pre-GWAS-SSF**이고 genomic REF/ALT 열 대신 `effect_allele`, `other_allele`를 제공한다. 분석용 normalized 파일의 역사적 `variant_id`는 주로 `chr:pos:other:effect` 형태였다. 이를 단순히 `REF:ALT`로 해석하면 첫 번째 allele가 FASTA REF가 아닌 경우 상당수의 유효 변이를 잘못 제거한다.

개별 SNP에서 **FASTA REF가 반드시 원본 효과대립유전자 또는 기타 대립유전자 중 하나와 일치할 때만** 그 REF를 채택했고, 다른 하나를 ALT로 지정했다. Effect allele가 REF라면 beta에 -1을 곱하고 EAF를 1-EAF로 변환했다. 원본 역사적 ID가 `other:effect` 또는 이미 `REF:ALT`라는 구조로 인정되는지 독립 검증하여 불일치·불완전·다염기 SNP는 제외했다.

**원본 출처 검증**
- EBI raw source MD5: `904e12e33ec3c2eeefe4dc93bec4995c` **일치**.
- 메타데이터: GCST90104545, GRCh37, AIS, East Asian, N=256,274, non-harmonized.
- **원본 GWAS 6,725,216행과 정규화 GWAS 6,725,216행을 각각 전체 스캔**했다.
- 선택한 QTL SNP 위치의 **23,893개 GIGASTROKE 원본–정규화 SNP 쌍 전부**에서 EA/OA, beta, SE, P값, EAF의 수치 일치, 역사적 ID 규칙 및 독립 FASTA REF 기준이 통과했다.
- FASTA REF:ALT가 이미 일치하는 역사적 variant ID: **13,358개**.
- 역사적 other:effect allele 순서를 FASTA REF:ALT로 재정렬한 ID: **10,535개**.
- 원본·정규화 GWAS 파일은 수정하지 않았다. 재정렬된 결과는 별도 derivative에만 기록했다.

## 4. GWAS–QTL 교집합: 원래, 좁은 ID 판독, 출처 검증 후

| 잠정 IS 영역 | 원래 2Mb nominal QTL 교집합 | 전체 cis + 엄격한 역사적 ID 해석 | 전체 cis + FASTA REF 정규화 |
|---|---:|---:|---:|
| G0015 (T1) | 32,576 | 43,150 | **79,297** |
| G0023 (T1) | 6,522 | 11,881 | **20,762** |
| G0025 (T4) | 26,990 | 54,343 | **95,150** |
| G0018 (T4) | 20,877 | 35,893 | **60,544** |
| **합계** | **86,965** | **145,267** | **255,753** |

숫자의 단위는 **gene×molecular-QTL study×non-palindromic GRCh37 SNP** 기록이다. 여러 유전자·조직·LD 상관 SNP이 같은 위치에서 반복되므로 독립 유전체 변이 수와 다르다. T1 G0007은 세 QTL 자료에서 연구별 BH q≤0.05인 유전자–조직 조합이 없어 이 표에서는 0건이다.

이전 145,267건은 오염되거나 거짓인 SNP가 아닌데도 GIGASTROKE의 역사적 ID를 잘못 REF:ALT로 가정하여 **보수적으로 과소평가했던 결과**다. 공식 원본 및 FASTA 검증 후 **255,753건을 이 분석 범위의 canonical technical match count로 채택한다.**

- [세 상태의 SNP 교집합 Figure](../figures/is/IS_EAS5_GIGA_HISTORICAL_VARIANT_RECODING_THREE_WAY_COMPARISON.svg)
- [전체 cis 복구 전후 Figure](../figures/is/IS_EAS5_98_EGENE_FULLCIS_BBJ_AIS_ALLELE_EXACT_SNP_OVERLAP.svg)

## 5. 새 탐색 결과: COL4A2–ROSMAP DLPFC

실제 같은 *비-palindromic GRCh37 SNP*에서 동시에 `QTL P≤1e−4`와 `GIGASTROKE EAS AIS P<5e−8`인 기록은 **98개 유전자–조직 assay 중 1개, 2개 SNP**로 확인되었다. 이 조합은 **COL4A2**, `ENSG00000134871`, **ROSMAP Brain DLPFC** (`QTD000434`)이다.

| 정확한 GRCh37 REF:ALT SNP | QTL beta_ALT | QTL P | GIGA EAS AIS beta_ALT / P | BBJ beta_ALT / P | QTL vs EAS AIS ALT EAF 차이 |
|---|---:|---:|---:|---:|---:|
| `13:110954237:A:G` | -0.249388 | 6.19611e−5 | +0.0694 / 4.64e−8 | +0.06856 / 9.67322e−9 | 0.04023 |
| `13:110954510:T:C` | -0.253888 | 2.18498e−5 | +0.0693 / 1.007e−8 | +0.06839 / 2.91645e−9 | 0.03819 |

QTL `p_beta=0.00627568`, 전체 ROSMAP 분자 QTL 연구 유전자 BH `q=0.0171972`. COL4A2의 공식 cis 검정 변이 **8,764개** 가운데 GRCh37 FASTA 검증·bbj·EAS AIS 두 연구에서 공통 매칭된 non-palindromic SNP은 **3,947개**이다. 이 두 SNP의 GIGASTROKE 역사적 variant ID는 `OTHER_EFFECT_ORDER`로, 이번 FASTA 기반 재정렬을 거치면서 정확하게 검출됐다.

이 정보는 **좌위 수준에서 뇌졸중 GWAS와 COL4A2 eQTL가 같은 태그 SNP에서 동시 연관**됨을 기술적으로 지지한다. 다만 같은 LD 블록의 두 SNP일 수 있어 두 독립 신호가 아니며, 두 형질이 **단일 인과변이를 공유한다는 증거도 아직 아니다**. 역방향 연관 부호가 나타나지만 이것만으로 COL4A2의 발현 감소가 AIS를 유발한다고 해석하지 않는다. 유전자 발현이 환자 뇌조직에서 변화한다는 주장도 금지한다.

### 선행 분석과 직접 비교

COL4A2는 기존 **646 GTEx–BBJ ABF 분석에 이미 포함된 유전자(16 tissue assays)**이므로 *신규 발굴 유전자*가 아니다. 이전 GTEx–BBJ `p12` 사전확률 민감도 테이블(`legacy_coloc_p12_sensitivity_20261010_v1`)에서 16개 COL4A2 assay의 **최대 H4**는:

| coloc p12 | 16 assays 중 최대 PP.H4 |
|---|---:|
| 1e−6 | 0.0471 |
| 1e−5 | **0.3308** |
| 1e−4 | 0.8318 |

높은 사전확률에서는 H4가 크게 증가했지만, 이 이전 값들은 **다른 GTEx tissue 및 이전 BBJ 좌위**의 posterior 재가중 결과이며, 현재 ROSMAP DLPFC의 독립적인 coloc 계산값이 아니다. 높은 p12에서의 H4만 선택해 결론 내리면 안 된다.

## 6. 남아 있는 과학적 차단 기준

| 검증 게이트 | 이번 실행 결과 | 논문 해석 |
|---|---|---|
| eQTL Catalogue 원본 SNP·effect allele 규칙 | **VERIFIED** | ALT=효과대립유전자 |
| QTL 전체 cis 변이수 + 공식 permutation 대표 ID | **VERIFIED: 393/393** | 원 출판 SNP 목록의 완전 동일성 입증과는 다름 |
| 98 eGene assays 전체 cis SNP → GRCh37 FASTA | **VERIFIED for included records** | 11개의 비일대일 변환 키 격리 |
| GIGA 원본–canonical 23,893 SNP 전체 수치·FASTA 검증 | **VERIFIED** | 10,535개 역사적 ID 순서 재정렬 |
| BBJ/GIGASTROKE 공통 nonpalindromic SNP | **255,753 association records** | 중복 gene/tissue SNP 기록이며 독립 변이가 아님 |
| COL4A2 같은 SNP의 QTL/AIS 동시 유의 | **EXPLORATORY: 2 SNP / 1 assay** | 공동위치화 또는 인과성 미확정 |
| 원본 GIGA EAS AIS case/control ratio | **SECONDARY SOURCE CROSS-CHECK ONLY** | 19,032 / 237,242 / total 256,274 |
| 원본 QTL 유전형 cohort-matched ancestry/LD | **NOT VERIFIED** | SuSiE conditional/coloc 진단 제한 |
| BBJ·GIGASTROKE AIS 개인 단위 코호트 중복 | **NOT RESOLVED** | 독립 복제로 부르면 안 됨 |
| 신규 eGene coloc.abf/coloc.susie 및 MR | **NOT RUN** | 신규 H4 유전자 0개 |
| AIS 환자조직 대조 DGE | **NOT RUN** | 환자 발현 기전 주장 불가 |

### 이어서 진행할 우선 분석

1. **COL4A2의 3,947개 비-palindromic SNP**에서 유전자 발현 QTL와 BBJ/EAS AIS의 LD 일치, 잔여 변이 결측, multi-signal 구조를 검토한다.
2. 가능하면 GIGASTROKE EAS AIS 개별 cohort 분해와 원자료 환자/대조군 비율을 확인하고, EUR 중심 GTEx/ROSMAP QTL vs EAS AIS의 LD 교차 조상 차이를 계량한다.
3. 같은 locus에 있는 COL4A1 등 인접 유전자를 포함해 공유 변이/공유 LD 신호에 의한 유전자 후보 분리를 검토한다.
4. 선행 COL4A2의 16개 GTEx–BBJ ABF 분석과 **신규 ROSMAP DLPFC (EAS AIS outcome) coloc**의 phenotype, genotype, LD, SNP 커버리지 차이를 분리해서 제시한다. 향후 H4는 `p12` 민감도와 H3 대비 우위가 일관될 때만 탐색적 후보로 사용한다.
5. **현재 2,225개 광범위 후보집합이나 기존 646개 coloc posterior는 절대 자동 갱신하지 않는다.**

## 7. 재현성: GitHub 소스와 서버 원본

분석 코드:
- `scripts/is/audit_is_eas5_98_egene_fullcis_gwas_exact.py`
- `scripts/is/audit_is_eas5_giga_raw_canonical_fastaref.py`
- `scripts/is/gigastroke_fasta_variant_contract.py` (기존 공용 규칙 재사용)
- `server/is_run_eas5_98_egene_fullcis_gwas_overlap.sh` (원시 historical ID 해석을 엄격히 적용하는 비교군)
- `server/is_run_eas5_98_egene_fullcis_gwas_fastaanchored.sh` (FASTA 원본 검증 후 수정된 기술적 정규화)
- `scripts/is/render_is_eas5_98_egene_fullcis_overlap.R`
- `scripts/is/render_is_eas5_giga_fastaref_recoding_coverage.R`
- `tests/test_is_eas5_98_egene_fullcis_gwas_exact.py`
- `tests/test_is_eas5_giga_raw_canonical_fastaref.py`

서버 파생 결과:
- `/srv/is-analysis/results/is/audits/is_eas5_98_egene_fullcis_gwas_overlap_20261011_v2/` (strict historical ID baseline)
- `/srv/is-analysis/results/is/audits/is_eas5_98_egene_fullcis_gwas_overlap_20261011_v3_fastaanchored/` (FASTA 검증 갱신 결과; **연구상 우선 채택할 source QC 산출물**)
- `/srv/is-analysis/results/is/audits/is_eas5_98_egene_giga_original_source_20261011_v1/` (원본/derived 전수 대조)

그림·작은 유전자별 TSV는 GitHub에도 포함되지만, 대용량 20MB 압축 유전자–SNP 전체 중간 산출물은 서버에만 보관한다. 이 문서의 `VERIFIED` 판정은 *입력 데이터의 기술적 검증*에 한정된다. 논문 내 과학적 인과 결론은 **BLOCKED** 상태로 남는다.
