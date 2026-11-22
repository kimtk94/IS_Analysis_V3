# IS P3A chr10 — 실제 36개 QTL assay 전체 GWAS 원본·cis 커버리지 일괄 검증 (2026-11-22)

**범위:** MASTER DEGREE Ischemic Stroke. **분석 단위 구별:** P3A 후속 검토 `36 gene×QTL dataset groups` = `20 unique gene IDs` = **하나의 유전체 영역 `IS_XDATA_G0015`**. 절대로 36개 독립 GWAS locus 또는 36번의 독립 복제라고 표현하지 않는다.

## 1. 실제 검증한 원본 자료와 분석 대상

- 식별자: 불변 커밋 `0bf646fd6f29322d17da878253f37def24c01c74`에 보존된 393행 원본 QTL readiness에서 36개 추출; **20개 유전자**, 분포 QTD000131(10), QTD000171(11), QTD000434(15). 원본 393/98/393행 Git blob 해시 검증과 manifest 대조: [PR #303](https://github.com/kimtk94/IS_Analysis_V3/pull/303), [CI #45019987189](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45019987189).
- 유전체 위치는 기존 GWAS 연구의 [30-regions-six-outcomes 원본](https://github.com/kimtk94/IS_Analysis_V3/blob/0bf646fd6f29322d17da878253f37def24c01c74/docs/data/is/IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv)의 `G0015` 그룹으로 고정: **GRCh37 chr10:104,663,676–105,720,351**, 기술적 GWAS context 스캔만 **chr10:104,163,676–106,220,351**로 제한. **2Mb context는 모든 gene-specific full cis windows의 완전한 대체물이 아니다.**
- GIGASTROKE publisher GRCh37 full original ancestry GWAS: **EAS `GCST90104545`**, **EUR `GCST90104540`**. [앞서 whole-file publisher MD5 + full SHA256 + bytes 실제 검증 #45018300939](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45018300939). 이번 신규 CI는 원본 압축 파일을 임시 파일로 받아 다시 해시하고, **원본 `beta`, `standard_error`, `p_value`, `effect_allele`, `other_allele`, `effect_allele_frequency`**를 직접 읽어 G0015 영역 SNP-level **숫자형 필드 존재와 기술적 분포만** 확인한다. 개별 SNP 효과크기 행 자체는 GitHub에 게시하지 않는다.
- 해당 원본 GWAS에는 SE가 있지만, 과거 **P3A 서버의 QTL+GWAS 합본 24-column 파일에는 BBJ/EAS-AIS의 관측 SE/N/case fraction 필드가 없었다.** 두 파일의 통계값이 동일 코호트·표현형·allele 정의인지 **SNP별 원본 일치 검증 없이** GWAS SE를 옛 서버 P3A 입력에 조합하지 않는다.

## 2. 이미 완료된 원본 25개 public BH-FDR eGene 그룹의 full-cis 커버리지

별도 불변 98-assay official eGene BH05 원본과 교집합으로 잡힌 **25/36개**(원래 서버의 *source-stat-positive 25개*와 동일한 그룹이라고 검증되지 않음)의 기술적 합계:

| 출처에서 정해진 단위 | gene×QTL×SNP 기록 행 수 |
|---|---:|
| Public studywide BH q≤0.05로 교집합인 25개 그룹의 완전한 nominal cis 검사 수 | **142,196** |
| BBJ와 EAS에서 둘 다 매칭됐다는 기존 소스 보고 비회문형 SNP 행 | **59,090** |
| 그중 기존 QTL–EAS EAF 차이 ≤ 0.1 행 | **31,313** |

비회문형 교집합은 nominal cis 검사 대비 **41.56%**. 이는 **25개 유전자–조직별 검사 기록의 합계**, 중복 genomic variant를 제거한 값이나 실제 원본으로 새로 독립 FASTA/allele을 맞춰본 결과가 아니다. QTL source에 존재하지 않는 11개 그룹을 원천적으로 연구 불가 또는 음성 결과로 분류하면 안 된다. 이 자료의 25개 BH05 그룹 모두 이전 EAS 같은 SNP에서 `QTL p≤1e−4, AIS p<5e−8` 결합 기준 기록 **0건**이나, 그 필터만으로 유전적 Coloc의 유무를 판정할 수 없다. 전체 GWAS 효과크기·표본 N·LD의 적절성은 별개다.

## 3. 변이 효과 방향의 중요한 문제 — 기존 lead가 C/G palindromic

동일 G0015 영역의 기존 `GWS_EAS_AIS` 대표 변이 `chr10:105,596,641:G:C`, 기록된 ALT β `−0.0762`, p `6.214×10⁻¹¹`는 **C/G 회문형 SNP**다. 2026-11-22 작성한 최초 검사에서 이를 일반 비회문형 SNP와 동일하게 계산하는 **테스트 설계 오류**를 발견, 원본 읽기 전에 수정했다. 코드는 현재:
- 실제 출판물에 같은 위치·표기 대립유전자 쌍·수치가 있더라도 **수치 일치 여부만 검토**, strand/ALT 방향의 독립 검증은 **HOLD**.
- 이 변이가 수치상 같더라도 **공식 GWAS source identity를 독립적으로 승인하는 SNP 개수에 0을 가산**. 실제 코호트 승인에는 여러 개의 **비회문형** SNP에서 원본과 archived signed β/p·FASTA·빌드 비교 및 독립 study metadata가 필요하다.
- EUR/GWAS와 EAS/GWAS는 동일 컨소시엄의 *서로 다른 ancestry summary*일 수 있으며, 이를 BBJ-excluded **독립 재현**이라고 부르지 않는다.

## 4. 실행 모델, 원본 신뢰도, 미충족 승인 기준

- `scripts/is/audit_p3a_chr10_original_publisher_gwas.py` + `.github/workflows/is-p3a-chr10-official-gwas.yml`: EAS/EUR **두 원본을 병렬 작업**으로 처리, 원본 압축 파일 전체 SHA-256/MD5/bytes 검증, GRCh37 사전 지정 chr10 2.05Mb context 원본 행 필드/SE/MAF/lead 비교, 실제 36개 그룹에 **공통 GWAS-region 기술적 QC** 연결. raw SNP 파일은 임시 디스크에서 사용한 뒤 폐기. 
- `scripts/is/join_p3a_chr10_publisher_source.py`: 두 GWAS 결과를 합칠 때 `1 genomic region, 36 assays, 20 genes`를 유지하고, 빈 공식 SE/N이나 QTL donor metadata를 임의로 채우거나 `ELIGIBLE`로 승격하는 것을 차단.
- 원본 molecular-QTL 자료: QTD000131 (GTEx artery aorta, catalogue cohort n387), QTD000171 (GTEx brain cortex, n205), QTD000434 (ROSMAP DLPFC, n560). **Study-level catalog numbers are contextual, not approved per-assay genotype sample Ns**. QTL source ancestry, matched LD, effects/SE full region SNP-by-SNP, BBJ/GIGASTROKE EAS precise cohort overlap, exact GWAS phenotype Ncases/Ncontrols and case fraction remain unresolved.

**중요한 연구 판단:** 36개 후보에서 Coloc을 생물학적 인과 결과로 만들기 전에 실제 원본 위상·대립유전자·표본 수·LD·원본 QTL과 GWAS 연결을 해결해야 한다. 기술적 GWAS region PASS가 발생하더라도 **신규 P3A coloc.abf 0/36, coloc.susie 0/36**은 유지한다. 후보 eGene/QTL 유의성은 허혈성 뇌졸중 인과 표적 입증이 아니다.

## 5. 2026-11-22 실제 Github Actions 수치 및 원본 대조 결과

**전체 원본 스캔 완료:** [Actions #45020605713](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45020605713), 4/4 Job SUCCESS, 17/17 fail-closed synthetic tests PASS. EUR/EAS 각 발행기관 전체 원본의 MD5/SHA-256/bytes 일치 후 파싱. 사전 정해진 동일 `chr10:104163676–106220351` 영역의 **유효** 변이 행:

| GWAS original | chr10 context valid variant rows | Publisher whole-file | Scientific |
|---|---:|---|---|
| `GCST90104545` EAS | **3,913** | MD5+SHA256+bytes PASS | HOLD |
| `GCST90104540` EUR | **4,397** | MD5+SHA256+bytes PASS | HOLD |

모든 36개 실제 QTL 분석 그룹이 한 chr10 locus를 공유하며, 20개 gene ID와 동일한 source context에 연결된다. 이 숫자는 **36개 각각의 전체 cis overlap SNP 수가 아니며**, 실제 SNP effect별 QTL–GWAS source matching을 완료했다는 뜻도 아니다.

**기존 30개 locus 원본 기록과의 재검산:** [Actions #45021159785](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45021159785) **SUCCESS**, 6/6 regression tests PASS. 이 경량 CI는 원본 compressed GWAS를 다시 다운로드하지 않고 2026-11-22 실제 원본 scan의 검증된 artifacts와 불변 old-branch 180행 GWAS context TSV를 결합했다.

| 이전 `GWS_EAS_AIS` chr10 연구 기록 | 이번 동일 EAS 원본 새 검사 | 결론 |
|---|---|---|
| Context 유효 행 **3,913** | **3,913** | 기술적 일치 |
| `p<5e-8` 행 **69** | **69** | 기술적 일치 |
| 기존 best p `5.584e-11` | log10 오차 ≤0.05 | 기술적 일치 |
| `10:105596641:G:C` ALT β `−0.0762`, p `6.214e-11` | β·p 숫자 비교 일치 | **C/G palindrome 때문에 strand unresolved** |

**통계적 해석:** 위의 서로 다른 몇 가지 요약 값이 일치하는 것은 아카이브와 발행기관 GWAS 간 **높은 기술적 일관성**의 근거이나, full per-variant source matching의 엄격한 승인 요건인 **2개 이상 영역의 20개 이상 독립된 비회문형 allele-signed variant**를 검증한 결과는 아니다. 해당 lead의 원본 SNP 방향은 FASTA, allele frequency 및 alt orientation 추가 검증 전 **승인 불가**.

또한 25개 public QTL BH05 assay의 142,196/59,090/31,313 기록은 **과거 자체 원본 집계**이고, 이번 공식 GWAS 전체 원본에서 25개 assay에 대해 새로 SNP별 조인한 숫자가 아님을 주의한다.

**다음 과학적 작업:** 불변 30-region 원본의 24개 비회문형 `GWS_EAS_AIS` lead 및 (가능하면) 더 많은 사전등록 SNP를 최소 두 개의 유전체 영역에서 공식 EAS publisher 원본에 직접 비교, effect-allele 방향·beta·p값과 관측 SE를 각각 QC. 출처가 통과하더라도 실제 EAS-IS의 정확한 Ncase/Ncontrol/phenotype 정의와 BBJ inclusion, QTL donor N/alleles 및 full cis QTL 통계가 승인을 받을 때까지 Coloc/MR는 HOLD.
