# P3A 36개 QTL 후보의 GWAS 출처 검증 — 공식 EAS 원본 24개 다염색체 SNP 일괄 재현

**확인일: 2026-11-22 | MASTER DEGREE / IS 연구 | 상태: 원본 GWAS 수치 재현 PASS, 신규 P3A Coloc HOLD**

## 핵심 확인

불변 연구 커밋 `0bf646fd6f29322d17da878253f37def24c01c74`의 **30개 GWAS 영역×6개 outcome(180행)** 원본에서 GIGASTROKE `GWS_EAS_AIS`의 지역별 대표 SNP 30개를 회수했다. **24개는 비회문형 SNP**, 4개는 A/T 또는 C/G 회문형, 2개는 해당 outcome에 대표 변이 미기록이다. 24개 비회문형 대표 SNP는 **16개 염색체**에 분포한다.

공식 GWAS Catalog 발행기관 원본 `GCST90104545_buildGRCh37.tsv.gz` **전체 260,689,208바이트**를 직접 다운로드하여 발행기관 MD5 `904e12e33ec3c2eeefe4dc93bec4995c`, 과거에 검증한 SHA-256 `c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3`, 크기를 확인한 다음 **6,725,216개 원본 변이 행**을 읽었다.

- **24/24** 비회문형 실제 발행기관 변이 위치·대립유전자 쌍 매칭.
- **24/24** 기존 GWS EAS AIS record의 ALT 방향 β, p, EAF가 사전 지정된 수치 허용차 내에서 재현됨.
- 재현된 SNP의 염색체는 **16개**, 그러나 **16개 독립 코호트나 24개의 독립 유전학적 검증을 뜻하지 않는다**.
- **4개 회문형 SNP**는 strand/FASTA/allele 주의가 필요해 수치가 일치하더라도 source-identity 근거에서 제외. 이 중 `IS_XDATA_G0015`의 대표 변이 `10:105596641:G:C`도 포함.
- **2개 대표 SNP 미기록**은 이번 24개 대조에서 제외됐고, 음성 GWAS 또는 no-association으로 판정하지 않음.

## 실제 CI 검증

| CI 근거 | 결과 |
|---|---|
| [원본 260MB GWAS 직접 재조회·24개 signed SNP 비교 #45022078752](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45022078752) | 2/2 jobs SUCCESS; 안전성 테스트 13/13 PASS |
| [24개 전체 결과 독립 검토 #45022398163](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45022398163) | SUCCESS; 8/8 tests PASS; 24/24 signed-β/p/EAF 재현, 불일치 0 |
| [P3A 36개 복원 #45019987189](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45019987189) | 실제 36 gene×QTL groups / 20 genes PASS, 신규 Coloc 실행 0 |

**수치 비교 기준 (원본 비교용 QC):** ALT 방향으로 원본 β를 뒤집어 맞춘 후 `abs(β 차이) ≤ max(0.00006, 0.015×abs(β_archive))`, `abs(log10(p_original)-log10(p_archive)) ≤ 0.05`, `abs(EAF_original_ALT-EAF_archive_ALT) ≤ 0.025`. 위는 서로 다른 포맷의 저장·반올림 차이를 점검하기 위한 원본 서명 일치 기준이며, Coloc 유의성이나 sample ascertainment 검증 기준이 아니다.

**중요한 제한:** 과거 GWAS 영역 표와 공식 원본은 **동일 GIGASTROKE 데이터 릴리스에서 비롯됐을 수 있다**. 본 결과는 `external study replication`이 아닌 **original-source numerical reproduction**이다. 24개의 signed-allele/key/β/p/EAF 일치 결과는 source consistency의 기술적 증거지만, 검증된 별도 BBJ 원본, source-specific effective Ncase/Ncontrol, genomic FASTA REF, 전 cis 구간의 SNP-by-SNP matched QTL/GWAS effect 증명을 대체하지 않는다.

## 실제 30개 provisional regions 전체 감사 기록

| Region | 역사적 EAS AIS 대표 SNP (GRCh37) | 독립 원본 감사 결과 |
|---|---|---|
| IS_XDATA_G0001 | 1:10796866:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0002 | 1:113042822:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0003 | 2:2645097:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0004 | 2:26929282:G:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0005 | 2:161641655:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0006 | 3:76418451:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0007 | 4:81181072:T:A | 회문형 HOLD |
| IS_XDATA_G0008 | 4:111704199:G:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0009 | — | 미기록 |
| IS_XDATA_G0010 | 7:27316640:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0011 | 7:125956297:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0012 | 9:94337776:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0013 | 9:96105375:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0014 | 10:36122291:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0015 | 10:105596641:G:C | 회문형 HOLD |
| IS_XDATA_G0016 | — | 미기록 |
| IS_XDATA_G0017 | 10:124227624:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0018 | 11:1898664:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0019 | 12:20583021:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0020 | 12:24770775:G:C | 회문형 HOLD |
| IS_XDATA_G0021 | 12:90019178:G:C | 회문형 HOLD |
| IS_XDATA_G0022 | 12:112241766:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0023 | 13:110959643:C:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0024 | 16:81808816:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0025 | 17:7128727:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0026 | 17:43076102:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0027 | 18:42710102:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0028 | 19:10893074:G:A | beta·p·EAF 재현 PASS |
| IS_XDATA_G0029 | 20:50108963:C:T | beta·p·EAF 재현 PASS |
| IS_XDATA_G0030 | 22:17434253:C:T | beta·p·EAF 재현 PASS |

## P3A 후보 유전자 전체에 주는 의미

복원된 **36개 QTL study×gene 분석 그룹 / 20개 gene ID**는 전부 하나의 `IS_XDATA_G0015` 구간에 있다. 이 유전체 영역은 기존 [chr10 원본 검사 #45020605713](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45020605713)에서 실제 발행기관 EAS 원본 유효 SNP 3,913행, EUR 원본 4,397행을 직접 읽었고, 과거 EAS G0015 원본의 3,913행 / p<5e-8 69행 / best p가 [#45021159785](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45021159785)에서 일치했다. 그러나 해당 지역의 대표 변이 자체가 C/G 회문형이므로 이번 비회문형 24개 재현 배치에는 포함되지 않았다.

25개 official public BH05 QTL assay의 full-cis 원본 집계는 **142,196 gene-assay SNP 기록**, BBJ/EAS 양쪽에서 기존 matching이라고 보고된 비회문형 **59,090**, 그중 QTL–EAS 빈도차 ≤0.1 **31,313**이다. 이것들은 **유전자별 검사 기록의 합**이고 유일한 SNP 수가 아니며, 새로 공식 원본과 실제 full-cis 36개 전체를 검증한 수치가 아니다. 공개 BH05 25개와 과거 서버의 source-stat-positive 25개 배정이 정확히 동일한지 미확인이다.

## 다음 과학적 승인 기준

**현재 정식 P3A `coloc.abf` 승인 0/36, 실행 0/36, `coloc.susie` 실행 0/36**을 유지한다. 이후의 주요 작업은 유전자별 하나씩 산발적으로 처리하지 않고 다음을 묶는다.

1. **원본 GWAS ↔ archived 36개 group full-cis 서명 대조:** 실제 chr10 variant key, effect allele, β/p/관측 SE, GRCh37 REF/ALT/FASTA를 전체 36 QTL study×gene 그룹과 연결한다. SNP 후보만 선별해서 H4를 계산하지 않는다.
2. **GWAS/QTL 표본 모형 확인:** 출판된 정확한 `Ncase/Ncontrol`, IS/AIS phenotype definition, BBJ participant overlap; QTD000131/171/434 원본 tissue별 sample N, ancestry, AF, summary SE와 전체 cis window의 유효 tested variants 검증. 사전 임의 N_eff 대입 금지.
3. **실제 `ELIGIBLE` 그룹만 R `coloc.abf` 일괄 실행:** region-wide PP.H0–H4, H3/H4, prior sensitivity 보고. 서명 및 sample-size 승인 없는 그룹은 HOLD로 구분. 다중 신호 `coloc.susie`에는 추가로 SNP-order/ancestry/allele이 일치하는 signed LD와 SuSiE 모델 QC 필요. 20개 후보유전자의 분석은 중복 신호와 여러 조직이 한 locus를 공유한다는 사실을 고려해 정리한다.

**저장 산출물:** 위 두 CI의 Actions artifacts에는 24개 `region_id`별 boolean/status QC 결과가 남으며, 원본 SNP 전체 데이터는 GitHub에 업로드하지 않았다.
