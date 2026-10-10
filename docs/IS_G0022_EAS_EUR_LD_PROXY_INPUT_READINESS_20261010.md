# IS G0022 동아시아–유럽인 GWAS 변이 커버리지 및 rs671 EAS LD proxy 감사

**2026-10-10 KST | 검증 수준: SOURCE_QC_VERIFIED / MULTI_ANCESTRY_FINE_MAPPING_BLOCKED**

## 1. 질문 및 입력

동아시아 AIS full-locus 95% credible set 4개가 EUR AIS에도 있고, rs671과 높은 EAS LD를 가지는 공통 SNP이 EUR summary statistics에서 관측되는지 확인했다. 두 결과는 유럽인에서의 rs671 LD를 계산한 결과가 아니라 **정확한 position + unordered allele pair 기반 source coverage 검사**이다.

- EAS AIS GWAS: GIGASTROKE GCST90104545, 기존 allele-verified 4,649 SNP variants.tsv와 signed EAS LD matrix 4,649×4,649 (1000G EAS 504명; rank≤503).
- EUR AIS GWAS: GIGASTROKE GCST90104540 GRCh37 original /srv/is-analysis/data/is/reference/gigastroke/additional/GCST90104540_buildGRCh37.tsv.gz, 7,482,032 rows.
- EUR에는 독립적인 정확한 SNP REF/ALT column이 없으므로 기존 EAS verified REF/ALT에 대해 (effect,other) pair가 직·역 매칭되는 것만 인정한다. 역일치 시 beta=-beta, EAF=1-EAF. Genome reference와 GWAS meta allele QC를 완전히 대체하지 않음.
- EAS 변이는 좌표 순서와 매트릭스 row 순서를 보존. rs671 row index=2356 (zero-index). 파일 크기 및 LD diagonal 1.0 확인.
- 각 source, GWAS summary, outputs는 canonical 원본을 변경하지 않고 read-only로 분석.

## 2. 4,649개 vs EUR AIS 실제 교집합 (실측)

| metric | observed |
|---|---:|
| EUR AIS original complete GWAS variants | 7,482,032 |
| EUR chr12 rows | 360,431 |
| EUR G0022 EAS LD analysis window (min–max SNP span) rows | 6,954 |
| EAS full-locus reference-QC SNP | 4,649 |
| EUR GWS position+allele matched | **3,303 (71.0475%)** |
| EUR source exact position absent | **1,346** |
| position present but alleles unmatched | 0 |
| multiple exact matches (ambiguous) | 0 |
| EAS 95% CS SNPs represented | **0/4** |

4/4 absent SNP identifiers (GRCh37):
- ACAD10 rs11066015: 12:112168009:G:A
- ALDH2 rs671: 12:112241766:G:A
- NAA25 gene-body annotation: 12:112468206:C:T
- HECTD4 gene-body annotation: 12:112736118:A:G

**원본에서 SNP가 없다는 것은 해당 EUR 코호트가 이 대립유전자를 가지고 있지 않다는 증거가 아니다.** 원본 QC, 빈도 기준, 임퓨테이션 품질, 각 메타분석 참여코호트 coverage가 여전히 미상이다.

## 3. rs671 기준 EAS signed LD와 EUR source availability

rs671과의 EAS genotype reference LD를 이용해 threshold 이상인 marker set을 만든 후 그 각 variant에 대한 EUR GWAS 원본 존재 여부를 검사한다. EAS r²를 EUR LD처럼 대체하지 않는다.

| r²(rs671, proxy) in 1000G EAS n504 | EAS SNP denominator | EUR exact-match source rows | available (%) |
|---|---:|---:|---:|
| >= 0.95 | 5 | 0 | 0.0 |
| >= 0.90 | 7 | 0 | 0.0 |
| >= 0.80 | 9 | 0 | 0.0 |
| >= 0.50 | **19** | **0** | **0.0** |
| >= 0.20 | 125 | 75 | 60.0 |
| >= 0.10 | 622 | 294 | 47.27 |

**새 핵심 관찰:** 전체 원본 중 EUR SNP 교집합이 71.0%라도, rs671을 EAS에서 강하게 태깅하는 SNP 19개는 EUR AIS source에서 전부 미수록이다. EUR을 추가했다는 사실만으로 rs671 signal resolution 개선을 주장할 수 없다. 고LD SNP 부재가 유전적 희소성 때문일 개연성이 있지만, source-specific filter 감사 없이 인과 원인 확정 불가.

## 4. 다민족 fine-mapping의 판단

**현재 BLOCKED.** 이유는 (a) EUR source CS 0/4 및 strongest EAS LD proxy 0/19, (b) ancestry-matched signed EUR LD 미준비, (c) GIGASTROKE meta-study cohort overlap/variant-specific N 미확인, (d) 1000G EAS n504 vs SNP4,649 저랭크 LD, (e) true causal signal과 ancestry-specific genetic architecture 미규명이다.

- MultiSuSiE (Rossen et al., Nat Genet 2026)는 ancestry 간 effect correlation을 허용하지만 cohort별 signed LD 및 GWAS summary input이 필요. 논문은 가능한 in-sample LD를 권고하고, 외부 LD라면 각 GWAS target 표본의 10% 정도 이상 규모의 reference를 대안으로 언급한다. 본 연구 n504는 G0022 GWAS sample-level N/MAF 구조와 불일치할 가능성이 높아 finite LD uncertainty가 중요하다. **N_eff≈70,474는 GWAS case/control의 대체근사로 10% 조건을 단순 검증하는 직접 분모로 삼지 않는다.**
- MultiSuSiE/SuSiEx의 공통 causal-variant 가정은 rs671 EUR rare region에서 신중해야 한다. MESuSiE는 ancestry-specific causal variant 허용하지만, 데이터 누락이나 EUR LD 부재를 해결해주는 도구가 아니다.
- 후속 연구는 EAS 중심 robustness/fine mapping, JCTF/MAEEA complete cis-QTL 접근, 독립 EAS AIS cohort 검증에 우선 배정하고 EUR 통합은 variant/LD provenance gate 이후 reconsider.
- 3,303개 EUR overlap의 실제 효과 방향/인구간 replication은 **분석하지 않았음**. 교집합 71%라는 수치를 shared effect replication으로 쓰지 않는다.
- 치료표적 방향/알코올 mediating effect/MR causal effect도 새롭게 추정하지 않았다.

## 5. 완전 재현(서버)

Research code in /srv/is-analysis/worktrees/is-broad-discovery-20261009
- scripts/is/audit_g0022_cross_ancestry_readiness.py
- scripts/is/audit_g0022_rs671_ld_proxy_ancestry_coverage.py
- tests/test_is_g0022_cross_ancestry_readiness.py
- tests/test_is_g0022_rs671_ld_proxy_coverage.py

Output root:
  /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/cross_ancestry_readiness_v1

Output files:
- G0022_EAS_EUR_SOURCE_VARIANT_COVERAGE.tsv (4,649 EAS rows; source matched only)
- G0022_EAS_EUR_READINESS_SUMMARY.json (full original source SHA256)
- G0022_RS671_EAS_LD_PROXY_EUR_COVERAGE.tsv (4,649 LD row annotations)
- G0022_RS671_EAS_LD_PROXY_SUMMARY.json (6 threshold cohorts)

Default exact input paths above are passed using command-line flags, not embedded in source Python.
Verification:
  python3 -m unittest discover -s tests -p 'test_is_g0022_cross_ancestry_readiness.py' -v
  python3 -m unittest discover -s tests -p 'test_is_g0022_rs671_ld_proxy_coverage.py' -v
  bash scripts/codex_smoke_test.sh

## Sources
- GIGASTROKE genome-wide IS loci: https://www.nature.com/articles/s41586-022-05165-3
- MultiSuSiE (2026): https://www.nature.com/articles/s41588-025-02450-5
- MultiSuSiE code: https://github.com/jordanero/MultiSuSiE
- 1000G rs671 ancestry contrast: https://pmc.ncbi.nlm.nih.gov/articles/PMC8284761/
