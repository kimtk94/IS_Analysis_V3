# IS G0022 전체 locus AIS SuSiE-RSS 검증 보고서
**기록:** 2026-10-09 (KST)
**작업 브랜치:** `research/is-broad-discovery-20261009`
**과학적 상태:** EXPLORATORY — 논문에서 최종 causal variant / causal gene으로 결론 내릴 수 없음.

## 1. 왜 재분석했는가

기존 G0022(chr12) AIS 분석은 clump index 주변 ±250kb 세 구간을 각각 별도로 fine-mapping한 탐색적 분석임. ALDH2, PTPN11, IFT81 근처에서 국소 credible set이 하나씩 나타나, 동일 큰 locus에서 서로 독립적인 신호인지 검증할 필요가 있었음.

AS에서 발생한 높은 GWAS–EAS LD mismatch는 AS-only 78개 SNP가 AIS canonical에 전혀 포함되지 않는 원자료 marker coverage 차이와 연관되었음. AS를 임의로 삭제하거나 allele을 바꾸지 않고 **AS 두 구간 BLOCKED, AIS 세 구간 EXPLORATORY** 판정 유지. 별도 분석 문서: `docs/IS_G0022_AS_LD_MISMATCH_ROOT_CAUSE_AUDIT_V2.md`.

## 2. Full-locus 자료 구축과 사용된 수치

- GWAS: GIGASTROKE `GCST90104545`, East Asian AIS, GRCh37.
- Reference: 1000 Genomes Phase 3 EAS, **504명**, G0022 확장 regional genotype.
- 범위: G0022 약 4.05 Mb chr12.
- GWAS/reference allele-matched non-palindromic 분석 후보 **5,060개**.
- 추가 reference MAF ≥0.01 등 QC 후 signed genotype 기반 SNP **4,649개**. MAF 기준에서 411개 제외; 세 원래 AIS clump lead 모두 포함.
- LD matrix: 4,649×4,649 full signed Pearson correlation, 173MB row-major float64, 계산 결과 대칭성 오차 `1.72e-15`; 검증 SHA-256은 결과 JSON에 기록.
- **LD rank 상한 503** (n_ref=504). Reference SNP 수가 reference person 수를 크게 초과하므로 샘플 수/finite-LD uncertainty에 특히 주의.
- Case/control study-wide 근사 유효 표본 수 `N_eff=70,474` 사용; 정확한 SNP별 N_eff 아님.
- 구현: `susieR 0.14.2`, `susie_rss(..., L=8, max_iter=120, estimate_residual_variance=FALSE)`; CS는 반드시 `susie_get_cs(fit, Xcorr=R, coverage=0.95, min_abs_corr=0.5)`로 LD purity 필터 적용.

## 3. Full-locus 실제 결과

| 항목 | 값 |
|---|---:|
| SuSiE 수렴 | Yes |
| 전체 credible set 수 | **1** |
| Purity-filtered 95% CS 안의 SNP 수 | **4** |
| 최고 PIP | **0.37484** |
| 최고 PIP SNP | `12:112241766:G:A` |
| Published/validated causal signals | **0** |

4 SNP의 PIP·유전자 몸체 위치상 주석:

| SNP | PIP | 주석상 겹치는 유전자 |
|---|---:|---|
| `12:112241766:G:A` | 0.37484 | ALDH2 |
| `12:112468206:C:T` | 0.29776 | NAA25 |
| `12:112168009:G:A` | 0.27247 | ACAD10 |
| `12:112736118:A:G` | 0.03778 | HECTD4 |

GENCODE v19 GRCh37 **gene body overlap ≠ causal gene assignment**. Distal eQTL and regulatory targets may differ. PIP under selected model is not per-study replication evidence.

## 4. 기존 ±250kb pilot와 달라진 해석

- ALDH2 근처 lead `12:112241766:G:A`: local PIP≈0.397 → full PIP≈**0.375**.
- PTPN11 근처 lead `12:112930475:T:C`: local PIP≈0.984 → full PIP≈**0.0042**.
- IFT81 근처 lead `12:110675363:C:T`: local PIP≈0.334 → full PIP≈**0.0018**.
- AIS local 3개 분석에서 각각 CS 1개였으나 **full-locus joint fit은 CS 1개**. 각 local CS는 독립 causal signal 3개로 해석하면 안 됨.

## 5. Full-locus 모델 sensitivity (실제 계산)

| 민감도 설정 | Convergence | 95% purity-filtered CS | SNP 수 | 최고 PIP |
|---|---|---:|---:|---:|
| Baseline N_eff≈70,474, L=8 | PASS | 1 | 4 | 0.37484 |
| 낮은 N 가정 20,000, L=8 | PASS | 1 | 4 | 0.37413 |
| N_eff≈70,474, L=3 | PASS | 1 | 4 | 0.37763 |
| N_eff≈70,474, L=8, LD shrink 1% | PASS | 1 | 4 | 0.37429 |

1% diagonal LD shrink는 단순 정규화 가정이며, **정식 finite-reference LD correction 아님**. 최신 susieR 문서의 `R_finite` 및 `R_mismatch` 같은 보정 접근을 격리된 환경에서 추가 검증해야 함. `susieR 0.14.2` 구버전 분석을 확정이라고 보고하지 않음.

## 6. Legacy molecular-QTL provenance (BBJ; AIS 재사용 금지)

기존 `BBJ_IS_L003` + GTEx v8 coloc ABF 자료 감사:
- ALDH2: 16개 조직, best PP.H4=**0.07904** (Artery - Coronary).
- PTPN11: 16개 조직, best PP.H4=**0.11315** (Brain - Nucleus accumbens).
- IFT81: 이전 GTEx coloc 해당 유전자 검사 **0개**.
- ATP2A2: 이전 GTEx coloc 해당 유전자 검사 **0개**.

AIS full locus의 새로운 eQTL/sQTL/pQTL coloc 완료 **0**. 기존 BBJ PP.H4를 AIS에 이전하지 않음. 음성 결과와 검사하지 않은 대상을 엄격히 구분.

## 7. 재현 코드 및 파일

코드:
```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/is/prepare_g0022_full_locus_signed_ld.py
python3 scripts/is/prepare_g0022_full_locus_signed_ld.py --execute
Rscript scripts/is/run_g0022_full_locus_ais_susie.R
Rscript scripts/is/run_g0022_full_locus_sensitivity.R
python3 scripts/is/audit_g0022_full_locus_gene_qtl.py
python3 scripts/is/audit_g0022_full_locus_publication_gate.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

원본을 변경하지 않는 분석 폴더:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/`

핵심 파일:
- `FULL_LOCUS_INPUT_QC.json`
- `ld.rowmajor.f64` (약 173MB; raw genotype-derived LD)
- `variants.tsv`, `G0022_FULL_AIS_PIP.tsv`
- `G0022_FULL_AIS_SUSIE_SUMMARY.json`
- `G0022_FULL_AIS_SENSITIVITY.tsv`
- `G0022_FULL_CS_GENCODE_NEAREST_GENES.tsv`
- `G0022_CANDIDATE_BBJ_QTL_PROVENANCE.tsv`
- `G0022_FULL_AIS_PUBLICATION_GATE_SUMMARY.json`

## 8. 다음 단계 / publication gates

1. **정확한 per-variant effective sample size** 및 meta-study SNP-specific cohort inclusion 확인.
2. Reference n=504에서 발생할 수 있는 finite-LD sampling uncertainty를 모델링하고, 더 큰 ancestry-matched EAS genotype reference를 구해 비교.
3. GTEx v8와 human arterial/brain 및 vascular single-cell eQTL/sQTL/pQTL의 원본 요약통계를 확보하여 full-locus GWAS와 같은 allele·assembly·SNP universe에서 fine-mapping과 coloc 재실행.
4. ALDH2/ACAD10/NAA25/HECTD4뿐 아니라 전체 locus distal cis-regulatory genes 후보 유지.
5. EUR AIS GWAS를 독립 replication으로 주장하지 말고 cohort overlap과 EUR-matched LD를 명시적으로 처리.

**논문 작성 상태:** full-locus 증거는 기술적 탐색/QC 및 Results *candidate evidence*만 가능. Validated independent causal variant / gene / therapeutic target 주장 금지.

## 참고 방법론
- GIGASTROKE source: https://www.nature.com/articles/s41586-022-05165-3
- susieR summary-statistics fine-mapping: https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html
- finite external LD mismatch: https://stephenslab.github.io/susieR/articles/rss_mismatch.html
