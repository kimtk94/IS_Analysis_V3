# IS Harness P0.4 — 통합 근거 평가표 인수인계

**작성:** 2026-10-18 KST  
**상태:** `CODE_AND_TEST_COMMITTED / SCIENTIFIC_SERVER_EXECUTION_PENDING`

## 목적
P0.1(유의 CALHM2 QTL SNP 3개), P0.2(개인별 EAS/JPT 유전형의 SNP 간 pairwise LD), P0.3(CALHM2·NEURL·FGF5·C4orf22의 동일 기준 비교)을 **하나의 locus/gene evidence scorecard**로 묶는다.

**핵심 규칙:** 입력이 없는 단계를 완료로 간주하지 않는다. 모든 결과는 `DESCRIPTIVE_QTL_GWAS_PUBLIC_REFERENCE_ONLY`이며 원인 유전자 확정, SuSiE credible set, multi-signal coloc 결과가 아니다. 이 파일은 점수 가중치를 임의로 설정하는 '인과 랭킹'이 아니라 검증된 근거와 결측 상태를 분리하는 장부다.

## 실행 입력 및 스크립트
- `scripts/is/build_is_priority4_evidence_scorecard.py`
- `tests/test_is_priority4_evidence_scorecard.py`
- P0.1 ledger: `/srv/is-analysis/results/is/audits/is_calhm2_three_qtl_snp_ledger_20261018_v1/CALHM2_THREE_QTL_SNP_LEDGER.tsv`
- P0.2: `/srv/is-analysis/results/is/audits/is_calhm2_pairwise_ld_20261018_v1/CALHM2_THREE_SNP_PAIRWISE_LD.tsv`
- P0.3: `/srv/is-analysis/results/is/audits/is_priority4_qtl_gene_comparison_20261018_v1/IS_QTL_GENE_COMPARATOR_SUMMARY.tsv`
- 새 산출: `/srv/is-analysis/results/is/audits/is_priority4_evidence_scorecard_20261018_v1/`

P0.3가 없다면 스크립트는 오류로 중단한다. P0.1/P0.2가 없으면 해당 파일이 없는 상태로 P0.4 장부를 만들 수 있지만, CALHM2의 검증 단계는 `NOT_EXECUTED_OR_MISSING`으로 유지한다. 이 경우에는 연구상 통합 분석 완료가 **아니다**.

```bash
set -euo pipefail
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
git status --short
git fetch origin research/is-broad-discovery-20261009
# Concurrent worktrees/untracked files를 보존하면서 변경 통합
python3 -m unittest discover -s tests -p 'test_is_priority4_evidence_scorecard.py' -v
python3 scripts/is/build_is_priority4_evidence_scorecard.py \
  --comparison /srv/is-analysis/results/is/audits/is_priority4_qtl_gene_comparison_20261018_v1/IS_QTL_GENE_COMPARATOR_SUMMARY.tsv \
  --ledger /srv/is-analysis/results/is/audits/is_calhm2_three_qtl_snp_ledger_20261018_v1/CALHM2_THREE_QTL_SNP_LEDGER.tsv \
  --pairwise /srv/is-analysis/results/is/audits/is_calhm2_pairwise_ld_20261018_v1/CALHM2_THREE_SNP_PAIRWISE_LD.tsv \
  --out /srv/is-analysis/results/is/audits/is_priority4_evidence_scorecard_20261018_v1
python3 -m unittest discover -s tests -p 'test_is_*.py' -q
```

## 출력·완료 판정
- `IS_PRIORITY4_INTEGRATED_EVIDENCE_SCORECARD.tsv`: 4개 유전자 × QTL/GWAS 명목근거, P0.1/P0.2/P0.3 상태, CALHM2 3 SNP가 EAS/JPT에서 서로 모두 `r²≥0.8`인지의 상태값, formal causal inference block.
- `IS_PRIORITY4_INTEGRATED_EVIDENCE_MANIFEST.json`: 모든 소스 경로·존재 여부·SHA256, 입력 유전자의 명시적인 QC cardinality.
- `P0.4=COMPLETE`는 세 입력이 **모두 실제 실행됐고**, 상태가 확인됐으며, 각 LD pair가 근거를 갖추고, 테스트가 통과하고, 산출물이 커밋된 경우에만 허용.
- P0.1/P0.2가 빠졌다면 분석 상태는 **PARTIAL_LEDGER_ONLY**, 원인 유전자 선정도 **NOT_ESTABLISHED**.

## 별도 잔여 과제
1. CALHM2 3 SNP의 QTL·GWAS 효과방향과 EAS/JPT pairwise (3\times3) LD heatmap.
2. NEURL 비교 및 FGF5 양성 비교군 결과 해석.
3. QTL effective sample N과 donor-specific LD 검증(여전히 없음).
4. 기존 2,225 discovery 유전자 및 646 coloc, 120 ABF 원본 보존.
5. 분석 코드가 GitHub에 있는 것과 **연구 서버에서 실제로 실행됐다는 것**을 구분할 것.

**실행 제약:** Remote Desktop Commander 월간 호출 제한으로 이번 턴 서버 실행·전체 회귀 테스트는 수행하지 않았다. GitHub 코드 작성만 확인했다.
