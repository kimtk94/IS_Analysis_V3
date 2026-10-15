# 석사 논문 Figure 및 Evidence Freeze — IS / CKD

> **2026-10-16 (Asia/Singapore) | SCIENCE CLAIMS NOT YET APPROVED**  
> 목표는 원홍희 교수님 연구 문법과 Nature/Nature Medicine급 **관련 workflow**를 충족하는 것이다. 새 분석법, 신규 후보 발굴 및 controlled data 획득을 필수 조건으로 확대하지 않는다.

## 0. 출처 상태를 3단계로 분리한다

1. **SERVER FILE VERIFIED:** 2026-10-15 17:22 UTC, GitHub Actions [read-only SSH inventory 38776762164](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/38776762164)에서 서버 파일 10개 존재, 크기/행수/SHA256을 확인했다. SHA 확인된 파일의 **존재만** 과학적 검증을 의미하지 않는다.
2. **ACTUAL SOURCE FROZEN:** 이어진 [run 38777071191](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/38777071191)에서 원본 네 파일을 암호화된 SSH로 *임시 GitHub Runner*에 읽어와 원본 SHA256 4개 모두 검증하고 기존 `audit_thesis_evidence_freeze.py` 검증을 통과시켰다. 서버 데이터 변경 없음; 연구 분석 재계산 없음.
3. **NUMERIC SCIENCE SNAPSHOT:** 아래 CKD 유전자별 수치는 2026-10-05 **Google Drive 후보별 MD**에서 추출한 값이다. 최신 서버 통합 TSV에 숫자가 존재함은 확인했으나 **9개 후보의 항목별 값까지 서로 비교해 확인한 것은 아니다.** 이 단계는 해석 초안이며 재현한 연구 결론이 아니다.

## 1. 서버 출처 / Provenance (SHA256)

| 데이터 | 행수 | SHA256 | 출처 |
|---|---:|---|---|
| IS source gene universe | 2,225 | `591cc54583022cc7e09626cf6343b06a1534928f4d584d90e70bee43ade2ff94` | `results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2225_GENE_EVIDENCE_MATRIX.tsv` |
| IS source gene–region | 2,425 | `3ad82780744c4ab571d44fbcf8c421860c0ea0719c6bd38c3a32b2bceb3f9560` | 같은 폴더 `IS_2425_GENE_REGION_EVIDENCE.tsv` |
| IS universe manifest | JSON | `2753b22a8a6633bb00302083aeb55ee4e3b2493a4a09f7af9a7ead958d5ece94` | 같은 폴더 `IS_GENE_EVIDENCE_MATRIX_MANIFEST.json` |
| CKD integrated nine-gene table | 9 | `32fa0475461b206c48295833acfc3075bfb24d89c322b5f9be554b1883a2f3d3` | `results/ckd/integrated_evidence/CKD_9GENE_EVIDENCE_MATRIX.tsv` |

추가 server row count: CKD SuSiE 9; kidney Stage3A 9; cell Stage3B 9; EAS anchor MR 18; EAS eGFR coloc 9. 본 점검에 사용된 server Git SHA `e79564dca3af61040f0bab8143fc9144090fe077`, 브랜치 `research/ckd-thesis-20260927`. 이 server worktree를 사용자가 확인하지 않은 채 reset/merge하지 않는다.

## 2. CKD 9개 유전자 결과 해석 매트릭스 (Drive 2026-10-05 snapshot)

모든 숫자는 분석 결과의 사후확률이며, 실험적 인과성 확증이 아니다. 원문에서 **coloc.abf**는 단일 causal signal을 가정하고, **coloc.susie**는 multi-signal LD와 수렴 진단의 영향을 받는다.

| Gene | EUR eGFR MR FDR | EUR ABF H4 | EUR SuSiE max H4 | EAS eGFR Wald P | Kidney annotated dominant cell | 잠정 해석 |
|---|---:|---:|---:|---:|---|---|
| **ACP1** | 2.54e-8 | 0.925 | 0.909 | 0.0333 | Collecting duct principal | EUR molecular support, EAS nominal direction support |
| **CPVL** | 0.00209 | 0.662 | NA | 0.221 | cDC | Regional signal unresolved; immune localization |
| **F12** | 0.000712 | 2.62e-92 | 1.12e-5 | 1.93e-14 | Proximal tubule (low expression) | Strong EAS anchor association but not shared-signal evidence |
| **GSTA1** | 7.06e-7 | 0.816 | 0.770 | 0.546 | Proximal tubule | ABF support, multi-signal weaker |
| **GSTA3** | 1.44e-7 | 0.952 | 0.965 | 0.584 | Loop of Henle (very low expression) | EUR statistical concordance; weak cell specificity and no nominal EAS |
| **HLA-E** | 1.40e-5 | 1.15e-7 | 0.000210 | 0.555 | Vascular endothelial | MHC/LD and nonshared signal; avoid causal claim |
| **INHBC** | 4.42e-18 | 0.951 | 3.44e-14 | 0.0644 | Proximal tubule | Strong ABF/SuSiE discordance; independent kidney pQTL precedent |
| **SDCCAG8** | 8.56e-14 | 0.979 | 0.979 | 0.0872 | Macrophage (broad) | Best internally concordant EUR molecular signal, but EAS only direction support in this snapshot |
| **UMOD** | 1.62e-218 | 5.54e-41 | NA / no eGFR convergence | 2.26e-47 | Loop of Henle | Strong marginal kidney genetics; colocalization **NOT** supported by current reference-LD analysis |

**Important corrections to prevent overclaiming:**
- IS/CKD results should not say `MR significant -> colocalized -> causal target`. Wald associations can reflect LD, pleiotropy, clearance and reverse renal physiology.
- UMOD is biologically well-established but **the tested pQTL and eGFR region may reflect distinct signals**; no forced-SuSiE completion.
- INHBC's Hirohama kidney pQTL support is *external literature evidence*; it does not resolve the internal ABF-vs-SuSiE conflict.
- EAS p<0.05 at one anchor SNP is **not** full EAS regional colocalization nor independent Korean KoGES longitudinal validation.
- HPA dominant cell reflects **expression/localization**, not cell-type-specific QTL causality.
- Gene novelty must be compared against prior CKD studies: Si et al. (2024), Hirohama et al. (2025). Do not describe UMOD/SDCCAG8/INHBC/HLA-E as wholly novel.

**Source:** CKD candidate notes under [Drive CKD_MASTER](https://drive.google.com/file/d/1fP51Xte0oBWnAvGUsajszZaxWEi1KCuk/view), generated 2026-10-05. Numerical entries will be promoted to latest-source VERIFIED only after matched-item readback from the SHA-pinned server TSV and relevant raw/provenance files.

## 3. IS manuscript workflow freeze

Existing **2,225 gene and 2,425 gene–region** source matrices have been authenticated by SHA and schema, not post-hoc filtered to the four historical named genes. Public nominal QTL coverage had 149 rows with targeted nominal source records and separately 629 not run, 91 targeted-slice no nominal rows, 1,356 outside that EAS pilot; no direct causal gene proof arises from counts.

IS main mechanistic branches remain:
- FGF5: GWAS/eQTL/protein, but qualify effective sample N and source LD.
- ALDH2 rs671: EAS coding/protein/metabolic; do not force eQTL mediation.
- SH3PXD2A: immune/vascular candidate needing matched-cell annotation.
- COL4A2: mural/perivascular candidate; maintain COL4A1 as competing locus gene.

The chr4 / chr10 GWAS-QTL multitrack SVGs from the broad research branch are already available ([original figure report](https://github.com/kimtk94/IS_Analysis_V3/blob/research/is-broad-discovery-20261009/docs/IS_PRIORITY4_GWAS_QTL_MULTITRACK_FIGURES_20261011.md)), but must not be labeled formal SuSiE coloc or conditional LD.

## 4. Final Figure mapping — no new methods

| Figure | IS | CKD | Source / acceptance gate |
|---|---|---|---|
| **Fig. 1** | BBJ/EAS → full eligible 2,225 → regulatory/cell narrowing | UKB-PPP instruments → 9 candidates → tissue/EAS | original input manifest, attrition denominator, methods and cohort metadata |
| **Fig. 2** | GWAS allele/QC/ancestry regional association, coding ALDH2 | EUR eGFR/CKD MR estimates with FDR and support outcomes | exact harmonized beta/SE/P, uncertainty, no p-value only |
| **Fig. 3** | FGF5 and competing locus QTL tracks plus ABF/SuSiE supported vs unavailable | ACP1/GSTA3/SDCCAG8 concordance and INHBC/UMOD disagreement diagnostics | original PP.H3/H4, priors, SNP overlap and LD source; unresolved clearly labelled |
| **Fig. 4** | brain vascular immune/mural cell evidence + prioritized genes | kidney tissue pQTL/eQTL, cell localization and EUR/EAS evidence | expression-only vs QTL vs functional link separate; gene novelty audit |

**Supplementary:** full IS candidate denominator, QTL availability, SNP harmonization, missingness; all CKD 9 candidates with allele/LD/variant counts, MR and regional co-localization, method discrepancies, and source/version hashes.

## 5. Scope-specific stop conditions

- **Finish / STOP analytic expansion** when Figure 1–4 and supplements reproduce from frozen dated outputs, failing/unknown fields are retained, and statistical claims do not exceed the original study data.
- **If BBJ large cohort-matched LD inaccessible:** report credible-set and COJO limits and reference-based sensitivity. Do not imply independent conditional signals.
- **If KoGES participant genotype inaccessible:** do not use `validated in longitudinal Korean cohort` in title/abstract. Public summary-data CKD study can still be completed.
- **No new PRS, Genomic SEM, PheWAS, de novo single-cell atlas, or deep AI scoring** unless the advisor explicitly changes the research question.
- Successful CI is **software QA**, not paper acceptance.

### Next only necessary scientific verification

1. Compare **9 numeric rows from the SHA-pinned latest CKD matrix** against the Drive snapshot above, especially UMOD EUR SuSiE and INHBC discordance.
2. Check Stage2B/Stage2C actual LD provenance, signed-allele orientation, prior sensitivity and any high-H4 signal with differing GWAS source variants.
3. Plot four manuscript figure groups *from frozen analysis outputs*, not from screenshots or approximate values.
4. Advisor review of study question and figure/claim boundaries, then **analysis freeze**.

See [runbook](./THESIS_EVIDENCE_FREEZE_RUNBOOK_20261015.md), [audit](./RESEARCH_SCOPE_COVERAGE_AUDIT_20261015.md), and [Won-guided research direction](./RESEARCH_DIRECTION_WON_HIGH_IMPACT_20261015.md).
