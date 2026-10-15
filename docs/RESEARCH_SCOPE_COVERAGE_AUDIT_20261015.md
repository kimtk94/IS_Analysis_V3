# IS / CKD — High-impact workflow coverage audit (2026-10-15)

**Overall classification:** GitHub evidence audit DONE; server artifact and scientific-result verification NOT DONE.

**Research scope policy:** Match question-relevant workflow of Mishra (Nature 2022) and Wälchli (Nature 2024) for IS; Hirohama (Nature Medicine 2025) and Sheng (Nature Genetics 2021) for CKD. Align the multi-omic evidence stack and reporting style with Hong-Hee Won laboratory work (Cell Reports Medicine 2023 dyslipidemia; eBioMedicine 2025 problematic alcohol use; Nature Genetics 2026 SingleBrain). Do not add extraneous PRS, genomic SEM, large PheWAS, wet-lab or new cell atlas just to imitate every experiment.

See [research direction](./RESEARCH_DIRECTION_WON_HIGH_IMPACT_20261015.md) and [minimal workflow](./HIGH_IMPACT_WORKFLOW_BASELINE_20261015.md).

## 1. Evidence and safety boundary

This audit used authenticated GitHub repo and branch contents, Actions run status/logs and Drive file discovery. Remote Desktop Commander listed taehoonserver ONLINE, but a monthly usage-limit error occurred before any server command ran. Therefore this document is NOT a live server execution audit; no server file hashes were rechecked. The 2026-10-05 main-branch IS_MASTER and Drive CKD_MASTER are historical snapshots and do not include later research-branch work.

GitHub inspected:
- main: CKD Stage1, Stage2B, Stage2C, Stage3A documentation and IS_MASTER.
- research/is-broad-discovery-20261009: 2026-10-11 2,225-gene readiness, allele, LD, COJO, ABF and multitrack Figure audits.
- research/is-evidence-audit-20261014: branch-reconciliation and evidence-audit documents, both expressly documentation/plan rather than server validation.
- research/ckd-thesis-20260927: Stage4 EAS scripts, KoGES controlled-input contract, 9-gene integrated evidence script.
- main vs research branches: all diverged; DO NOT blindly merge, reset or fast-forward.

Interpretation: CODE PRESENT != RUN COMPLETED != QC VERIFIED != BIOLOGICAL CLAIM SUPPORTED.

## 2. IS workflow gate assessment

| Gate | Evidence from versioned files | Provisional assessment |
|---|---|---|
| IS-01 full discovery universe | 2,225 stable Ensembl IDs; 2,425 gene-region combinations (EAS/Japanese 869 and EUR 1,556). Operational readiness A14, B29, C826, D1356; tiers are NOT biological priority | DOCUMENTED; freeze and inspect actual input manifest |
| IS-02 GWAS region/allele QC | Six EAS/Japanese GWAS; 45 phenotype clusters and 30 coordinate-overlap components. 15/45 clusters and 7/30 merged components have genome-wide significant variants. Among 30 lead positions: 3 exact, 1 swapped, 26 absent from existing regional PVAR | PARTIAL: unmatched leads cannot be promoted |
| IS-03 molecular QTL/ABF | Historical 646 assay records: 30 PASS and 616 MISSING_INPUT. The 30 replayed ABF records represent 14 stable genes, not 30 genes. A successful numerical replay is NOT full causal confirmation | PARTIAL |
| IS-04 credible sets / conditional | Existing EAS504/JPT104/EUR503 LD are reference-based sensitivities. Zero of 30 ABF assays had validated GTEx donor signed QTL LD for publication-ready multi-signal coloc in the 2026-10-11 readiness audit. GCTA COJO had only format-QC inputs; formal run not done, EAS504 below recommended reference size | UNRESOLVED / DATA LIMITED. Do not require untrustworthy COJO/SuSiE PASS |
| IS-05 cell context and figures | Vascular/immune functional framework documented; chr4 FGF5/C4orf22 and chr10 CALHM2/NEURL GWAS-eQTL multitrack figures recorded. They plot marginal GWAS and nominal eQTL, not mediation or conditional causal result | PARTIAL |
| IS-06 manuscript evidence freeze | Gene-level output manifest, allele QC, source-matched inference, full attrition, evidence-status table, result-to-figure trace need readback | SERVER NOT VERIFIED |

**Practical IS conclusion:** The candidate universe is already broad; DO NOT keep extending discovery or discard candidates just for absent QTL source data. Prioritize manuscript tables, qualifying the four illustrative loci FGF5, ALDH2, SH3PXD2A, COL4A2 appropriately. ALDH2 rs671 is a coding/protein/metabolic branch; expression absence cannot falsify it. A BBJ-vs-GIGASTROKE marginal-lead difference is not proof of independent causality.

Relevant tracked sources:
- [full universe audit](https://github.com/kimtk94/IS_Analysis_V3/blob/research/is-broad-discovery-20261009/docs/IS_2225_COMPLETE_UNIVERSE_WORKFLOW_READINESS_AUDIT_20261011.md)
- [formal COJO readiness](https://github.com/kimtk94/IS_Analysis_V3/blob/research/is-broad-discovery-20261009/docs/IS_PRIORITY4_BBJ_EAS_COJO_READINESS_20261011.md)
- [30-ABF multimodal readiness](https://github.com/kimtk94/IS_Analysis_V3/blob/research/is-broad-discovery-20261009/docs/IS_30_ABF_ANCESTRY_MULTISIGNAL_READINESS_AUDIT_20261011.md)
- [branch reconciliation](https://github.com/kimtk94/IS_Analysis_V3/blob/research/is-evidence-audit-20261014/docs/is/IS_BRANCH_RECONCILIATION_20261014.md)

## 3. CKD workflow gate assessment

| Gate | Evidence from versioned code/docs | Provisional assessment |
|---|---|---|
| CKD-01 cis-pQTL and GWAS QC | UKB-PPP ST16 ALT orientation, F>=10, strongest-F primary Wald, EUR eGFRcrea BH FDR, outcome/allele harmonization are specified in docs/ckd_stage1_mr.md | CODE/DOCUMENTED, output not reread |
| CKD-02 regional MR/ABF | Stage2B GRCh38-to-37 liftover, UKB-PPP alleles, coloc prior grid p12=1e-6, 1e-5, 1e-4 specified | CODE/DOCUMENTED, output not reread |
| CKD-03 multi-signal | Stage2C 1KG EUR reference LD, signed orientation, coloc.susie and convergence cap specified. UMOD remains an especially difficult unresolved locus in prior project state | PARTIAL, do not force convergence |
| CKD-04 kidney tissue and cell | Hirohama 2025 supplementary plus kidney meta/tubule/glomerulus eQTL mapping documented, with significant-only and workbook-retrieval caveats. Stage3B cell annotation present in older project snapshot | PARTIAL / need 9-gene readback |
| CKD-05 EAS support | EAS regional matching and integrated evidence scripts confirmed on research/ckd-thesis-20260927 branch, not necessarily merged into main | CODE IN BRANCH, runtime not verified |
| CKD-06 longitudinal KoGES | Controlled-data analysis contract and scripts exist, but nothing confirms actual controlled individual genotype availability or completed longitudinal inference | BLOCKED IF CLAIMING KOGES VALIDATION |
| CKD-07 thesis tables/figures | Integrated 9-gene script enumerates ACP1 CPVL F12 GSTA1 GSTA3 HLA-E INHBC SDCCAG8 UMOD. Complete candidate table and corresponding final figure output not read directly today | SERVER NOT VERIFIED |

**Practical CKD conclusion:** The core high-impact analysis structure already exists in code. Auditing the 9-row gene matrix and actual output provenance is more useful than implementing extra models. Label pre-reported INHBC, UMOD, SDCCAG8 and HLA-E properly. A public-data-only inference paper is viable in scope; KoGES longitudinal replication must not be claimed before individual-level analysis.

## 4. GitHub CI: concrete defect, repair, and verification

**Run with fault**: [38565876304](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/38565876304), 2026-10-14. IS broad discovery regression reported 343 unittest cases, 4 test module import errors and 30 skips. All four import errors traced to missing pandas. This was CI environment configuration, not evidence the four genomic methods were invalid.

**Minimal repair**: Added pandas==2.2.3 to the dependency-install step in research/is-broad-discovery-20261009 / .github/workflows/is-broad-discovery-ci.yml, alongside existing PyYAML==6.0.2, numpy==2.1.3 and requests==2.32.3. Commit [c07eb93](https://github.com/kimtk94/IS_Analysis_V3/commit/c07eb93bfd2f03959d6464ce11fc042a4a03922a). No genomic algorithm or live server setting changed.

**Verified re-run**: [38703317020](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/38703317020), 2026-10-15, GitHub Actions status COMPLETED / SUCCESS. Installation, Python compile, and synthetic IS regression steps all SUCCESS. This validates CI collection/execution, NOT scientific reproducibility of actual GWAS results.

## 5. Minimal remaining work; strictly no analysis sprawl

1. **P0 read-only server evidence inventory** when authorized tool access resumes: git SHA/branch/worktree status; IS candidate universe manifest; CKD 9-gene evidence and QC; exact output timestamps and SHA256. Do not overwrite, git reset, merge or rerun science solely for a cleaner report.
2. **P0 IS paper table freeze**: 2,225 universe/attrition + ancestry/gene/source/provenance and pass/partial/unresolved flags. Four core locus examples are illustrations, not universe definition.
3. **P0 CKD nine-gene evidence freeze**: original result row for Wald/IVW, ABF PP.H3/H4, SuSiE sets/convergence, kidney tissue pQTL/eQTL, cell labels, EAS and prior-paper overlap. Missing rows should be marked MISSING, not treated as negative.
4. **P1 final manuscript figures**: IS study flow, ancestry+locus, molecular evidence, cell annotation; CKD instrument/MR flow, candidate results, representative coloc, kidney tissue/cell/EAS evidence. Supplementary tables preserve full discovery denominator and failed/unresolved tests.
5. **STOP once evidence QC+plots+methods are reproducible** and any unavailable LD/KoGES controlled data is transparently reported as a limitation. Do not add unrelated TWAS, PGS, PheWAS, SEM, new atlas or wet-lab modules.

## 6. Completion and interpretation policy

A stage is PASS only when an original result and QC provenance are present; reference ancestry/build/alleles/outcome match; statistical interpretation and multiple-testing are correct; and a manuscript table or figure can be regenerated. If data are inaccessible, mark BLOCKED and explicitly restrict the claim rather than inflating a sensitivity result. **High journal impact factor does not replace scientific validity.**

**Audit dated 2026-10-15.** GitHub readback and CI re-run are verified; scientific server outputs remain to be inspected.
