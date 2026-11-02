# Ocular TM DEX — 2026-11-02 thesis/publication readiness and immutable source index

**Decision:** computational master’s thesis Methods/Results **GO for a carefully bounded working draft**; journal submission **HOLD**, without demanding investigator contact, expensive raw-data download or false causal interpretation.

## Reproducible evidence index (actual completed GitHub Actions, not plans)

| Evidence tier | Actual source CI / archived artifact | Reproducible objective | Current judgment |
|---|---|---|---|
| Primary GSE313736 paired model | [#41074586446](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446), [actual aggregate freeze #41383491819](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41383491819) | 57 paired donors; 16,050 genes tested; 13,717 FDR significant; high-effect internally defined **83 (67↑/16↓)** | **GO same-source RNA association only** |
| BeamA classification threshold | [#41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259) | same 57 donors, confidence ≥0.75 vs ≥0.50, 81/83 high-effect shared among score-0.75 internal set | **GO cutoff sensitivity** |
| Alternative cell annotation | [#41207583296](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296) | same **55 paired donors**, BeamA vs scPred; 13,975 matched tested genes, 12,634 sign-concordant | **GO within-cohort annotation sensitivity** |
| Effect/abundance stratification | [#41389639386](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41389639386), [Figure/TSV artifact #13116138830](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41389639386/artifacts/13116138830) | weak baseline abs(log2FC)<0.25 genes' direction stability just **78.258% among 5,349** in fixed-55 cross-annotation; full 16-row strata | **GO exploratory sensitivity, not replication** |
| Three-source directional comparison | [#41186736618](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41186736618) | candidate **MYOC↑/TGFBI↓/MMP2↓** in original, independent Mehrotra published 100nM 48h, and GSE301517 500nM 7d; different estimands | **GO signed descriptive synthesis** |
| Strict 7d strain audit | [#41423487666](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41423487666), [per-strain archive #13130831706](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41423487666/artifacts/13130831706) | three core genes same direction in **TM62/TM70/TM72**, 3/3 leave-one-strain-out; EMT 3-study aggregate **124/194**, strict all-3-strain **90/194** | **GO directional specificity; no 7d FDR** |
| Original library confounding / cell selection | [#41430965326](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41430965326), [design/weight artifact #13133928461](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41430965326/artifacts/13133928461) | nine source paired pools, 18 condition-specific 10x libraries, unrestricted library+treatment **19 columns/rank18**; pool selection +0.868638 pp within +1.537701 pp weight = 2.406339 pp | **GO identifiability demonstration; HOLD drug-only causality** |
| Main Figures 1–3 (original) | [#41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011), [PNG+SVG+S1 archive #13045363678](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011/artifacts/13045363678) | original signed gene heatmap, annotation effects, EMT mean study directions | **GO asset/source data; final visual editorial review HOLD** |
| **Final four-artifact machine QA** | [#41432752451](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41432752451), [cross-source audit manifest #13134578619](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41432752451/artifacts/13134578619) | 3 PNG CRC + 3 SVG XML + original source SHA pins + 16-gene S1 + all four source stats + updated abstract | **PASS 2/2 CI jobs** |

### Figure/caption crosswalk (latest scientific claims)

| Figure or table | Cite in Results | Required caption limitation | Disposition |
|---|---|---|---|
| Fig 1 main cross-study candidates | Results 3.2 | Three settings, **500nM 48h 57 donors**, published **100nM 48h**, **500nM 7d 3 TM strains**; CPM log2ratio ≠ edgeR log2FC; no direct dose/time effect test | **GO** existing valid Fig1 |
| Fig 2 main cell annotation | Results 3.3–3.4 | 57 and separately matched 55 same-source donor sets; annotations share biological donors and 10x libraries; ERFE filtered under scPred | **GO** existing valid Fig2 |
| Fig 3 main EMT | Results 3.5 | **124/194** is **three aggregate study-level direction signs**, NOT 3-strain unanimity or pathway inference | **GO** existing Fig3 only with explicit supplementary 90/194 clarification |
| New supplementary Fig — effect and expression strata | Results 3.4a | Baseline-defined bins, 78.258% among 5,349 low-effect genes, not independent validation | **GO** existing #13116138830 figure |
| Supplementary 7d strain gene/EMT table | Results 3.5a | Core genes 3/3 strain signs, strict **90/194** includes 3-study mean and **all 3 strains**; 3 strains are not proven distinct donors | **GO** existing #13130831706 table |
| Supplementary design/proportion table | Methods 2.3a, Results 3.4b, Discussion 4.3a | **19 columns/rank18** does not quantify *fraction* drug-vs-batch. Pool weight term +1.538pp is descriptive composition | **GO** #13133928461 |
| Supplementary Table S1 | Fig 2/results 3.4 | 16 legacy candidates; cross-study estimates present only for 4-core where verified, untested ERFE must not be assigned FDR | **GO** existing original S1 |

## Scope of novelty and submission decision

This work does **not** discover MYOC, ECM changes, a novel drug pathway or a patient glaucoma risk score. Those prior observations and related high-impact works constrain novelty. The differentiated contribution is **fully sourced donor/pool-aware reanalysis, within-cohort label/threshold robustness by effect size, independent published-study strain-level sign scrutiny, and exact library-treatment identifiability and capture-weight accounting**. These are methodological/reporting contributions.

**THESIS GO with caveats:** finalize a computational reproducibility/sensitivity thesis, cite source works and acknowledge unresolved causal and original-80 source limitations. The manuscript currently exists as a working draft; no claim of journal publication or journal-ready peer-review should be made.

**JOURNAL SUBMISSION HOLD** until these two immediate human/editorial checkpoints are closed:

1. **Figure editorial pass:** visually inspect each final PNG/SVG and supplement at expected print size; check text collisions, font sizes, symbol labels, legends and distinguish 124 from 90, since CRC/XML parsing cannot detect visual overlaps. Verify final captions and supplement numbering in the assembled document.
2. **Full manuscript scientific/editorial pass:** confirm bibliography and all source citations, full methods and thresholds, inferential vs descriptive FDR distinctions, external strain provenance/unknown independence, no misleading exposure-time or clinical phenotype claims, and journal formatting. The paper's novelty relative to the 2022–2026 source literature must remain appropriately limited.

A **separate, optional experimental dependency** exists only for **causal/functional** claims: condition-balanced independently randomized 10x/culture processing plus protein/ECM/outflow experiments. Absence of these experiments should restrict claims, not block a limited computational thesis.

**Non-blocking HOLD:** original investigators' 63↑/17↓ complete DEG80 identity table. The original authors will not be contacted as part of this workflow. No exact author80-versus-our83 identity overlap has been established.

**Storage/retention caveat:** GitHub Actions artifacts are subject to expiry. The exact CI IDs, source-file SHA256 and source code are tracked in Git; original Figure PNG/SVG, supplementary PNG and per-gene TSV remain in their linked run archives until retention and should be copied to durable project storage before expiry if needed for submission. This evidence index does not assert such a copy has already happened.
