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

**JOURNAL SUBMISSION HOLD**: visual inspection of the four current original/fixed figure PNGs has now been completed (see evidence below). One main outstanding human/editorial checkpoint remains for bibliography/source/caption assembly and final formatted thesis proof:

1. **Figure actual visual pass: completed for four PNGs.** Original Fig 1–3 visually GO. Initial effect-stratified supplementary Figure had overlapping \`n=\` labels (VISUAL FAIL despite prior CRC pass). **Final fixed supplement [GitHub Actions #41435268815](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41435268815) 2/2 SUCCESS, [artifact #13135314010](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41435268815/artifacts/13135314010)**; re-opened actual PNG visual GO with n tables and footer separated. Both underlying statistical TSVs are **byte-for-byte identical** old versus fixed. Physical print-size and assembled thesis PDF layout remain downstream editorial checks.
2. **Full manuscript scientific/editorial pass:** confirm bibliography and all source citations, full methods and thresholds, inferential vs descriptive FDR distinctions, external strain provenance/unknown independence, no misleading exposure-time or clinical phenotype claims, and journal formatting. The paper's novelty relative to the 2022–2026 source literature must remain appropriately limited.

A **separate, optional experimental dependency** exists only for **causal/functional** claims: condition-balanced independently randomized 10x/culture processing plus protein/ECM/outflow experiments. Absence of these experiments should restrict claims, not block a limited computational thesis.

**Non-blocking HOLD:** original investigators' 63↑/17↓ complete DEG80 identity table. The original authors will not be contacted as part of this workflow. No exact author80-versus-our83 identity overlap has been established.

**Storage/retention update (actual completed backup):** Original PNG/SVG main Figures, 16-gene S1, corrected Fig S1, superseded original Fig S1 (audit history only), 7d strain audit, library-batch model and cross-archive QA have now been uploaded as SIX verified ZIPs into [MASTER DEGREE / GSE313736 TM Thesis Figures Evidence](https://drive.google.com/drive/folders/1LjrbN9Lzr6DUXX0h8yaEXQNPLpHTm7Ye). Five original backups were read back through Drive metadata (original figure ZIP size 530,797 B), and the corrected supplemental ZIP was also uploaded successfully. The resulting Drive files are user-managed and do not prove public journal archiving. Use the **corrected** supplement, not the superseded original. The [image-level visual inspection and Drive archive ledger](FIGURE_VISUAL_REVIEW_AND_DURABLE_DRIVE_ARCHIVE_20261102.md) contains individual links and SHA256.


## 2026-11-02 completed actual figure visual check and manuscript link/source audit

The completed **image-level review** of main Fig1–3 and the corrected Fig S1 is documented in [visual inspection & Drive storage ledger](FIGURE_VISUAL_REVIEW_AND_DURABLE_DRIVE_ARCHIVE_20261102.md). The corrected [Fig S1 Actions #41435268815](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41435268815) **2/2 SUCCESS** visually eliminates the original supplemental sample-n label collision. The original result and corrected result **both have byte-identical gene-level statistics TSV files**. Corrected image is stored in [Drive](https://drive.google.com/file/d/1Zc_0p4l-HI7digyzQmaCIu2jWBKzT6tG/view?usp=drivesdk).

**New independent manuscript CI [#41436138351](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41436138351) 2/2 SUCCESS**, with evidence [artifact #13136980302](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41436138351/artifacts/13136980302) and [Drive backup](https://drive.google.com/file/d/1L2BvTyauXoW5d85XBek10IP7fTXB2YMT/view?usp=drivesdk). This checked **7 required publication/research docs, 15 local Markdown relative links (all exist), 90 external-link URL syntaxes**, the six required scientific source citations in the bibliography, exact 124 vs 90/194 and 19-column/rank18 caveats, corrected Figure S1 link and removal of author-contact dependency. **External websites were NOT fetched live** by this deterministic check; URLs passing syntax and bibliography inclusion is *not* an independent publisher proof.

### Actual remaining release item: ONE composed-manuscript proof

- **Completed:** four figure PNGs visually reviewed; Fig S1 originally FAILED and corrected re-render PASS; Github source QA; internal document link and citation-presence gate; original four-archive numeric cross-check; original and fixed figure plus model evidence backed up to Drive.
- **Remaining: 1 editorial package gate:** render the actual final thesis/journal document with selected publication format, verify all figure/table insertion positions and print-size caption typography, check reference-manager output/journal DOI details and human scientific prose. Do NOT claim submission-ready based only on source file CI. No new experimental work and **no original author outreach** required for the *restricted computational thesis*.
