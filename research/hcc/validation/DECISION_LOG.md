# HCC Research Decision Log (immutable decisions, append-only)

**Purpose:** preserve decisions after each analysis, especially why a checkpoint moved or a manuscript claim changed. Never rewrite older decisions; append a dated correction.

| UTC/KST date | Decision ID | Gate(s) | Evidence | Decision | Meaning for thesis |
|---|---|---|---|---|---|
| 2026-10-30 | D-001 | V00–V02 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41053421959 | Ng 51 triple-omic assay-flag biopsies = 49 distinct patients. Actual gene×sample matrix join is **not verified**. | Matched quantitative mechanistic claims are blocked; original paper's N must not substitute for a quantified join. |
| 2026-10-30 | D-002 | V04–V05 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41053421959 | Within 49 flagged patients: ALD-only 13, HCV-only 11, NAFLD-only 7, HBV-only 4; mixed etiologies separate. | Etiology effect modification is exploratory. Predeclare interaction contrast and assess BCLC, treatment, sex/age and platform confounding. |
| 2026-10-30 | D-003 | V01–V02 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41055576655 | PDC Protein Assembly proteome/phosphoproteome files and PRIDE SWATH/Progenesis candidates identified. | File presence is not proven RNA–protein join; avoid promoting P02 and V02. |
| 2026-10-30 | D-004 | V12–V13 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41052039696 | GSE14520 has 242 samples with raw Recurr status/time columns; unique patient rows and event code not audited. | No independent recurrence model/clinical significance claim; preserve RQ3 as conditional. |
| 2026-10-30 | D-005 | V06–V09 | https://doi.org/10.1038/s41586-019-0987-8 ; https://doi.org/10.1016/j.cell.2019.08.052 ; https://doi.org/10.1038/s41467-022-29960-8 | SOAT1/cholesterol, PYCR2/ADH1A, and HCC metabolic RNA/protein differences are **published prior findings**. | These are positive controls; novelty requires independent cross-etiology portability or biological source evidence. |
| 2026-10-30 | D-006 | All | research/hcc/validation/HCC_VALIDATION_MASTER.md | Separate method feasibility, effect significance, biological meaning, novelty and clinical utility. | Never claim high research significance because q<0.05 or CI is green. |

## Entry template

Date / ID / repo commit / checkpoint / dataset / source SHA256 / patient N / event N (if relevant) / outcome definition / analysis and plot artifact path / effect estimate / 95% CI / multiple-test correction / independent validation / evidence sufficiency / remaining bias / status change / exact allowed claim / next step.

## Automatic continuation rule

For a future request such as “다음”, first inspect the protocol, checkpoint registry and **latest CI**. Complete the first actionable public-data blocker. Save real quantitative evidence and update the ledger. If the hypothesis turns out null, negative, confounded, or restricted by access, log the limitation rather than manufacturing a PASS.

Study-governance PR is independent from the ongoing code/data-feasibility PR #251; merge or review both before treating them as available on main.
