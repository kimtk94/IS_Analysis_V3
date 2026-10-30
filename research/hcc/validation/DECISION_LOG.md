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

## 2026-10-30 — D-007: Evidence-gate hardening and Drive handoff

- Protocol governance PR: https://github.com/kimtk94/IS_Analysis_V3/pull/252
- Automated verification: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41057719470 (12 offline tests PASS)
- Google Drive reference: https://docs.google.com/document/d/1Ipk7kS6ukGYqwgwhSWHICs42EBqpu3eRFd82mBOkaAE/edit
- Rule update: `V12`/`V13` clinical event/count gates alone do **not** authorize a biology-guided clinical significance claim when the biological validation prerequisites are incomplete. Full biological prerequisite set includes V00, V01, V02, V03, V04, V05, V06, V07, V08, V09, V14 and V15.
- To mark any gate PASS, record explicit reviewer approval and reproducible analysis artifact path under `results/hcc/`, SHA256 of source (64 hex), exact code commit (40 hex), and checked date. V02/V04/V10 additionally require independent patient N; V12/V13 require verified event N.
- The static evidence checker validates *format and declared evidence*, **not** remote artifact contents, statistical effects, new discovery, or prospective clinical usefulness. A scientific reviewer still needs to inspect outputs. Even every workflow gate being PASS only makes a claim eligible for review.
- The user-facing handoff document resides with the original feasibility audit; the Markdown protocol and machine-readable status remain canonical in GitHub. Existing conclusions: molecular and clinical findings unproven.

## 2026-10-30 — D-008: PDC tumor-normal lineage and 11-plex analysis, RNA–Protein gate still BLOCKED

- Evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41059030994 and code on feasibility PR #251.
- Public PDC `paginatedCasesSamplesAliquots` on PDC000198 (Proteome) and PDC000199 (Phosphoproteome) yielded an identical set of **171 case, 331 sample, 331 aliquot identifiers** in both studies. Output saves only aggregates, not case/aliquot identifiers.
- Per study: **165 Primary Tumor samples, 165 Solid Tissue Normal samples, 1 Not Reported sample**. Distinct cases with **both tumor and normal: 160**; 5 tumor-only, 5 normal-only, 1 unclassified case. Publication's 159 paired analysis subjects should not be replaced with 160 without adjudicating inclusion/QC exclusions.
- Experimental metadata includes **33 TMT experiment numbers × 11 labels = 363 unique channel slots**, with 726 run metadata records; not 363 patients or 726 quantification observations. Header 667 (Protein) / 334 (Phosphosite) columns cannot be joined by direct full case/sample/aliquot string matching (0 matches).
- Scientific implication: this resolves which PDC cases are linked at the biospecimen layer, not the quantitative-matrix value layer, and **does not contain a BioSino RNA expression join**. This strengthens V01's source-level evidence but leaves **V02 BLOCKED**, all biological/clinical primary claims ineligible.
- Next required test: aliquot-to-(experiment, label) uniqueness and bridge/reference channels; replicate vs biological sample distinction; quantitative TSV header parsing; BioSino OEP000321 RNA gene-by-sample mapping and compatible patient keys.

## 2026-10-30 — D-009: Gao TMT analyte collision resolved in public metadata, V02 NOT promoted

- Source evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41062185538 ; official PDC `paginatedCasesSamplesAliquots` and HCC feasibility PR #251.
- In each PDC study, 726 `aliquot_run_metadata` rows span 33 experiment numbers. If keyed only by experiment+TMT label, 330 of 363 channel slots are ambiguously assigned to two aliquots. **The ambiguity is due to mixing two assay analytes:** `Proteome` and `Phosphoproteome`.
- Split by analyte: each class has **363 records, 363 experiment×label slots, 331 distinct aliquots, and 0 ambiguous slots**. `replicate_number` is blank in all 726 records; `fraction` alone leaves 60 ambiguous slots. The study field, analyte, experiment and label must not be conflated.
- Protein quantitative TSV has 330 Log Ratio + 330 Unshared Log Ratio columns (667 total including metadata); phosphosite has 330 Log Ratio columns (334 total including metadata). PDC values are **relative TMT ratios**, not absolute abundance. No public assay-specific TSV column→case mapping has yet been verified.
- NODE `OEP000321` was probed; HTTP 200 returned only a 906-byte public HTML shell, no source-accession text and no processed RNA matrix URL. Secondary publications report RSEM-quantified RNA from OEP000321, but actual file/permission/IDs were not verified.
- **Scientific gate implication:** PDC Proteome/Phosphoproteome experimental lineage is materially better resolved. Nonetheless cross-modal quantitative **V02=BLOCKED**, RQ1 inferential signal not testable, recurrence clinical claim ineligible. Next: match protein TSV normalized columns to `(analyte,experiment,TMT label)` mapping; locate and legally obtain processed OEP000321 RNA and verify same-patient overlap.
