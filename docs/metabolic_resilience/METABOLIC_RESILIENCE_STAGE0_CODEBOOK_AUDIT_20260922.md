# Metabolic Resilience — KoGES Stage 0 Codebook Audit

**Date:** 2026-09-22  
**Source audited:** `KoGES 지역사회기반코호트 반복 추적조사 통합자료 코드북(ver5.0)_202412.xlsx`  
**Status:** **PROVISIONAL GO — participant-level raw-data audit still required**

---

## 1. Executive verdict

The longitudinal KoGES codebook is strongly compatible with the proposed **Metabolic Resilience** thesis.

The five core domains required for the proposed Metabolic Burden Index (MBI) are represented repeatedly from baseline through the 10th follow-up:

- Waist circumference
- Blood pressure
- Fasting glucose
- Triglycerides
- HDL cholesterol

In addition:

- HbA1c is represented across all 11 survey waves.
- Hypertension, diabetes, and dyslipidemia diagnosis/treatment variables are longitudinal.
- Height and weight are repeated, so BMI can be derived even though no exact integrated `BMI` variable was identified.
- Alcohol status is repeated.
- Physical activity is usable, but **the original plan must be modified** because `MET` is baseline-only.

The study should therefore proceed to participant-level feasibility testing.

---

## 2. Wave mapping

The integrated codebook uses `차수 = 1–11`.

| Codebook wave | Study wave | Approximate period |
|---:|---|---|
| 1 | Baseline | 2001–2002 |
| 2 | 1st follow-up | 2003–2004 |
| 3 | 2nd follow-up | 2005–2006 |
| 4 | 3rd follow-up | 2007–2008 |
| 5 | 4th follow-up | 2009–2010 |
| 6 | 5th follow-up | 2011–2012 |
| 7 | 6th follow-up | 2013–2014 |
| 8 | 7th follow-up | 2015–2016 |
| 9 | 8th follow-up | 2017–2018 |
| 10 | 9th follow-up | 2019–2020 |
| 11 | 10th follow-up | 2021–2022 |

`A01_*` therefore means baseline, while `A11_*` means the 10th follow-up.

---

## 3. Primary metabolic phenotype feasibility

### 3.1 Waist circumference

Canonical integrated variable:

`WAIST`

Wave-specific variables:

`A01_WAIST ... A11_WAIST`

Description:

> 허리둘레 평균

**Decision:** retain as a core MBI and MetS variable.

---

### 3.2 Blood pressure

Available repeatedly:

- `SBP_L`
- `DBP_L`
- `SBP_R`
- `DBP_R`

for all 11 waves.

There is no separate integrated variable named simply `SBP` or `DBP` in the audited codebook.

### Recommended derivation

Participant-level raw data should be used to compare:

1. left-arm values,
2. right-arm values,
3. mean of both arms when both are available,
4. single available arm when only one is valid.

A prespecified harmonized rule should then be fixed before outcome modeling.

**Do not select the rule solely from association strength.**

---

## 4. Critical laboratory harmonization finding

A major issue for 20-year trajectory analysis is laboratory instrument change.

The codebook documents:

- through 2002-08: **HITACHI 7600**
- 2002-09 to 2010-12: **ADVIA 1650**
- from 2011: **ADVIA 1800**

For fasting glucose, HDL-C, and triglycerides, KoGES provides original and official transformed values.

### Fasting glucose

- `GLU0_ORI`: original measurement
- `GLU0_TR`: transformed value

Official conversion equations documented in the codebook include:

For pre-2002-09 measurements:

`Glu_ADVIA1650 = 1.757 + 0.907 × Glu_HITACHI7600`

For measurements from 2011 onward:

`Glu_ADVIA1650 = -0.9 + 1.000 × Glu_ADVIA1800`

### HDL-C

Pre-2002-09:

`HDL_ADVIA1650 = 3.265 + 0.806 × HDL_HITACHI7600`

From 2011 onward:

`HDL_ADVIA1650 = -1.1 + 0.932 × HDL_ADVIA1800`

### Triglycerides

Pre-2002-09:

`TG_ADVIA1650 = 21.073 + 0.942 × TG_HITACHI7600`

From 2011 onward:

`TG_ADVIA1650 = 1.7 + 0.988 × TG_ADVIA1800`

### Recommended harmonized series

Subject to raw-data verification:

```text
GLU0_HARM   = coalesce(GLU0_TR,   GLU0_ORI)
HDL_HARM    = coalesce(HDL_TR,    HDL_ORI)
TRIGLY_HARM = coalesce(TRIGLY_TR, TRIGLY_ORI)
```

Why this is attractive:

- waves that were already measured using the target ADVIA 1650 scale do not require a transformed variable;
- early/late instrument periods have official transformation variables;
- it minimizes artificial longitudinal jumps due to assay platform changes.

### Required raw-data checks

Before final adoption:

- missingness of ORI and TR by wave
- whether TR is populated for all participants requiring conversion
- distribution continuity around 2010/2011
- median/IQR and mean/SD by wave
- site-specific distributions
- implausible values and special missing codes

The final manuscript should report this laboratory harmonization explicitly.

---

## 5. Recommended primary MBI v0.2

The original score remains suitable, but use harmonized laboratory values:

\[
MBI_{it}=
mean[
z(WAIST_{it}),
z(SBP_{it}),
z(GLU0\_HARM_{it}),
z(log(TRIGLY\_HARM_{it})),
-z(HDL\_HARM_{it})
]
\]

### Standardization recommendation

Primary candidate:

- sex-specific z-scores
- standardized within wave

Reason:

Repeated biomarker distributions can shift because of aging, laboratory changes, and secular trends.

Sensitivity analyses:

1. baseline-reference z-score
2. pooled z-score with wave fixed effects
3. raw clinical thresholds / MetS count

Using multiple phenotype definitions reduces dependence on an arbitrary scale.

---

## 6. Primary clinical transition outcome

Baseline metabolically healthy population:

```text
MetS component count <= 1
```

Primary event:

```text
first transition to >=3 MetS components
```

Component construction:

1. Waist circumference
2. Elevated BP or antihypertensive treatment
3. Elevated fasting glucose or diabetes treatment
4. TG >=150 mg/dL
5. Low HDL-C

Medication-support variables confirmed in the codebook include:

- `TREATD1`: current hypertension treatment
- `TREATD2`: current diabetes treatment
- `TREATD9`: current dyslipidemia treatment
- `DRUGHT`: antihypertensive drug use
- `DRUGDM`: oral diabetes drug use
- `DRUGLP`: lipid-lowering drug use

Exact medication hierarchy should be determined after participant-level completeness is known.

---

## 7. BMI

No exact integrated variable named `BMI` was identified in the audited codebook.

However both are repeated at every wave:

- `HEIGHT`
- `WEIGHT`

Therefore derive:

\[
BMI = rac{WEIGHT(kg)}{[HEIGHT(m)]^2}
\]

Use BMI as:

- adiposity covariate
- MHO/MUO sub-study variable
- sensitivity phenotype

Do **not** automatically include BMI inside the primary MBI because waist is already a central adiposity component and including both may overweight adiposity.

---

## 8. Physical activity — required redesign

### Important finding

`MET` is only available at baseline in the integrated codebook.

Therefore the previous idea of using repeated `MET-min/week` as the main lifestyle interaction cannot be implemented across 20 years from this integrated variable.

### Repeated alternatives

#### `EXERCUR`
Current exercise status.

Available from wave 2 through wave 11.

#### `EXER`
Regular exercise sufficient to produce sweating.

Available from wave 3 through wave 11.

#### `EXERFQ`
Exercise frequency:

- 1 = 1–2 times/week
- 2 = 3–4 times/week
- 3 = 5–6 times/week
- 4 = almost daily

Available waves 3–11.

#### `EXERDU`
Average duration per exercise session, minutes.

Available waves 3–11.

### Recommended strategy

Primary longitudinal physical-activity modifier:

```text
EXER + EXERFQ + EXERDU
```

Construct a repeated exercise-dose proxy, for example:

```text
weekly_exercise_minutes_approx =
frequency_midpoint(EXERFQ) × EXERDU
```

This is not formally identical to accelerometer-derived MVPA or MET-min/week, so terminology must remain conservative.

Baseline `MET` can be used as:

- baseline physical-activity covariate
- sensitivity stratification
- validation against the exercise proxy where overlap is possible

---

## 9. Alcohol

Longitudinal variables are strong enough for covariate use:

- `DRINK`: drinking status, waves 1–11
- `TOTALC`: total alcohol intake, g/day, waves 1–11 in the integrated codebook

Raw missingness still needs confirmation.

---

## 10. Sleep

No dedicated repeated sleep-duration / insomnia / bedtime variable was identified in the integrated codebook using:

- 수면
- 취침
- 기상
- 불면
- 잠
- sleep

The only “수면” occurrence located was explanatory text stating that resting-time activity excludes sleep.

### Decision

Sleep should **not** be a planned core modifier for the initial thesis based on the current integrated codebook.

If an auxiliary/non-integrated sleep table is later identified, sleep can be restored as exploratory.

---

## 11. What the Drive files actually contain

A folder contains:

```text
BASE_2001_2002.bin
F01_2003_2004.bin
...
F10_2021_2022.bin
```

However inspection shows these files are stored as XLSX and their sizes/sheets match the official wave-specific **codebooks**, not participant-level cohort data.

Example `F10_2021_2022.bin` contains codebook-style sheets rather than ~5,500 participant observations.

### Consequence

Current Drive material is sufficient for:

- variable mapping
- assay harmonization design
- wave availability
- analysis specification

but **not** sufficient to calculate:

- actual missingness
- eligible participant N
- incident MetS event count
- longitudinal retention after phenotype restrictions
- actual distribution of laboratory values

These require participant-level KoGES files.

---

## 12. Provisional GO/NO-GO assessment

### GO

**Core longitudinal phenotype:** GO  
All five core metabolic domains are represented longitudinally.

**Clinical MetS transition:** GO  
Diagnosis/treatment support variables exist.

**Continuous MBI:** GO  
Required repeated continuous measures exist.

**MHO→MUO:** GO  
BMI can be derived and metabolic components are present.

**Alcohol adjustment:** GO.

### MODIFY

**Physical activity interaction:** MODIFY  
Use repeated exercise variables instead of repeated MET.

### HOLD

**Sleep interaction:** HOLD  
No suitable dedicated repeated variable identified yet.

### PENDING RAW DATA

- actual analysis N
- missingness
- event rate
- genotype intersection
- laboratory distribution continuity
- exact BP harmonization rule
- exercise-dose completeness

Overall verdict:

> **PROVISIONAL GO**

The main thesis concept survives the codebook audit and is more feasible than initially expected.

---

## 13. Stage 0.5 raw-data audit requirements

The next pass should produce:

```text
KOGES_FILE_INVENTORY.tsv
KOGES_RAW_SCHEMA.tsv
KOGES_VARIABLE_AVAILABILITY_BY_WAVE.tsv
KOGES_CORE_MISSINGNESS_BY_WAVE.tsv
KOGES_LONGITUDINAL_SAMPLE_FLOW.tsv
KOGES_LAB_DISTRIBUTION_BY_WAVE.tsv
KOGES_BP_HARMONIZATION_AUDIT.tsv
KOGES_PA_AVAILABILITY_BY_WAVE.tsv
STAGE0_GO_NO_GO.md
```

### Minimum sample-flow counts

```text
N baseline
N with genotype
N with >=2 valid metabolic visits
N with >=3 valid metabolic visits
N with >=5 valid metabolic visits
N baseline healthy
N incident MetS
N stable healthy
N incident T2D
```

---

## 14. Raw-data phenotype rules to test first

### Valid missing codes

The codebook defines special values:

```text
55555 = did not participate in follow-up
66666 = not surveyed
77777 = not applicable
99999 = unknown / nonresponse / unmeasured
```

These must be converted to missing before numeric calculations.

### Initial biological QC ranges

These are audit flags, not automatic exclusion thresholds:

```text
HEIGHT: 120–210 cm
WEIGHT: 30–200 kg
WAIST: 40–180 cm
SBP: 60–260 mmHg
DBP: 30–160 mmHg
FPG: 40–600 mg/dL
TG: 20–2000 mg/dL
HDL: 10–150 mg/dL
HbA1c: 3–20 %
```

Outliers should be inspected before exclusion.

---

## 15. Analysis hierarchy after raw-data audit

### Primary

```text
Baseline healthy
→ incident MetS
~ candidate protein genetic proxy
```

### Secondary

```text
Repeated MBI
~ time + protein_GRS + protein_GRS × time + random effects
```

### Lifestyle

```text
Repeated exercise proxy
× protein_GRS
→ metabolic deterioration
```

### Exploratory

```text
MHO → MUO
trajectory classes
HbA1c-based glycemic trajectory
alcohol interaction
```

---

## 16. Immediate next step

Run the lightweight input-audit script supplied with this report against:

```text
/srv/is-analysis/data/koges
```

It intentionally inventories files and schemas before any full-table loading.

This is important because actual participant-level data appear not to be present in the Drive codebook archive.

Once the raw table names and columns are known, the next script can calculate actual wave-level missingness and participant eligibility without guessing file structure.
