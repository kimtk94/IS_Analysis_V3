# Metabolic Resilience — KoGES Public Training Prototype

**Date:** 2026-09-22  
**Purpose:** Validate the phenotype/pipeline logic before controlled-access individual-level KoGES analysis.  
**Data status:** **Educational/public training data only — not suitable for scientific inference or manuscript effect estimates.**

---

## 1. Data used

Five repeated public training files in Google Drive `KoGES_CKM_Thesis/01_Public_Data`:

1. `follow_01_base__FOLLOW_01_DATA.txt`
2. `follow_02_F1__FOLLOW_02_DATA.txt`
3. `follow_03_F2__FOLLOW_03_DATA.txt`
4. `follow_04_F3__FOLLOW_04_DATA.txt`
5. `follow_05_F4__FOLLOW_05_DATA.txt`

These files share `T00_ID`, so they can be linked longitudinally.

Important: the separate `base_data*` educational files use a different identifier family and must not be joined to these follow files by assuming ID compatibility.

---

## 2. Participant retention

| Wave | Rows | Complete 5-component MetS | Complete MBI |
|---|---:|---:|---:|
| BASE | 1,000 | 999 | 985 |
| F1 | 985 | 982 | 939 |
| F2 | 953 | 949 | 902 |
| F3 | 899 | 899 | 897 |
| F4 | 698 | 697 | 697 |

Unique IDs across the five files: **1,000**.

Number with complete MetS phenotype at at least:

- 1 visit: 1,000
- 2 visits: 1,000
- 3 visits: 999
- 4 visits: 995
- all 5 visits: 532

### Interpretation

The educational dataset confirms that the proposed analysis architecture can tolerate incomplete follow-up and is well suited to mixed models / time-to-event logic rather than complete-case-only analysis.

---

## 3. Prototype metabolic syndrome definition

For pipeline validation only, five components were constructed:

1. **Central obesity**
   - male waist ≥90 cm
   - female waist ≥85 cm

2. **Elevated BP**
   - SBP ≥130 mmHg OR DBP ≥85 mmHg
   - OR self-reported hypertension diagnosis (`HTN=2`)

3. **Elevated glucose**
   - fasting glucose ≥100 mg/dL
   - OR diabetes diagnosis (`DM=2`)

4. **Elevated TG**
   - TG ≥150 mg/dL

5. **Low HDL-C**
   - male <40 mg/dL
   - female <50 mg/dL

MetS prototype:

```text
component_count >= 3
```

The public training files do not include the full medication hierarchy planned for the controlled-access analysis, so the training phenotype is intentionally simplified.

Special missing codes converted to NA before calculations:

```text
55555  = follow-up nonparticipant
66666  = not surveyed
77777  = not applicable
99999  = unknown/nonresponse/unmeasured
```

---

## 4. MetS frequency in the training sample

| Wave | Evaluable N | MetS prevalence |
|---|---:|---:|
| BASE | 999 | 33.8% |
| F1 | 982 | 29.9% |
| F2 | 949 | 29.1% |
| F3 | 899 | 30.0% |
| F4 | 697 | 33.1% |

These values are **not epidemiologic estimates**. They only demonstrate that the endpoint generation code functions across repeated files.

---

## 5. Primary resilience transition prototype

### Recommended primary baseline definition

```text
baseline MetS component count <= 1
```

Training result:

- baseline healthy: **407 / 999 evaluable**
- incident MetS during F1–F4: **89 / 407 (21.9%)**
- stable ≤1 component with ≥3 valid visits: **179**

First incident MetS wave among the 89 events:

| First event | N |
|---|---:|
| F1 | 24 |
| F2 | 22 |
| F3 | 21 |
| F4 | 22 |

### Why `<=1` remains the preferred primary definition

Sensitivity comparison in the training sample:

| Baseline threshold | Baseline N | Incident MetS | Event rate* |
|---|---:|---:|---:|
| exactly 0 components | 163 | 19 | 11.7% |
| ≤1 component | 407 | 89 | 21.9% |
| ≤2 components | 661 | 224 | 33.9% |

\*Training-data event rate only, not a population estimate.

`<=1` provides a useful balance between a clinically favorable baseline state and adequate event count. `0 components` should be retained as a strict sensitivity definition.

---

## 6. Continuous Metabolic Burden Index prototype

The prototype used:

```text
MBI = mean[
  z(WAIST),
  z(SBP),
  z(FPG),
  z(log(TG)),
  -z(HDL)
]
```

Z-scores were calculated **within wave and sex** for this pipeline test.

Participants with ≥3 valid MBI observations and usable time values: **983 / 1,000**.

OLS individual-slope distribution:

- median slope/year: **0.0016**
- IQR: **-0.0294 to 0.0300**

Again, these numeric slopes are not scientific findings; the important result is that the repeated data support individual trajectory estimation.

### Controlled-access analysis upgrade

The real analysis should compare:

1. sex × wave standardized MBI
2. baseline-reference standardized MBI
3. raw-value mixed model with wave/site adjustment
4. clinical MetS component count

The conclusion should only be considered robust if candidate-protein findings are consistent across reasonable phenotype definitions.

---

## 7. MHO secondary analysis prototype

At training baseline:

- BMI ≥25 kg/m²: **438**
- baseline healthy (`components <=1`): **407**
- MHO (`BMI >=25` AND `components <=1`): **97**
- evaluable MHO with follow-up: **97**
- later incident MetS: **37 / 97**

### Decision

MHO→MUO is feasible as a **secondary/exploratory analysis**, but should not replace the broader metabolic-resilience primary outcome because subgroup restriction reduces sample size substantially.

---

## 8. Physical activity prototype

The educational follow files contain a repeated binary `Txx_EXER` field.

Valid exercise observations:

| Wave | Valid N | `EXER=2` N |
|---|---:|---:|
| BASE | 1,000 | 273 |
| F1 | 983 | 386 |
| F2 | 950 | 382 |
| F3 | 895 | 362 |
| F4 | 698 | 293 |

Among baseline-healthy training participants:

- baseline `EXER=1`: 69 incident events / 298 (23.2%)
- baseline `EXER=2`: 20 / 109 (18.3%)

**Do not interpret this as an exercise effect.** The educational data are being used only to prove that lifestyle interaction code can be constructed.

### Important difference from the 2024 integrated research codebook

The current integrated codebook indicates:

- `MET`: baseline only
- `EXERCUR`: waves 2–11
- `EXER`, `EXERFQ`, `EXERDU`: waves 3–11

Therefore the final controlled-access thesis should follow the integrated research codebook rather than the educational file layout.

---

## 9. Laboratory harmonization in the real KoGES analysis

The public training files expose simplified fields (`Txx_GLU0`, `Txx_HDL`, `Txx_TG`).

The actual integrated codebook provides:

- `GLU0_ORI`, `GLU0_TR`
- `HDL_ORI`, `HDL_TR`
- `TRIGLY_ORI`, `TRIGLY_TR`

and documents assay transitions:

```text
HITACHI 7600
→ ADVIA 1650
→ ADVIA 1800
```

The primary controlled-data candidate series is therefore:

```text
GLU0_HARM   = coalesce(GLU0_TR, GLU0_ORI)
HDL_HARM    = coalesce(HDL_TR, HDL_ORI)
TRIGLY_HARM = coalesce(TRIGLY_TR, TRIGLY_ORI)
```

subject to raw-distribution validation.

---

## 10. What this prototype proves

### Confirmed

1. Shared longitudinal IDs can be joined safely across the five follow files.
2. MetS component generation works.
3. `baseline healthy → incident MetS` works as an event phenotype.
4. MBI repeated-measures construction works.
5. Individual MBI slopes can be estimated for nearly all training participants with ≥3 waves.
6. MHO→MUO can be implemented as a secondary endpoint.
7. Time-varying/baseline exercise variables can be incorporated into the data model.
8. The analysis can be performed with a very small memory footprint by selecting only needed columns.

### Not yet confirmed

1. controlled-access KoGES actual missingness
2. actual baseline-healthy N
3. actual incident MetS event count
4. genotype × phenotype overlap
5. pQTL instrument availability in the genotyped longitudinal subset
6. final BP harmonization rule
7. real repeated physical-activity completeness
8. 10th follow-up availability in the eventual CODA workspace

---

## 11. Revised thesis hierarchy after this prototype

### Primary endpoint

```text
Baseline metabolic health (<=1 component)
→ incident metabolic syndrome
```

### Primary continuous secondary endpoint

```text
Repeated MBI trajectory
```

### Genetic validation model

```text
MBI_it ~ time + protein_GRS + protein_GRS×time + covariates + random effects
```

### Wellness layer

```text
protein_GRS × repeated exercise exposure
→ metabolic deterioration
```

### Exploratory

- strict baseline health = 0 components
- MHO→MUO
- trajectory classes
- HbA1c trajectory
- alcohol interaction

---

## 12. Stage 0 verdict

### Codebook feasibility

**GO**

### Public training prototype

**GO**

### Controlled-access participant-level feasibility

**PENDING**

### Overall

> **PROVISIONAL GO, with strong technical feasibility.**

The remaining gate is no longer whether the phenotype can be constructed. The main unresolved issue is whether the controlled-access KoGES sample provides adequate event counts and genotype overlap for the planned proteogenomic validation.

---

## 13. Next controlled-data outputs

When participant-level research data become available, generate in this order:

```text
01 KOGES_RAW_SCHEMA.tsv
02 KOGES_CORE_MISSINGNESS_BY_WAVE.tsv
03 KOGES_LAB_DISTRIBUTION_BY_WAVE.tsv
04 KOGES_METABOLIC_COMPONENT_COUNTS.tsv
05 KOGES_LONGITUDINAL_SAMPLE_FLOW.tsv
06 KOGES_MBI_TRAJECTORY_QC.tsv
07 KOGES_PA_AVAILABILITY_BY_WAVE.tsv
08 KOGES_GENOTYPE_INTERSECTION.tsv
09 STAGE0_GO_NO_GO.md
```

Do not begin proteome-wide candidate validation until items 1–8 establish the usable longitudinal genetic sample and event counts.
