# TPMI 433.21 targeted SNP and GIGASTROKE gate — 2026-11-21

## Verified evidence (do not promote to scientific GO)

- Original 58,302,393-byte PLINK2 file, SHA256 `ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62`, 1,739,048 GWAS records.
- Previous original-source run [#44872892525](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44872892525): offline-negative-tests and full-original-genomewide-source jobs completed successfully; original aggregate artifact ID 14694881855 still unexpired when inspected.
- Independent OBS_CT audit run [#44873247155](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873247155): artifact ID 14695131212 still unexpired when inspected.
- One source variant meets P <= 5e-8, but its identity and GRCh37 mapping are **not yet verified in an observed run**.
- PLINK2 OBS_CT = 204099 at every recorded variant, while phenotype page lists 9249 cases + 304660 controls = 313909. The case/control breakdown among the 204099 regression observations is UNKNOWN; no binary N_eff is approved.

## Current executable pipeline

Workflow: `.github/workflows/ckm-tpmi-43321-lead-grch37.yml` in branch `research/ckm-stroke-staged-actions-20261111` (PR #284, not merged).

1. PR and push: test source one-hit parsing, bidirectional assembly mapping, and GIGASTROKE lookup request offline. This is **code contract validation**, not primary-source result validation.
2. Explicit original-source execution (workflow_dispatch `run_original_source=true`, or a push whose commit message includes `[tpmi-lead-original]`): verify the two previous source-run pins, download original TPMI 58MB **once**, SHA-check, scan all 1,739,048 rows, save ONLY the unique genomewide-significant variant.
3. Using validated single SNP: query official Ensembl forward and reverse assembly maps and original/mapped reference bases; reject ambiguous mapping, invalid position, unequal REF, and unsupported indels.
4. Publish GRCh37 lookup **target only**; GIGASTROKE effect is not downloaded or inferred and scientific result remains HOLD.

## Proof needed to move to gate 3 PASS

- A new GitHub Actions run ID with SUCCESS for `extract-pinned-original-lead`, `bidirectional-build-and-allele-map` and `prepare-gigastroke-grch37-target`.
- Archived outputs `tpmi_43321_GRCh38_single_lead_public_GWAS_QC.json`, `tpmi_43321_GRCh37_Ensembl_bidirectional_map_HOLD.json` and `tpmi_43321_GRCh37_GIGASTROKE_target_only_HOLD.json`.
- Record exact chr/pos/REF/ALT/A1, source P and A1 logOR, Ensembl roundtrip, mapped REF and target position; do not infer replication from a coordinate alone.

## Proof needed for gate 4 scientific consideration

- Independently sourced GIGASTROKE *original* target-row chr/pos/alleles, effect allele, beta or logOR, SE and P; handle flips, complements, and strand-palindromic SNPs with frequency evidence.
- Document each cohort's ascertainment and participant overlap. EUR and EAS effects must not be called independent simply because ancestry labels differ.
- Obtain actual PLINK2 433.21 per-analysis cases and controls; do not equate `OBS_CT` to effective sample size or use 2025 PheWeb/2026 study sample counts interchangeably.
- No LDSC, MR, coloc or independent replication claim until each corresponding analysis-specific gate is met.

## Current disposition

- Technical full-file QC: **PASS** from prior pinned run.
- New one-SNP extraction, Ensembl liftover, GIGASTROKE matched effect: **UNVERIFIED, HOLD**.
- Scientific independent replication and effective-N: **HOLD**.
- New workflow run: **not confirmed** from available GitHub commit workflow-run/status lookups (which do not enumerate all events). Do not interpret an empty lookup as proof of zero runs.

## Additional targeted AIS source matching implementation

Code: `scripts/ckm_stroke/match_tpmi_43321_gigastroke_ais.py`, with offline cases `tests/test_ckm_tpmi_43321_gigastroke_ais.py`.

The main TPMI workflow now includes an **optional** manual input `run_gigastroke_ais=true`, guarded by `run_original_source=true`. It downloads the official 283,482,604-byte GCST90104540 GRCh37 gzip **once**, verifies publisher MD5 `1bfe4ae8a5042fb24cb0a562b05b0d2f` after scanning all rows, enforces exact mapped coordinates and uniqueness, and extracts only one target association and raw effect orientation.

- On no target locus, ambiguous allele pair, duplicates, wrong published hash, or invalid numeric input, the job returns HOLD or fails closed.
- On a valid source coordinate plus non-palindromic allele match, the single-locus **technical** gate can PASS; sample independence, valid binary N_eff, independent replication and LDSC remain HOLD.
- No successful run, target rsID, numeric original AIS beta or allele direction has been observed from this new code at the time this section was committed.
- A draft/feature-branch workflow may not be available through GitHub's manual Run workflow button until registered on the default branch. This is a deployment constraint, not proof the code passed its tests.

## Verified executed TPMI source and GRCh37 mapping — 2026-11-21

**Source-backed GitHub Actions run:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44878942971

Individual GitHub job API conclusions:
- `offline-allele-assembly-tests`: **SUCCESS** (job 136769756519).
- `extract-pinned-original-lead`: **SUCCESS** (job 136769807292), artifact 14697201895, full pinned original SHA verified, 1,739,048 records, one source P <= 5e-8.
- `bidirectional-build-and-allele-map`: **SUCCESS** (job 136770356497), artifact 14697022641, Ensembl forward/reverse and both assembly REF bases checked.
- `prepare-gigastroke-grch37-target`: **SUCCESS** (job 136770422991), artifact 14696224277.
- `source-verified-gigastroke-ais-one-locus`: **SKIPPED** (not yet harmonized with GIGASTROKE original).
- `report-ci-status`: **FAILED** (job 136770469771) with `gh: Resource not accessible by integration (HTTP 403)` when attempting a PR #284 comment, unrelated to original GWAS parsing or Ensembl mapping. The overall run is therefore failure even though the four scientific-data preprocessing jobs passed. The main workflow was corrected in follow-up commit `288c109a4cde7ec9772e658d8a72b75c68fd5dd4` to write only to `GITHUB_STEP_SUMMARY`.

Original pinned lead extracted from real official TPMI PLINK2 GWAS:
- `rs187759247`, source `chr19_15462568_T_C`; GRCh38 `19:15462568:T:C`; A1 effect allele C.
- Source A1 frequency `0.0212101`, OR `1.68467`, log(OR) `0.5215696989267228`, log(OR) SE `0.0662893`, Z `7.86808`, P `3.60115e-15`, OBS_CT `204099` (not effective N).
- Successfully mapped to **GRCh37 `19:15573379:T:C`**, effect A1 C, reverse liftover and source/mapped reference alleles verified.
- No SNP effect from GIGASTROKE AIS GCST90104540 has yet been checked at the locus, and independence/actual case/control N remain unverified. **Scientific HOLD** for replication, LDSC, MR and coloc.
- Updated narrow milestone: original source integrity and coordinate mapping TECHNICAL PASS; GIGASTROKE comparison/outcome causal interpretation HOLD.
