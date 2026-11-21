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
