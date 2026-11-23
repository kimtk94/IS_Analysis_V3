# IS P3A — original GIGASTROKE EAS/EUR beta, SE, p, effect scale and case fraction

**2026-11-23 KST**. Ischemic stroke only, 36 genuine QTL gene×tissue groups, 20 genes, one chr10 locus. All formal Coloc **0/36 HOLD**.

## Prior primary evidence, already verified

- Original GIGASTROKE 2022 publisher Supplementary Table 1 independently verified ([Actions #45035576990](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45035576990)): EAS AIS 19,032 cases / 237,242 controls; EUR AIS 62,100 cases / 1,234,808 controls. **These are paper cohort totals and not certified SNP-specific or file-level regression N**.
- BBJ is an EAS AIS discovery contributor (17,651 cases, 159,610 controls); do not label an ordinary BBJ stroke GWAS as independently replicated against GIGASTROKE EAS.
- Previously full publisher-source integrity checked: EAS GCST90104545: 260,689,208 bytes; MD5 904e12e33ec3c2eeefe4dc93bec4995c; SHA256 c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3. EUR GCST90104540: 283,482,604 bytes; MD5 1bfe4ae8a5042fb24cb0a562b05b0d2f; SHA256 e5d9809b608650e36a2808ca27bbb9de5b883ff9661d6d6056d02b3f792eacc0.
- Existing [SNP matching audit](IS_P3A_36_SNP_ORIENTATION_GWAS_AUDIT_20261122.md): 89,597 EAS and 128,282 EUR allele-signed SNP×gene-assay matches, but records are not independent SNPs or fully approved coloc inputs.
- Official GRCh37 source YAML metadata (both accessions) were previously confirmed **without sample_size, case_count, control_count**, in [CI #45033915115](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45033915115).

## New primary-publisher evidence

The [Mishra et al. original 2022 Nature publication](https://www.nature.com/articles/s41586-022-05165-3) states individual GWAS sources used logistic or Cox regression (latter for some longitudinal cohorts), mentions log[OR] QC filtering, and documents meta-analysis. This supports a **log-odds methodological context**, not an incontrovertible certified scale for the beta column of these exact two files.

The CI retrieves primary original European PMC full-text XML, both exact publisher source YAMLs and only the initial 256KiB HTTP 206 compressed byte range of each original GRCh37 GWAS file. The full original GWAS was already MD5/SHA256 verified in prior CI; it is **not redownloaded in this step**.

Audit points:
1. Assert article DOI 10.1038/s41586-022-05165-3 and logistic regression / Cox regression / log[OR] terms in original primary publication.
2. Publisher YAML: exact original accession, trait Ischemic stroke, GRCh37, original published compressed file name and pinned MD5; classify missing case counts and effect-unit labels explicitly.
3. Original GWAS compressed real rows: GRCh37 chromosome and positions, effect and other allele, EAF, numeric beta, SE and source p-value. Check effect allele orientation and reject malformed or nonfinite values.
4. Compare two-sided normal P derived from beta/SE with source P using per-ancestry median, p95 and discordant-row aggregate counts. **Agreement identifies mathematical consistency, not the numerical unit of beta**.
5. Publish only aggregate evidence, source hashes and 36 scientific HOLD group decisions; no variant IDs, beta rows, sample IDs or genotypes.
6. Keep original publisher cohort totals as context only, no silently approved model case fraction s, effective sample size or log-OR beta transformation. BBJ as independent replication remains prohibited.

## Research acceptance

Scientific status remains HOLD until the precise original GCST90104545 / GCST90104540 source effect-size scale and case fraction mapping to the original summary files is justified, alongside the model-N / sdY approximation for ROSMAP QTL. An exact signed-allele and genome-build harmonization plus appropriate binary phenotype coloc (type cc, case fraction s) is then required.

- Script: scripts/is/audit_p3a_gigastroke_original_beta_se_p_scale.py
- Tests: tests/test_is_p3a_gigastroke_original_beta_se_p_scale.py
- Workflow: .github/workflows/is-p3a-gigastroke-beta-se-p-scale.yml
- Output: artifacts/is/IS_P3A_36_GIGASTROKE_original_beta_scale_QC_HOLD.json

**CI technical success does not equal scientific Coloc eligibility.**

## 2026-11-23 new direct original-file OR-scale source result

The original uploaded EAS GCST90104545 and EUR GCST90104540 GRCh37 compressed files both contain **populated** original odds_ratio, ci_lower and ci_upper columns **in addition** to beta, standard_error and p_value. The source audit of 750 real original SNPs **per ancestry** found:

| Original source | Real rows beta, SE, p | Original OR populated | exp(beta) = original OR | 95% CI = exp(beta ± 1.96 SE) | Maximum absolute |ln(OR)-beta| |
| --- | ---: | ---: | ---: | ---: | ---: |
| EAS AIS GCST90104545 | 750 | 750 | 750/750 | 750/750 | 4.9730e-15 |
| EUR AIS GCST90104540 | 750 | 750 | 750/750 | 750/750 | 4.9730e-15 |

No OR/CI source rows were missing in the 1,500 checked rows; p-values also agree with z=beta/SE to median absolute delta log10(p) 0.00090 EAS and 0.00124 EUR. This is **direct source evidence that the uploaded column beta is the logarithm of the uploaded odds_ratio in these tested rows**, not merely a statistical assumption. The updated source runner requires the exact source OR/CI agreement and rejects a drift in original data.

**Scientific boundary:** This 750-row per-ancestry prefix is **not** the full chr10 candidate region, the complete files, or a proof of original per-variant sample N. Original YAML remains missing case_control N; EAS discovery BBJ overlap remains present. The next approval step must audit all candidate chr10 original source rows or prove that the original file's OR column is globally generated from beta under the exact file format; then justify case fraction and quantitative QTL model inputs. Formal coloc remains 0/36 HOLD.

Evidence: original branch CI #45161635002, artifact is-p3a-gigastroke-EAS-EUR-beta-SE-P-original-source-HOLD. Original paper XML context from prior verified CI SHA256 is labeled as previously verified, not reread live in this audit.
