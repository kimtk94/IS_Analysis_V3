# IS P3A — entire original GIGASTROKE EAS/EUR chr10 log(OR) source verification

**2026-11-23 KST.** Research only; 36 true QTL gene×tissue assays, one chr10 ischemic stroke locus. Formal coloc.abf / coloc.susie remains **0/36 HOLD**.

## Motivation

Earlier [PR #324](https://github.com/kimtk94/IS_Analysis_V3/pull/324) demonstrated original uploaded GWAS EAS+EUR odds_ratio = exp(beta) and original 95% CI = exp(beta ± 1.96·SE) in the first 750 source variants for each ancestry. These rows were not an entire original chr10 candidate region.

This workflow downloads original compressed publisher files **into disposable runner disk** and checks entire published MD5, previously audited SHA-256, full byte count and all original data row counts before scanning every chr10 record.

- EAS AIS GCST90104545 GRCh37: 260,689,208 compressed bytes; MD5 904e12e33ec3c2eeefe4dc93bec4995c; SHA-256 c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3, expected 6,725,216 rows.
- EUR AIS GCST90104540 GRCh37: 283,482,604 bytes; MD5 1bfe4ae8a5042fb24cb0a562b05b0d2f; SHA-256 e5d9809b608650e36a2808ca27bbb9de5b883ff9661d6d6056d02b3f792eacc0, expected 7,482,032 rows.

Standard gzip in these original files does not offer the documented BGZF/tabix random chromosome lookup. Full ephemeral source download is needed to obtain all chr10 rows; the file is never saved into any artifact.

## Statistical and source-safety gates

Each original chr10 row is checked for actual GWAS GRCh37 chromosome, position, signed effect allele, second allele, EAF, beta, SE, p-value, odds_ratio and confidence bounds. For valid numerical rows, compare log(original OR) with beta, and log(original 95% CI) with beta ±1.96SE. Report all source rows, nonfinite/missing/discordant counts and maximum absolute original log-space errors, not original SNP effect values.

The previous exact source 2Mb chr10 context in GRCh37 positions 104,163,676–106,220,351 is cross-checked against original EAS 3,913 and EUR 4,397 valid rows from earlier full-publisher source CI [#45020605713](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45020605713). This 2Mb region is NOT the entire QTL cis region after GRCh38→GRCh37 liftover.

Use original GWAS beta/SE to compute p-value as a diagnostic only; p-agreement alone cannot identify odds-ratio units. Source provenance pass is distinct from scientific OR-scale pass. Report discrepancies rather than silently dropping them.

## Scientific approval limits

Even if every chr10 original GWAS variant shows OR=exp(beta), original GWAS file-specific case fractions and SNP-level effective N are NOT proven. The exact publisher YAML lacks original sample_size/case_count/control_count fields; the separate original publication discovery cohort totals are not approved file-level N or case fraction s for coloc type=cc.

BBJ contributed to EAS AIS discovery, so a BBJ-only GWAS cannot count as independent replication. For QTL, original ROSMAP 15 gene groups yielded 89,844 deduplicated gene-wise variant records with AN/2=560, but this is allele count proxy, not direct measured normalized-expression sdY nor verified per-gene expression regression N.

**No invented N_eff, sdY=1, case fraction, GWAS allele flip or PP.H4.** All 36 scientific Coloc groups remain HOLD.

## Reproducibility and next analysis

- Script: scripts/is/audit_p3a_gigastroke_entire_chr10_logOR_scale.py
- Synthetic tests: tests/test_is_p3a_gigastroke_entire_chr10_logOR_scale.py
- GitHub Actions: .github/workflows/is-p3a-gigastroke-entire-chr10-logOR.yml
- Output artifacts: per ancestry original source hashes, aggregated entire chr10 beta–OR/CI consistency, study counts and 36 original HOLD decisions; never raw original GWAS SNP rows.

Next: establish source-justified case/control model s and appropriate study-based N, validate QTL model N/MAF/SE sdY estimate sensitivity, and construct formally harmonized original cis variants for 36 candidate QTL groups. Formal Coloc still requires independent approval.
