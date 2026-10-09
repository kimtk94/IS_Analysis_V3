# IS broad discovery — next acquisition and LD gate (2026-10-09)

## Verified GIGASTROKE ancestry from original GWAS Catalog YAML
**30 studies = 5 phenotype files × 6 study families**:
- MULTI_ANCESTRY GCST90104534–38
- EUR GCST90104539–43
- EAS GCST90104544–48 (canonical present)
- AFR GCST90104549–53
- HIS GCST90104554–58
- SAS GCST90104559–63

The older metadata registry had 15 UNKNOWN ancestry records; source YAML resolves those records. Five multi-ancestry studies were previously misclassified as AFR/AMR because only the first samples block was read. Source sample sizes, per-population breakdown, GRCh37 assembly and MD5 are available in `GIGASTROKE_SOURCE_VERIFIED_V3.tsv`. Published source: Mishra et al., Nature 2022, doi:10.1038/s41586-022-05165-3.

All population-specific or combined analyses may share contributing participants; they must not be treated as independent replication cohorts.

## Actual genotype-based EAS LD audit
Tool: `scripts/is/build_anchor_ld_proxies.py`

This server's 2023 PLINK2 `--r2-unphased` command reports "under development", so the script uses PLINK2 `--export Av` and computes pairwise-complete (r^2) from variant-major sample dosage for 504 1000G EAS individuals. It handles missing data and zero-variance variants and reports allele-orientation status. For 500-kb windows around the selected lead variants:

| BBJ legacy locus | Expanded group | r²≥0.2 | r²≥0.8 | Orientation |
| --- | --- | ---: | ---: | --- |
| L001 | G0007 | 20 | 4 | Reverse REF/ALT |
| L002 | G0015 | 174 | 14 | Exact |
| L003 | G0022 | 30 | 10 | Exact |
| L004 | G0023 | 26 | 1 | Exact |

Reference-only LD **does not establish** statistical independence, 95% credible sets or variant causality. The remaining 26 regions are not genotype ready.

## Acquisition: explicit execution only
Status: 25 GIGASTROKE accessions beyond the 5 EAS canonical datasets are not available as local normalized EAS inputs. Priority next: EUR AIS `GCST90104540` for broader discovery, followed by EUR AS/CES/LAS/SVS and cross-ancestry meta-analysis with cohort overlap caveats.

Planning the first download (no network side effects):

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
python3 scripts/is/verify_gigastroke_ancestry_v3.py
python3 scripts/is/acquire_gigastroke_verified.py --accession GCST90104540
```

Optional manual download with resumable curl and integrity enforcement:

```bash
python3 scripts/is/acquire_gigastroke_verified.py --accession GCST90104540 --execute
```

**No automatic giant transfers**, no overwriting canonical GWAS, and no results presented as finished before file size and MD5 checks.

## Reproducibility
```bash
python3 scripts/is/run_broad_phase2_evidence.py
python3 -m unittest discover -s tests -p 'test_is_*' -v
```

## Next science gate
1. Normalize EUR + trans-ancestry GWAS to per-ancestry validated GRCh37 variant schema with REF/ALT/effect-allele and EAF checks.
2. Do LD clumping using ancestry-matched reference genome-wide, not geographic proximity.
3. Fine-map credible sets across all independent signals with SuSiE; restore missing or distant gene regulatory targets.
4. Extend QTL and human stroke cell-atlas evidence to new loci before equal-coverage ranking.
