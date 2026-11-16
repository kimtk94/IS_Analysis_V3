# GLGC EUR source-SNP export readiness and scientific boundaries

Status: **code prepared; export regression CI queued; zero live full-chromosome source-SNP exports attempted by this new pipeline.** Evidence cut: 2026-11-16 KST.

## Why the earlier technical PASS was not an analysis dataset

The successful chr22/chr21/chr20 GLGC indexed-source QC jobs wrote a compact 7-column scratch file (`CHR,POS,EA,NEA,BETA,SE,P`) and deleted it at job exit, preserving only JSON QC results. That is appropriate for bounded source tests, but insufficient for a future MR/LDSC workflow: per-variant `N`, effect allele frequency and source `rsID` were unavailable and the normalized rows were not retained.

Do **not** treat `ckm-glgc-batch-...json` as GWAS data, do **not** calculate LDSC or MR from row counts, and do **not** assert original full-file SHA256 integrity from the indexed shard checks.

## New optional export

`scripts/ckm_stroke/glgc_chr_shard.py` now accepts:

```bash
python -m scripts.ckm_stroke.glgc_chr_shard \
  --trait LDL_C --chromosome 22 --max-rows 5000000 \
  --export-source-snp-dir artifacts/source_snps \
  --out artifacts/glgc_chr22_ldl_qc.json
```

The optional normalized `*.tsv.gz` source-SNP export has:

| Output | Source interpretation |
| --- | --- |
| `CHR`, `POS` | Chromosome and `POS_b37`, GRCh37 |
| `EA`, `NEA` | GLGC ALT effect allele and REF non-effect allele, not externally harmonized |
| `BETA`, `SE`, `P` | Reported GLGC summary association columns |
| `RSID` | Raw `rsID` field; may be `NA` and must not be fabricated |
| `N` | Source per-variant `N`; missing/nonpositive values block optional export |
| `EAF` | Source `POOLED_ALT_AF`; not ancestry-adjusted to a new outcome |
| `NEG_LOG10P` | Original `pvalue_neg_log10`; preserves a higher dynamic range for later review |

This file is emitted only when the **entire selected chromosome** obtains technical QC PASS, has zero duplicate SNP keys and the new export-required `N`/`NEG_LOG10P` fields are finite and valid. JSON contains the exported SHA256 and byte length. The output is **UNHARMONIZED SOURCE SNPs**, so `analysis_ready_data_persisted=false` and `scientific_status=HOLD` remain mandatory.

Default use **without** `--export-source-snp-dir` still writes JSON only and avoids changing the prior large-data CI cost.

## Controlled GitHub Actions

- Workflow: [`ckm-stroke-glgc-source-snp-export.yml`](https://github.com/kimtk94/IS_Analysis_V3/blob/research/ckm-stroke-staged-actions-20261111/.github/workflows/ckm-stroke-glgc-source-snp-export.yml)
- Manual-only inputs: `trait`, `chromosome`, `source_qc_run_id`.
- Downloads the *earlier independent real-source QC artifact*, checks same trait, chromosome, source and technical PASS before downloading.
- Re-reads the official indexed BGZF URL over TLS, rechecks sample size, variant statistics, and compares source row/exclusion counts with the earlier QC artifact. The output gzip SHA256 is checked before upload.
- On success only, retains the public per-chromosome source-SNP gzip and JSON provenance in a GitHub Actions artifact for **7 days**. HOLD/failure uploads diagnostic JSON only (if available), not an unchecked data table.
- This workflow is not scheduled, and does **not** silently enqueue a full-source job on every push. GitHub may require the workflow to be present on the default branch to dispatch from the Actions UI; do not merge Draft PR #284 without approval.
- The separate lightweight regression: [`#43860195691`](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43860195691), **QUEUED as of this note**. Never infer tests PASS from a queued run.

## Essential scientific checks still absent

1. Full original GLGC genome-wide file byte-level integrity / checksum and stable source-version check.
2. Raw-source and outcome reference-build reconciliation, strict ALT/REF matching and strand complement, palindromic AF-aware filtering, variant-specific `N` QC.
3. Official regression-statistical scale (LDL/HDL inverse normalized and logTG transformation), corresponding outcome interpretation, European ancestry-matched LD reference panels, sample overlap checks.
4. LDSC, cross-trait enrichment, locus discovery, cis-pQTL MR, colocalization, and EAS replication.
5. `P=0` underflow variants must be separately evaluated from the official negative-log10 P data before MR or genome-wide significance thresholding; the existing SNP-only QC explicitly counts and excludes those rows. A zero count in tested chromosomes does not establish zero across all autosomes.

**Current observed evidence**: chr22/chr21/chr20 × LDL-C/HDL-C/logTG (9/66 source jobs) technical PASS, 6,457,067 retained SNP rows across phenotype/chromosome slices. chr19 queued; scientific completion remains 0/7, P1 HOLD.

Related: [Draft PR #284](https://github.com/kimtk94/IS_Analysis_V3/pull/284), [Issue #285](https://github.com/kimtk94/IS_Analysis_V3/issues/285).

## SHA reproducibility and extremely small P-value safety

The optional gzip writer explicitly fixes gzip `mtime=0` and omits the filename in the gzip header, so identical source rows yield a stable **export gzip SHA-256** on reruns. This still does not establish official upstream source file integrity. If any otherwise eligible source SNP has `pvalue == 0` due to floating-point underflow, the optional export is held until original `pvalue_neg_log10` can be interpreted correctly; truncating the strongest associations would bias downstream MR and LDSC. Technical SNP-only QC may still PASS while the **export gate is HOLD**.
