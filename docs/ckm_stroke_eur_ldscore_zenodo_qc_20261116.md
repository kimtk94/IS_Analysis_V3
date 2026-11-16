# CKM–Stroke: GRCh37 1000G EUR LDSC reference mirror integrity audit

Audit: 2026-11-16. **Technical/reference-material QC; no LDSC h2, rg, or causal inference performed.**

## Source pin

- LDSC's original precomputed EUR LD Score recommendation: [bulik/ldsc Heritability and Genetic Correlation wiki](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).
- Original Alkes tarball (`eur_w_ld_chr.tar.bz2`) was not reachable from the current machine at the legacy Broad URL on 2026-11-16.
- Accessible **third-party Zenodo mirror**: [Oliver Pain — European LD scores from 1000 Genomes, Zenodo 8182036](https://zenodo.org/records/8182036), DOI **10.5281/zenodo.8182036**, released 2023-07.
- DOI API original metadata: `https://zenodo.org/api/records/8182036`. Archived file: `eur_w_ld_chr.tar.gz`, **33,357,890 compressed bytes**, Zenodo-published MD5 **e2f16343c4cfaa76caa7d0c03d26b489**.
- Complete source gzip archive SHA-256 observed on Orca: **9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069**. The actual compressed bytes match the Zenodo published MD5 (not a separately authenticated checksum from the original Broad server).

## Actual archive components (streaming tar inspection; no extraction)

| Component | Coverage | Interpretation |
|---|---:|---|
| `{chr}.l2.ldscore.gz` | 22/22 | Full member gzip streams decompressed, gzip CRC verified, required `CHR SNP BP L2` columns present |
| `{chr}.l2.M_5_50` | 22/22 | Positive finite regression reference M values; these are the **default** M files for LDSC |
| `{chr}.l2.M` | 0/22 | Optional alternate M; `--not-M-5-50` must **not** be used with this mirror |
| `w_hm3.snplist` | 1 | A bundled reference list (~17,264,312 bytes); not identical byte count to the previously verified MRCIEU HM3 file (~17,170,716 bytes) |
| `6_old.l2.ldscore.gz`, `6_old.l2.M_5_50` | 2 | **Historical copies; do not include in chromosome 6 selector** |
| `README` | 1 | Additional bundled provenance |
| Archive member count | 49 | Included inside pinned full MD5 |

Script: `scripts/ckm_stroke/audit_eur_ldscore_zenodo.py`. It verifies Zenodo record DOI, fixed accession and checksum, source bytes, tar path safety, member count, LD score gzip headers and coverage, M numeric content. It writes only a small JSON evidence result to ignored `artifacts/`; it never extracts or commits the reference.

**Scope of technical PASS:** Zenodo mirror archive is intact **according to Zenodo's own MD5** and contains the standard 22 EUR LD Score and M_5_50 files. Does not prove bit-level identity to the original Alkes `.bz2` distribution, which may have been repacked to `.gz`. Does not establish that the bundled HM3 list is interchangeable with the previously pinned MRCIEU HM3 list. Compare marker IDs and release lineages explicitly.

**Remaining HOLD:** 9 duplicate GIANT BMI variants require effect-stat reconciliation; 44 zero P export values need a declared statistical treatment; use `BETA/SE` for signed statistics only after strand and scale QC, not a fabricated tiny P. Build and allele/LD reference harmonization, signed allele direction, UKB overlap with GIGASTROKE, LDSC h2/intercept/rg QC remain NOT_RUN. SBP/T2D authorization still blocks three-trait LDSC matrix.

## Executing without re-downloading unnecessarily

From the active Orca repository root:

```powershell
py -3 -m unittest discover -s tests -p test_ckm_stroke_eur_ldscore_zenodo.py -v
py -3 scripts/ckm_stroke/audit_eur_ldscore_zenodo.py --out artifacts/ckm_stroke/eur_ldscore_zenodo_20261116.json
```

Workflow `.github/workflows/ckm-stroke-eur-ld-reference-qc.yml` runs fixtures on changes, but the 33MB original archive streams **only** on deliberate `[eurld-full]` marker or manual workflow dispatch. No archive or original summary-stat files are uploaded to GitHub.
