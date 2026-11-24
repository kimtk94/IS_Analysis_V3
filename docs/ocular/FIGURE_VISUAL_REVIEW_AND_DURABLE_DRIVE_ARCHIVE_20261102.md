# Ocular TM thesis — actual visual PNG/SVG review and verified Google Drive archiving (2026-11-02)

**Figure audit stage:** VISUAL GO for the four inspected PNG images after the fixed supplementary build. **Overall scientific publication readiness:** limited computational thesis GO; journal submission requires complete final proof/bibliography/editing, irrespective of CI success. This document is an image-level visual review, **not** just a PNG CRC/SVG parse.

## 1. Actual image inspection findings (visualized at native original resolution)

| Asset | Original GitHub CI source | Actual image resolution | Visual outcome | Final recommended action |
|---|---|---|---|---|
| Fig 1: MYOC/TGFBI/MMP2/ACTA2 three-study sign matrix | [run #41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011) | **2,598 × 1,471 PNG**, SVG also stored | **GO**: glyphs, colored rows, source labels and bottom caption visible; no observed collision or clipping | Preserve. Final caption must explain 7d CPM+1 ratio ≠ 48h edgeR/DESeq2 and 3 strains ≠ 7 donors |
| Fig 2: same donor 57 confidence and matched-55 scPred | [run #41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011) | **2,956 × 1,424 PNG**, SVG | **GO**: both side-by-side panels have separated labels, legible log2FC and footnote, no observed overlap | Preserve. Clarify same cohort, original 10x libraries reused |
| Fig 3: Hallmark EMT 194 gene sign totals | [run #41208972011](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41208972011) | **2,517 × 1,125 PNG**, SVG | **GO with caption update**: bars and source denominator values readable, no observed clipping | The 124/194 bar is the three-study **aggregate** sign criterion, *not* strict 3-strain agreement. Supplement states **90/194** |
| Fig S1: baseline effect-size/abundance sign concordance | Original [#41389639386](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41389639386) | **2,775 × 1,081 PNG** | **ORIGINAL VISUAL FAIL**: `n=` point labels collide near y≈100% in both series and overlap in expression quartiles. CRC/source SHA alone previously did not catch this | Do not submit original asset |
| **Corrected Fig S1**: same statistics, count table below the plot | [fixed CI #41435268815](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41435268815) | **2,957 × 1,333 PNG** | **VISUAL GO**: sample-count tables placed independently below each panel; no observed count/footnote collision after second corrective layout change | **Use fixed archive #13135314010 exclusively** |

The corrected Fig S1 was directly re-opened visually. It separates legend, axes, bin labels, source sample counts and caveat text. This image-level review has a narrow scope: real appearance in these archived PNGs; journal-specific physical print sizing, color accessibility and the assembled final thesis PDF remain separate editorial checks.

### Exact source-integrity and no-statistics-change proof

- The two published aggregate tables in the original Fig S1 artifact and corrected Fig S1 artifact were compared **byte-for-byte**: `ocular_effect_stratum_sign_and_fdr.tsv` **SAME** and `ocular_sensitivity_panel_effects.tsv` **SAME**. Thus the corrective rendering did **not** change test universes, gene effects or evidence thresholds.
- Original main figure PNG SHA256: Fig1 `42b3a734234db6b777143ddc79f93bc52055c06feac99ef7b0f6eec5793a7582`, Fig2 `d5fd890969960354413a37acc0f2c52207fbeece6258f63af003c954320c0545`, Fig3 `61778f075a5bc288a6c77586e2c756bec20046cb8d684ce8b87928c42860d5`.
- Final corrected Fig S1 PNG SHA256: `e4d135fe6c663b3aaf423d049cee2451604bb0038c41af8727098696671cf478`.
- GitHub Actions job #41435019882 rendered a first version that resolved `n=` labels but its bottom caption was still too close to the new count table in local visual preview; `8a881e0be636236a95ad11f9cfd3a8f56b53e139` moved the footer below, and **fixed #41435268815 2/2 CI jobs SUCCESS**. The actual fixed PNG was visually inspected again.
- Original Fig1–3 source ZIP SHA256: `60d560b7f3704f81e809834de72654874b88a1e140088b9eba774bbcbb0a4a91`.
- Corrected Fig S1 ZIP SHA256: `c2e24c8cd6f6adb47f1222deb5281e4ff01381f863ad64aebd72008e380bd2a9`.

## 2. Completed cross-connector durable Google Drive backup (readback verified)

**Destination:** [MASTER DEGREE / GSE313736_TM_Thesis_Figures_Evidence_20261102](https://drive.google.com/drive/folders/1LjrbN9Lzr6DUXX0h8yaEXQNPLpHTm7Ye). All archives were downloaded from the actual source GitHub Actions via official artifact references and uploaded through the Google Drive connector; the **five original archive items** were subsequently verified by Drive metadata and reported exact stored byte sizes. The corrected replacement was also successfully uploaded.

| Material | Drive file (completed write) | Source |
|---|---|---|
| **Main Figures 1–3, 3 PNG, 3 SVG, 16-candidate S1 and legends** | [530,797-byte original ZIP](https://drive.google.com/file/d/1tUB9sljGAViCfHCcQMlhSYTG8Ui_l5JG/view?usp=drivesdk) | GitHub Actions #41208972011, artifact #13045363678 |
| **Corrected Fig S1** (use for submission, not original) | [final visually corrected supplementary ZIP](https://drive.google.com/file/d/1Zc_0p4l-HI7digyzQmaCIu2jWBKzT6tG/view?usp=drivesdk) | GitHub Actions #41435268815, artifact #13135314010 |
| Superseded original Fig S1 (audit history only) | [original overlapping caption ZIP](https://drive.google.com/file/d/1pOHlZ2nj8a271hQup1GZb9ZS9L8imj-4/view?usp=drivesdk) | GitHub Actions #41389639386, artifact #13116138830 |
| Three-strain GSE301517 evidence (16 candidates + 194 EMT) | [archived TSV/JSON](https://drive.google.com/file/d/1L0_p16yZTj7laf83kHZobDfEAxu2T1bb/view?usp=drivesdk) | GitHub Actions #41423487666 |
| Original GSE313736 design rank and BeamA pool decomposition | [source evidence ZIP](https://drive.google.com/file/d/12G9rJGvW2leetgxkImDLPOT6Gf8wUQd3/view?usp=drivesdk) | GitHub Actions #41430965326 |
| Four-source and original figure CRC provenance | [review manifest ZIP](https://drive.google.com/file/d/1OELBOyfWRUkfPYEUdfCARfKgZ_7RZYXv/view?usp=drivesdk) | GitHub Actions #41432752451 |

Original GitHub Actions artifact retention is finite. **Real Drive backups now exist**; nevertheless they are user-managed files, not an unchangeable public archive. This is not an assertion that the Google Drive files have been permanently shared or independently published.

## 3. Manuscript captions and submission gates

**Main Fig 3 required addendum:** "Concordance is assessed on per-study aggregate signed estimates. In an additional 3-strain GSE301517 test, 90 of 194 EMT genes retained the signed effect across all three study aggregates *and* the three individual seven-day TM strains; no independent gene-wise seven-day FDR is inferred."

**Corrected Fig S1 caption:** "In the same 57 paired donors, BeamA classification score thresholds 0.75 versus 0.50 showed genome-wide sign robustness. A matched-55 donor comparison of BeamA and scPred TM_Fibro displayed stronger sign changes in low baseline effect genes: 78.258% agreement for baseline |log2FC| <0.25, n=5,349. Counts below each panel represent common tested Ensembl genes in baseline-defined bins. These sensitivity comparisons reuse original donors and condition-specific 10x libraries and should not be interpreted as independent biological replication or a batch-corrected DEX effect."

**Visual gate:** PASS for the four actual fixed/original figure PNGs at the inspected screen resolution. **Machine-source and statistical gate:** original source 4-archive QA [#41432752451](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41432752451) PASS; after fixed Figure S1, its two gene-level tables are unchanged. **Full final submission:** still HOLD until bibliography/style/source table integration and final formatted thesis assembly.

**Research communication:** original authors will not be contacted. The original DEG80 gene table remains an explicitly non-blocking source limitation; clinical steroid-induced IOP causal and predictive claims remain STOP.
