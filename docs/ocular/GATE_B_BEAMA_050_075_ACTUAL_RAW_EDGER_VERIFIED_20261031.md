# Gate B — actual BeamA 0.50 vs 0.75 original RAW UMI edgeR re-fit: verified (2026-10-31)

## Evidence and execution provenance

- **Final successful GitHub Actions**: [run #41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259), conclusion **SUCCESS**, downloadable public [artifact #13041944372](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259/artifacts/13041944372). Artifact actually downloaded, decompressed and checked (summary JSON, P9 16-gene TSV, 83/88 union, 18-source QC). Earlier failed raw-count runs (#41198827259, #41198957708, #41199257286) are preserved in Actions. They failed on invalid/truncated GEO gzip transfers **before** completion. Resilient official retry + full gzip footer/CRC on ALL small and large source files resolved it.
- **Original comparator**: P5b [run #41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185), score>=0.75, 57 donor paired raw-count edgeR.
- **Re-analysis**: score>=0.50, same 57 donor pairs; use original 18 public 10x matrices streamed and donor UMI summed, edgeR \`~ pair_id+treatment\` with TMM, robust quasi-likelihood, same minimum 30 cells per condition and same software **edgeR 4.0.16 / R 4.3.3**.
- All 18 original library count matrices gzip CRC **PASS** and observed expected feature order. Runner-only 114-column pseudonymous donor-level counts and pairing metadata were removed and were **not uploaded**. Only aggregate gene-level FDR/effect estimates, anonymous MDS/MA figures, QC and manuscript summary are published. \`score0.50\` actual joined **48,782 vehicle cells + 55,048 DEX cells = 103,830 cells**, vs metadata-only 48,903 and 55,652 (some metadata barcodes absent from filtered 10x matrices). All 57 paired donors verified **after** joining actual barcodes. **Pool G DEX barcode recall 0.92222** (passes earlier >0.90 gate but retains a QC warning; targeted per-pool follow-up appropriate).
- Actual Gate B source code: \`scripts/ocular/ocular_gate_b_050_p4_retry.py\`, \`ocular_gate_b_050_edger.R\`, \`ocular_gate_b_050_vs_075_compare.py\`; reusable manual \`.github/workflows/ocular-gate-b-050-raw-edger.yml\`.

## Numerical reproducibility of threshold change

| Metric | BeamA score>=0.75 | BeamA score>=0.50 | Interpretable conclusion |
|---|---:|---:|---|
| Eligible independent donor pairs | 57 | 57 | Same donor identities |
| edgeR FDR-testable genes | **16,050** | **16,728** | +678 new filtered-in genes, so FDR universes differ |
| FDR<0.05 | **13,717** | **14,304** | Both unusually widespread DE; tech/batch confounding remains |
| Common Ensembl genes | \multicolumn{2}{c}{**16,050**} | Full baseline set present in 0.50 |
| Common-gene log2FC Pearson correlation | \multicolumn{2}{c}{**0.99911597**} | Internal threshold stability, not independent replication |
| Median shared |Δlog2FC| | \multicolumn{2}{c}{**0.0112462**} | Strong effect estimate stability |
| Same sign among 16,050 shared | \multicolumn{2}{c}{**15,851/16,050 (98.76%)**} | Includes near-zero effect genes, so rate depends on sign threshold |
| FDR<0.05 in both | \multicolumn{2}{c}{**13,563**} | Same donor / source cohort: non-independent |
| P9 candidates, same sign/FDR in both | \multicolumn{2}{c}{**16/16**} | Includes condition-dependent across other cohorts |
| 3-study stable candidates | \multicolumn{2}{c}{**3/3 MYOC/TGFBI/MMP2**} | **GO cutoff stability** |

Note: the markdown above presents cross-column metrics as combined values; neither side is a separate external biological replicate. *Statistical significance does not eliminate DEX/10x-library confounding.*

### RNA candidate effects verified from artifact

| Gene | 0.75 log2FC | 0.50 log2FC | Both adjusted FDR<0.05 | Overall evidence tier |
|---|---:|---:|---|---|
| **MYOC** | +4.025688 | +3.907811 | yes (0.50 FDR 8.78e-41) | 3-cohort signed candidate, steroid response reference |
| **TGFBI** | -1.869208 | -1.881365 | yes (5.35e-34) | 3-cohort signed candidate; RNA cannot establish ECM deposition |
| **MMP2** | -0.589067 | -0.595032 | yes (1.82e-17) | 3-cohort signed candidate; zymography needed for activity |
| **ACTA2** | +0.469496 | +0.419432 | yes (3.12e-5) | **Condition-dependent across chronic 7-day TM**, not universal core |
| **FKBP5** | +3.984015 | +3.962229 | yes (1.43e-46) | GR positive reference in GSE313736, external P9 value not verified |
| **SDC1** | -0.553539 | -0.546771 | yes (2.61e-9) | Signed discordance across studies persists |

All P9 16 candidate genes retain original sign across the two **within-cohort** threshold choices.

## High-effect genes: 83 vs 88 is NOT original-author 80 overlap

Fixed internal definition: \`BH-FDR<0.05 && abs(edgeR log2FC)>=3\`.

- **Score0.75**: 83.
- **Score0.50**: 88.
- **Shared by stable Ensembl ID**: **81**.
- **0.75 only**: **ACP5** (−3.0393 → −2.9568) and **HLF** (+3.0405 → +2.9604), both changed just across the absolute 3.0 threshold while FDR remained strongly significant.
- **0.50 only**: **SYT7**, **PCK1**, **ALDH1L1**, **KRT19**, **MARCH10**, **SMIM5**, **LINC01423**. Of these, **PCK1/ALDH1L1/MARCH10/SMIM5/LINC01423 were absent from the score0.75 _tested_ universe**, so their inclusion is not a gene-level effect reversal. SYT7 and KRT19 moved across the log2FC=3 threshold.
- This is **our own score0.75 vs score0.50 high-effect set comparison**, NOT independent source author high-effect list validation.
- GEO abstract says 63-up/17-down (total 80) from the original GSE313736 publication deposit; **individual 80 gene symbols/Ensembl IDs were not found in official GEO Series supplemental list**. Thus original author80-vs-our83 identity mapping remains **HOLD**. Do not claim the 83/88 overlap answers this distinct question.

## Scientific implications

**PASS/GO**: full count-matrix processing, matched 57 donor reaggregation, score0.50/0.75 \`gene-wise log2FC\` numerical concordance, all 16 previously frozen P9 genes and three-core stable sign/FDR, effective selection threshold robustness. This materially resolves one **classification-confidence sensitivity** question.

**HOLD**: ability to separate steroid pharmacological effects from condition-specific 10x library chemistry, sample processing or post-treatment BeamA state selection. A very large fraction of all genes are statistically altered in the same source cohort (85%+ FDR), which may be partly technical. Other cell-state annotations (\`scpred_vanzyl2022\` TM_Fibro or JCT) may have distinct biological composition and donor coverage; these must first be verified before claiming broader cell-state robustness.

**HOLD**: original 80-gene author-level source list, exact experiment phenotype-linked IOP, independent protein/functional validation.

**STOP**: treating cutoff stability or two conditions as randomized donor-independent 100-vs-500 nM dose comparison, clinical glaucoma classification or microbial mediation.

## Next focused gate

Perform **low-cost alternative-annotation metadata eligibility** to determine whether \`scpred_vanzyl2022/TM_Fibro\`, \`BeamA ∩ TM_Fibro\` and \`JCT\` have meaningful ≥30 cells per arm across independent donors. Only order additional **large** 18-matrix raw-count processing if sufficiently powered for three-core candidates; do not re-run 6GB GEO downloads automatically.

*Technical reproducibility note: all release artifacts are publicly downloadable via the GitHub Actions archive subject to retention, and source commit/actual execution IDs are fixed above. The source's public-access status does not grant permission to export the separate individual-level human genotype VCF.*
