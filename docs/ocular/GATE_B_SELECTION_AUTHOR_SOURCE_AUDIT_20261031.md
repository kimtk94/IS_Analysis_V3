# Ocular Gate B — cell-state selection vs original 80 DEG provenance (2026-10-31)

## Audited public evidence and actual runs

- **Metadata/official GEO source twin audits:** [GitHub Actions #41198598417](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198598417) **2/2 SUCCESS**, artifacts [eligibility #13041527282](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198598417/artifacts/13041527282) and [author-list #13040504794](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198598417/artifacts/13040504794). Both ZIP contents downloaded and inspected.
- Independent metadata implementation [#41198576725](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198576725) also **SUCCESS**; artifact #13041721835 confirms source thresholds and count totals, but is **not another biological cohort**.
- Public GEO source GSE313736: 169,398 cell metadata rows, **59 paired donors**; file SHA256 `41e4a6bf5792fe58e023b92f0129119dcf2b49877f62dd580b9dab5de48f5192`. Official Series SOFT SHA256 `f44a2ec6888e6a7e7ec8416f90d97427f71a8106e0dac9e348373dcce7d9a850`. No donor IDs, barcodes, private metadata, or donor-level expression exported.

## Gate B1 — post-treatment annotation and donor eligibility **PASS metadata only**

| Annotation | Minimum cells in each paired treatment | Eligible biological donors | Selected metadata Vehicle cells | Selected metadata DEX cells |
|---|---:|---:|---:|---:|
| BeamA score≥0.50 | 30 | **57** | **48,903** | **55,652** |
| BeamA score≥0.75 | 30 | **57** | **40,115** | **46,454** |
| BeamA score≥0.50 | 50 | 57 | 48,903 | 55,652 |
| BeamA score≥0.75 | 50 | 57 | 40,115 | 46,454 |
| JCT score≥0.50 | 30 | **30** | 1,872 | 2,312 |
| JCT score≥0.75 | 30 | **8** | 481 | 564 |
| scpred TM_Fibro (unscored) | 30 | **57** | 13,233 | 19,259 |
| scpred TM_Fibro (unscored) | 50 | **51** | 12,831 | 18,610 |

- Lowering BeamA confidence from 0.75 to 0.50 at ≥30 cells/arm **does not add donors**; it adds **8,788 Vehicle and 9,198 DEX metadata cells** among the same 57 donors.
- Only **55/57** strict BeamA donor identities overlap with donors selected under the alternative scpred TM_Fibro≥30 criterion; the alternative has 57 donors but a **different donor set and much smaller total cell numbers**. It must not be narrated as an interchangeable anatomical class.
- For the strict BeamA selection, within-donor median DEX minus Vehicle *selected-cell fraction* is **−0.42146 percentage points**, vs BeamA score≥0.50 **−0.63107 pp**. Pooled cell-count fractions and paired-donor fraction deltas differ; neither establishes a causal change in cell identity.
- All of these counts are **metadata-based only**, not 18-matrix barcode-recalled raw UMI joins. The successful metadata audit **did not rerun edgeR**. Marker/DEG signed stability must be evaluated from the original integer count matrices at the alternative threshold; donor count equality alone is not sufficient.

## Gate B2 — original paper 63 UP / 17 DOWN exact gene list: audit PASS; biological overlap HOLD

The official GSE313736 Series SOFT source contains **five** listed Series supplementary files: genotype VCF, scRNA cell metadata, classification/demux, donor metadata, and `GSE313736_RAW.tar`; it does **not** provide an explicitly published signed 80-DEG gene-level table through that file inventory. Do **NOT** download any individual donor genotype VCF for this audit. Full publication supplementary material or author-provided differential tables may still exist elsewhere.

| Comparison | Official study abstract | Our 57-donor score≥0.75 paired edgeR |
|---|---:|---:|
| High-effect increased | **63** reported | **67** FDR<0.05 & log2FC≥3 |
| High-effect decreased | **17** reported | **16** FDR<0.05 & log2FC≤−3 |
| Total | **80** reported | **83** |

**Critical**: 83 minus 80 equals 3 **in counts**, **not** 3 specific discordant genes. Different donor selection, QC, model, and potentially threshold definitions can produce far more than 3 non-overlapping gene identities. The GEO abstract's downregulated threshold notation appears ambiguous ("<3", not necessarily "<−3"). The author-side exact list has **not** been acquired or matched; exact gene overlap is **HOLD**.

Artifact `GateB_our_83_high_effect_genes_NOT_author_overlap.tsv` contains our 83 actual Ensembl IDs, symbols, effect/FDR and an explicit unverified-author-overlap flag. It is a reproducible input for later true published-list crosswalk, **not** a reconstructed author list.

## Gate B3 — true raw UMI cutoff effect check: separate higher-cost job

An initial [.50 raw-score run #41198827259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198827259) failed **at GEO transfer before any count analysis** with `BAD_GZIP_MAGIC`; no differential-analysis results were produced. An official GEO retry wrapper was added to recover legitimate files without relaxing source gzip/feature/barcode/NNZ integrity. A corrected [run #41198957708](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198957708) is executing as of the first writing of this report; **its scientific outcome remains unverified until its artifact is reviewed**.

Prespecified interpretation gates:
- Retain true 57 paired donors by actual 10x feature/barcode joins (do not assume metadata eligibility implies successful count joins).
- Identical robust edgeR NB-QL design to fixed score≥0.75 primary: donor as biological unit, not cell; compare common filtered gene universe and explicitly report changed FDR as well as logFC.
- Check MYOC, TGFBI, MMP2 (three cross-cohort stable transcript candidates), ACTA2 and other condition-dependent candidates (EGR3/LPAR1/TFPI/SDC1), all original 83 high-effect genes.
- Stop/hold if 10x source cannot be fetched faithfully or if cutoff-dependent effect directions flip. Do **not** fabricate true count refit from metadata counts or CPM.
- Treatment-specific 10x library processing and post-treatment cell annotation are unidentifiable causal confounders even if high concordance survives confidence-cutoff sensitivity.

## New priorities

**Gate A:** 3 published TM experiments have already been compared; MYOC↑, TGFBI↓, MMP2↓ remain signed consistent across 48h/7d and differing platforms, ACTA2 is context-dependent. This is independent **descriptive expression comparison**, not clinical, protein or causal DEX-dose proof.

**Gate B:** finish actual 0.50 raw-count result, quantify shared gene/score robustness and report author 80 exact-list **HOLD** until verifiable gene IDs are found. JCT 8 donors at 0.75 is too sparse for an analogous full same-power comparison.

**Gate C:** genomic/disease relevance or experiment only after Gate B limitations are explicit, and only with independently accessible locus/tissue information. Ocular surface/oral 16S remains a separate HOLD project.

**No student/user manual action is required for these public metadata/analysis stages.** Clinical phenotype genotype access or wet-lab experiments require distinct permissions and resources.
