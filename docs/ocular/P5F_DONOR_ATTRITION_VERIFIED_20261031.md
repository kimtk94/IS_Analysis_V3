# P5f: verified GSE313736 original-vs-final sample attrition (2026-10-31)

**Direct original-source validation:** GitHub Actions corrected [run 41091777826](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091777826), [evidence artifact 12994467026](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091777826/artifacts/12994467026). Source `scripts/ocular/ocular_p5f_attrition.py` reads the GSE313736 official family SOFT and downloadable cell metadata; comparisons use exact GSM accession mapping, not inferential pool name parsing.

## Verified chain and limitations
| Stage | Donors | Evidence |
|---|---:|---|
| GEO study abstract, all participants | **74** | GEO abstract, original 316,799 singlet cells |
| Across **18 official GSM sample descriptions** | **70** | 9 source pools, each Vehicle and DEX; 70 donors listed under both treatment arms |
| Final public `GSE313736_HTMC_scRNA_cell_metadata.tsv.gz` | **59** | All 59 appear in original GSM descriptions and are paired across treatments, 169,398 SNG rows |
| Strict BeamA score>=.75 and >=30 cells per arm | **57** | P4 actual barcode joins, 18 raw 10x matrices |

Described sample donor IDs absent from final processed metadata: **11**. Study abstract participants absent from GSM descriptions: **4**. Thus 74 to 59 is a **15-donor numerical difference**, broken into 4 abstract-to-GSM and 11 GSM-to-final-metadata. These comparisons do **not** explain why individuals were excluded; demographic attrition mechanism remains unverified.

The paper reports 316,799 singlets vs final downloaded 169,398 metadata rows, a numerical difference of **147,401 cells**. These counts may represent different QC stages or processing subsets; we cannot conclude all 147,401 were filtered by a particular cause.

## Pool-level sample description to final donor metadata
| Pool | GSM description donors | Final demultiplexed donors | Absent from final |
|---|---:|---:|---:|
| A | 8 | 8 | 0 |
| B | 8 | 7 | 1 |
| C | 8 | 7 | 1 |
| D | 8 | 7 | 1 |
| E | 8 | 7 | 1 |
| F | 8 | 6 | 2 |
| G | 8 | 7 | 1 |
| H | 7 | 6 | 1 |
| I | 7 | 4 | 3 |
| **Total** | **70** | **59** | **11** |

The sample-description donor sets contain all the final metadata donors in their corresponding pools, verified per pool. No individual donor ID, age, genotype, disease, PRS or clinical phenotype is saved in artifacts.

## Research implications
- Primary expression model uses 57 high-confidence BeamA donors, **not all 74 originally recruited**.
- Because exclusions may relate to capture efficiency, quality and cell states, missingness/selection bias remains possible. Pool I loses 3/7 described donors (42.9%), Pool F 2/8, more than most pools; report across-pool differences and sensitivity.
- **GO** numerical donor source traceability and 18/18 library-to-pool mapping. **HOLD** provenance of exclusion reasons and the original 316,799-to-169,398 singlet reduction.
- Clinical adverse-event effect remains STOP: these are cultured TM cells without confirmed clinical steroid-induced IOP outcomes.
