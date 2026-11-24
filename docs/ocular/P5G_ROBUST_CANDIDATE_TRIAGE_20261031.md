# Ocular P5g: statistical robustness × prior-literature triage (2026-10-31)

## Evidence
This is **post-hoc screening of already-computed gene-level aggregates**, not a new raw-count differential expression model. Verified primary/LOO [P5b run 41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185), [artifact 12993824564](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185/artifacts/12993824564); independent-experimental-unit **within same study** [P5d run 41090827791](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791), [artifact 12995175346](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791/artifacts/12995175346). Reproducible ranking code in `scripts/ocular/ocular_p5g_candidate_screen.py` takes both downloaded artifact ZIPs; no donor-level counts required.

### Prespecified *for this exploratory ranking* filters
Among only the original P5 top **100** genes with full 9-pool signed summaries: primary donor-model BH-FDR<0.05 and abs(log2FC)>=1.0 and logCPM>=2; conservative 9-pool library-unit BH-FDR<0.05; all nine leave-one-pool-out refits retain BH-FDR<0.05 and abs(log2FC)>=0.8 with same sign; all nine independent donor pools have the same signed response; pool delta standard deviation divided by the absolute pooled mean <=0.35.

- **86/100** pass all these *exploratory* robustness criteria, of which **13 have already-identified published responses in the manually reviewed literature set** and **73 remain UNREVIEWED**, **not** 73 genuinely novel mechanisms. The literature-flag list is a partial/manual reference map and does not comprehensively audit the entire published transcriptome.
- LOO fit used `robust=FALSE` for runtime while primary donor model uses `robust=TRUE`; signed log2FC agreement is more comparable than exact FDR.
- These are re-use of same source material and cannot establish independent biological replication or exclude treatment-correlated library processing effects.

## Prioritized high-effect unreviewed candidates (not novelty GO)
| Candidate | P5 57-donor log2FC | P5d 9-pool-unit log2FC | Pool effect SD / abs mean | Decision |
|---|---:|---:|---:|---|
| HIF3A | +6.35084 | +6.41175 | 0.04944 | Full publication and pathway crosswalk pending |
| SLC16A12 | +5.23109 | +5.27563 | 0.08737 | Full publication and ocular biology crosswalk pending |
| CORIN | +5.11300 | +5.23145 | 0.09107 | Full publication and ocular biology crosswalk pending |
| CITED1 | −4.79664 | −4.92348 | 0.16812 | Full publication and cell-state validation pending |
| NEBL | +4.38820 | +4.19860 | 0.06581 | Full publication and cell-state validation pending |
| OTULINL | +4.26147 | +4.24665 | 0.11585 | Full publication and cell-state validation pending |
| AFAP1L1 | +3.66681 | +3.72804 | 0.04821 | Full publication and cell-state validation pending |

**Clinical caveat:** any of these being a strong DEX-associated in-vitro transcript is not a glucocorticoid-induced IOP risk gene and does not prove mechanism, causation or therapeutic value.

## Existing published findings moved from candidate to **known response**
- **SAA1, SAA2, OCA2**: already significantly induced in short- and long-term primary human TM DEX RNAseq by Kathirvel et al., *Scientific Reports* 2022, DOI [10.1038/s41598-022-12443-7](https://doi.org/10.1038/s41598-022-12443-7). These were among the largest original P5g effects; **they cannot be claimed as new expression effects**.
- **BDKRB1**: already studied by independent DEX-TM RNA and protein experiments in *International Journal of Ophthalmology* 2022, [PMID 36017046](https://pubmed.ncbi.nlm.nih.gov/36017046/).
- **HIPK2, PTK2B, FKBP5, FAM107A, ZBTB16, ATP10A, EDIL3, MAOA, MIEF1, DNM1L, MYOC, ANGPTL7, TSC22D3** are manually flagged based on prior clinical/DEX/TM studies (see `P5B_LITERATURE_GENE_CROSSWALK_20261031.md`). The presence of genes in the known dictionary need not imply every gene passes the current top100 filter.

## Novelty and inferential GO/HOLD/STOP
1. **GO** stable in-vitro high-effect gene-level candidate shortlisting for subsequent evaluation.
2. **HOLD** systematic *exact* DEG overlap with full 2022/2026 supplementary tables, true independent experimental replication, cell-state/score threshold robustness and treatment-vs-library technical identification. Note 2026 Mehrotra et al., *IOVS*, 100nM DEX/2d study [PMID 41636427](https://pubmed.ncbi.nlm.nih.gov/41636427/) has 857 TM and 2,086 SCE significant genes and full Supplementary Table 4 available in publication; need exact overlap instead of named positive controls.
3. **HOLD** original GSE313736 paper's exact high-effect gene overlap and 74→70→59 donor attrition reasons.
4. **STOP** claims that any unreviewed DEG is a new causal glaucoma gene, clinically predictive adverse event or proven drug target without independent validation.

## Next implementation
Acquire the official 2026 IOVS supplemental TM complete DEG table (no guesswork in URLs), harmonize gene identifiers, evaluate all-significant-gene overlap, sign concordance and 857-TM-DEG replication fraction against the same gene universe. Then recompute the priority list **before** considering mechanistic experiments.
