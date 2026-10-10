# IS native BBJ allele provenance, existing GTEx ALT-major coverage and nine-gene context

**Executed:** 2026-10-11 KST.  
**Scope:** existing on-server BBJ ischemic-stroke original GWAS, saved eQTL Catalogue GTEx-v8 QTD interval files, previous 646 recovered gene–tissue SNP-level coloc inputs, and six-donor adult healthy brain reference data. **No new restricted, paid, CKB or TPMI GWAS data accessed.**

## Summary: science gate changes

| Check | Before | New directly tested outcome |
|---|---|---|
| Original BBJ GWAS effect allele beta/EAF vs archived GTEx-coloc GWAS fields | Prior 8,309 SNPs only | **8,956/8,956 PASS, zero allele/beta/EAF or chain conflicts** |
| Archived 646 SNP-level coloc numerical reproducibility | 646/646 direct ABF H0–H4 reproduced | **Unchanged: reproducibility, not causal validity** |
| eQTL Catalogue ALT-major absence in 646 cached inputs | 0 ALT-major AC/AN entries across 1,224,224 SNP–assay rows | **Direct raw GTEx QTD source shows ALT-major variants ARE present upstream** |
| Ten raw-QTD ALT-major SNPs missing from coloc at BBJ L003 | Reason unknown | **Nine absent from original BBJ native GWAS; one reported with low BBJ MAF and still absent from archived coloc** |
| Healthy adult control brain reference × nine candidate genes | Earlier descriptive gene-level summaries | **Nine-gene explicit effect-/feature- and prior-aware table + Figure; 7/9 present, 2/9 not assessable, no disease DGE** |

## 1. Independent source-level validation: 8,956 archived unique GWAS variants

**Native original**: `/srv/is-analysis/data/is/east_asia/japan/bbj/hum0197.v3.BBJ.IS.v1.zip`; original `GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz` member (nested gzip in ZIP).

- Fully streamed **13,435,541 original autosomal GWAS rows**.
- Collected all **1,224,224 gene–tissue–SNP records** from **646 independent *assay definitions***, representing **8,956 unique (BBJ locus, GRCh38 variant)** IDs and **8,956 unique native GRCh37 SNP IDs**.
- **8,956/8,956** IDs found once in original BBJ GWAS (no native duplicate IDs).
- Original `Allele2`, `AF_Allele2`, `BETA` matched archived genotype-oriented ALT allele, EAF and beta at all sites.
- Each native GRCh37 coordinate mapped through the existing **UCSC hg19→hg38 chain with a unique plus-strand mapped site**, and target GRCh38 `chr:pos:REF:ALT` agreed with the cached QTL join identity.
- All effect alleles for this set were original `Allele2=ALT`. The full-QTL cohort effect-allele beta has not been independently re-estimated from original genotypes; catalogue documentation states ALT is the effect allele.
- Separate actual *GRCh38 reference FASTA* validation remains a distinct boundary; chain agreement with QTL alleles does not itself certify the target reference base.
- **No causal gene claim**; this is a crucial *source-harmonization* gate, not shared causal signal validation.

Artifact:
`/srv/is-analysis/results/is/audits/is_8956_bbj_native_source_qc_20261011_v1/`
- `IS_8956_BBJ_NATIVE_EFFECT_ALLELE_CHAIN_QC.tsv`
- `IS_8956_BBJ_NATIVE_PROVENANCE_SUMMARY.json`

## 2. The GTEx ALT-major observation is an input selection/provenance question, not proof of flipped effects

The original 646-coloc archived inputs comprise **1,224,224 SNP-assay rows**, 8,956 unique locus SNPs (four BBJ regions). Every archived eQTL `ac/an` is ≤0.5 and `maf=min(ac/an,1-ac/an)`. The eQTL Catalogue **official field definition** specifies `ac = alternative allele count`, `maf = minor allele frequency`, and `alt = effect allele`. Its public data access FAQ also specifies ALT as the effect allele. Thus, minor AF alone cannot establish effect-beta orientation.

References:
- https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md
- https://www.ebi.ac.uk/eqtl/Data_access/
- https://www.gtexportal.org/home/faq

To differentiate catalog semantics from historical selection, the following **four existing original QTD regional extraction files** were read locally (no new QTL API calls): `QTD000131` artery aorta, `QTD000136` artery coronary, `QTD000141` artery tibial, `QTD000171` brain cortex, over GRCh38 chr12:111650000–112370000.

- **110,739 raw gene–variant records**, including **813 ALT-major records** with original `AC/AN>0.5`; their `MAF` values still correctly satisfy `min(AC/AN,1-AC/AN)`.
- Restricting to the gene IDs actually evaluated in the archived L003 coloc assays yields **549 original ALT-major gene–tissue–variant records representing 10 unique QTL sites**.
- None of these 549 exact original pairs (or any of the 10 corresponding GRCh38 loci) appeared in archived BBJ L003 coloc SNP-level source files.
- **This finding does not by itself prove a source pipeline bug.** The 10 may be absent from BBJ GWAS, fail upstream variant-level GWAS QC, fall outside the original included SNP lists, or be excluded by an explicit analysis design.

Artifacts:
`/srv/is-analysis/results/is/audits/is_gtex_qtd_alt_major_original_20261011_v2/`
- `IS_GTEX_QTD_L003_TISSUE_ORIGINAL_VS_COLOC.tsv`
- `IS_GTEX_QTD_ALT_MAJOR_ROW_GENE_MATCH.tsv`
- `IS_GTEX_QTD_SAVED_SOURCE_COVERAGE_MANIFEST.json`

A first intermediate `v1` directory had a partial TSV due to a summarization error when a category count was zero; the complete scientifically interpretable output is **`v2` only**.

## 3. Original BBJ GWAS resolves nine of the ten ALT-major omissions

All ten QTL variants were uniquely inverse-mapped via a **local UCSC chain** to GRCh37 (checked by the rs671 control `GRCh38 12:111803962 → GRCh37 12:112241766`).

A **fresh stream over all 13,435,541 native BBJ GWAS autosomal rows** then established:

- **9/10**: no original BBJ GWAS record at the mapped GRCh37 position. Their absence from BBJ/GTEx overlap may therefore be explained by native outcome GWAS coverage; it is **not proof of a harmonization error or intentional exclusion**.
- **1/10**: `GRCh38 12:111652218:C:A` ↔ `GRCh37 12:112090022:C:A` with directly verified native BBJ source REF/ALT, effect ALT. Original BBJ `AF_Allele2=0.997129291706149`, `p=0.161785997450046`. Its BBJ **minor allele frequency ~0.00287**, whereas GTEx ALT AF in the four saved original tissue files is **0.82687–0.86473**.
- **The last one-site exclusion is now resolved at the archived pipeline-input level.** Local 1KG EAS original `.pvar` includes `12:112090022:C:A` with `EAS_AF=0.999`, giving MAF ≈0.001. Historical PLINK `BBJ_IS_L003.QC.log` explicitly records `--geno 0.05 --maf 0.01`, zero variants removed for missingness, and 64,428 reference variants removed at the frequency threshold. This SNP is absent from the post-QC `.pvar` and the downstream `stage3_finemap/japan/bbj/susie_inputs_v3/BBJ_IS_L003.tsv` input (3,039 SNPs, minimum MAF ≈0.0100043). The original Colab preparation script uses this preselected Stage3-derived BBJ region input. This is an input-universe/MAF QC exclusion, **not evidence of an ALT direction flip**.
- **Never auto-insert this one marker to alter H4 post hoc.** An alternative lower-MAF analysis must be prespecified, independently quality-checked, and paired with suitable LD; it cannot be introduced selectively to increase H4.

Artifact:
`/srv/is-analysis/results/is/audits/is_gtex_10_altmajor_native_bbj_20261011_v1/`
- `IS_10_GTEX_ALT_MAJOR_SNP_NATIVE_BBJ_AVAILABILITY.tsv`
- `IS_10_GTEX_BBJ_COVERAGE_MANIFEST.json`

### Input lineage verification for the remaining native-BBJ marker

The existing Drive file [prepare_IS_COLAB_inputs.sh](https://drive.google.com/file/d/1QQP8oYzUD836WtPD3nbR7vOLtjpq_viw/view) points the Colab workflow to `/srv/is-analysis/results/is/stage3_finemap/japan/bbj/susie_inputs_v3/BBJ_IS_L003.tsv`; the existing [Phase10/11 repair Colab notebook v3](https://drive.google.com/file/d/1HM2j0kS_Twk5zBF0DWZPpafy-Xk4k0iu/view) consumes those preprocessed regional inputs. The **1KG EAS source** `BBJ_IS_L003.1KG_EAS.GRCh37.pvar` has the original site with ALT AF=0.999. In the saved Oct 2 PLINK QC log, `--maf 0.01` removes frequency-rare sites. The **post-QC reference pvar** and **3,039-row Stage3 summary-stat input** do not contain `12:112090022:C:A`, establishing a complete, existing-data-only **BBJ outcome → 1KG EAS LD MAF gate → Stage3 regional input → Colab coloc intersection** provenance trail.

Machine-audited derivative: `/srv/is-analysis/results/is/audits/is_gtEx_one_altmajor_maf_lineage_20261011_v1/IS_GTEX_ALTMAJOR_ORIGINAL_BBJ_LD_MAF_GATE.json` (`IS_GTEX_BBJ_LD_MAF_GATE_CONFIRMED`); source code `scripts/is/audit_is_altmajor_maf_gate_lineage.py` and fail-closed regression tests `tests/test_is_altmajor_maf_gate_lineage.py`.

The original input conversion [IS_Phase10_11_eQTL_Coloc_Master.ipynb](https://drive.google.com/file/d/1bHz26P9Y7BF2kCwO364OOH4t5kKVKrps/view) explicitly reads `INPUT/bbj_locus/{locus}.tsv` (the Drive export of `susie_inputs_v3`), performs GRCh37→38 conversion, then writes `WORK/{locus}.BBJ_GRCh38.tsv.gz`. The v3 repair notebook consumes these saved `WORK` intermediates and performs an exact-variant-key inner join with GTEx. This directly verifies the downstream path from the MAF-filtered precomputed SNP set to coloc, rather than assuming the repair notebook reads Stage3 files directly.

Do not generalize that ALL rare or missing GWAS variants were excluded by exactly the same filter. The 9 remaining variant identities are simply absent from the native BBJ outcome file; this 1 site has a positive PLINK ancestry-MAF QC provenance demonstration.

## 4. The nine legacy genes: do not confound healthy reference localization with colocalization proof

All nine historic gene symbols map uniquely to retained **stable Ensembl IDs** in the 2,225-gene positional universe, including `NEURL1` → `NEURL`, both `ENSG00000107954`. All 9 have at least one existing historical QTL assay; these are **not** nine causal genes and do not narrow the other **2,216** positional genes.

| Legacy gene | Max H4 at original p12=1e-5 | Adult healthy brain reference top annotated cell type | Human RNA feature |
|---|---:|---|---|
| CALHM2 | 0.79768 | Endothelial cells | Present |
| FGF5 | 0.77870 | Not assessable | **Not in matrix** |
| NEURL1 | 0.67759 | Astrocytes | Present |
| C4orf22 | 0.65994 | Not assessable | **Not in matrix** |
| INA | 0.53637 | Neurons | Present |
| SH3PXD2A | 0.35082 | Oligodendrocytes | Present |
| COL4A2 | 0.33078 | Smooth muscle cells | Present |
| COL4A1 | 0.11715 | Fibroblasts | Present |
| ALDH2 | 0.07904 | Endothelial cells | Present |

**7/9** target genes are measurable in the existing adult control brain reference; FGF5 and C4orf22 have `FEATURE_NOT_IN_MATRIX`, **not negative expression**. The reference uses **80,515 cells from six adult control donors**, tissue selected for vascular/perivascular contexts. This study does not test stroke patient vs control expression and is not a genetic eQTL disease mechanism experiment.

The nine-gene overlay includes H3/H4 original p12, H4 at strict/high prior for the SAME selected assay, original sample counts, top expression cell annotation and donor-pair QC counts. It does **not** average a GTEx brain tissue and a separate healthy temporal lobe reference into a causal score.

Artifacts:
`/srv/is-analysis/results/is/audits/is_nine_gene_gtEx_human_reference_20261011_v1/`
- `IS_NINE_GENE_QTL_PRIOR_HUMAN_DONOR_CONVERGENCE.tsv`
- `IS_NINE_GENE_MECHANISM_READINESS_MANIFEST.json`
- [Figure: nine-gene GTEx H4 prior sensitivity and healthy control reference](../../docs/figures/is/IS_NINE_GENE_PRIOR_H4_HEALTHY_CONTROL_REFERENCE.svg)
- [Reproducible SVG Figure source](../../scripts/is/render_is_nine_gene_qtl_healthy_reference_figure.R)

## 5. Original coloc reproducibility, prior fragility and retained broad scope

- **646/646** numeric source SNP-level `coloc.abf` replays passed; highest H0–H4 posterior discrepancy `6.72e-15`.
- On the original prior `p12=1e-5`, **0 of 646** assays reached `PP.H4>=0.8`; this does not depend on more generous prior post hoc selection.
- At `p12=1e-4`, **12 of 646** assays across **seven gene IDs** reached H4≥0.8. These are *prior-dependent*, not newly validated candidates.
- Human donor localization, GTEx coloc ABF, case/control GWAS association, matched-LD multi-signal causal inference and experimental functional validation are distinct layers.
- All **2,225 positional gene IDs / 2,425 gene-region relations remain present**, and archived input/canonical science gates have not been overwritten.
- ALDH2 **BLOCKED** in the alcohol fine-mapping branch; ADH1B **EXPLORATORY**. This GTEx-v8 ABF branch does not substitute for their required matched-LD scientific gates.

## 6. Code, regression and reproducibility

All code committed on GitHub branch `research/is-broad-discovery-20261009`:
- `scripts/is/audit_all_8956_bbj_native_allele_chain.py`
- `scripts/is/audit_saved_gtex_qtd_alt_major_overlap.py`
- `scripts/is/audit_gtex_alt_major_10_snp_native_bbj.py`
- `scripts/is/build_is_nine_gene_mechanism_readiness.py`
- `scripts/is/render_is_nine_gene_qtl_healthy_reference_figure.R`

New offline regression tests:
- `tests/test_all_8956_bbj_native_allele_chain.py`
- `tests/test_saved_gtex_qtd_alt_major_overlap.py`
- `tests/test_gtex_alt_major_native_bbj.py`
- `tests/test_is_nine_gene_mechanism_readiness.py`

**37/37 combined offline Python regression tests PASS** at commit `17e9dd8f7ee3b135a97b08852b4a00f939eb6f79`, including the four new exact-input MAF gate tests (server test report: `/tmp/is_646_8956_2225_all37_ci_20261011/RESULT.txt`). (includes earlier 646-coloc and 2,225-stable-ID guards). Latest GitHub Actions workflow is updated to run these. Remote GitHub Actions hosted execution has **not** been independently checked, and R replay runs require the original archived source files (not included in GitHub CI fixtures).

## Next steps that require no new databases

1. **Completed** the available source-lineage recovery: Drive original `prepare_IS_COLAB_inputs.sh` and Colab notebook identify the Stage3-derived input; original 1KG EAS PLINK QC log confirms `--maf 0.01`, the site MAF≈0.001, and its removal. Next evaluate possible implications of this prespecified MAF cutoff across loci without post hoc insertion or changing the original 646 replays.
2. Stratify **8,956 unique SNPs** by palindromic status and cross-population allele frequency without assuming reference-aware REF/ALT proof implies effect-allele source certainty across all molecular QTL.
3. Combine the nine-gene human donor evidence, coloc prior-grid and 2,225-gene broad discovery matrix into manuscript figures/results with clear `NOT_TESTED`, `EXPLORATORY`, `BLOCKED` states. Avoid ischemic-stroke disease DGE claims from adult control cells.
4. Separately inspect any externally reused BBJ/GIGASTROKE cohort overlap and conditional LD multi-signal requirements. Keep CKB and TPMI download blockers paused.

### Study interpretation

**Verified:** all 8,956 archived BBJ GWAS beta/EAF/ALT chain source identities, numerically reproduced 646 archival ABF tests, original QTL MAF vs AC/AN mathematics, 10-site native GWAS presence classification, and nine-gene descriptive control brain reference mappings.

**Not verified:** whether the GTEx eQTL source's ALT beta remains correctly referenced for all SNPs after prior processing, raw molecular-QTL genotype allele counts, harmonized cohort-matched molecular-QTL LD, independent East Asian GWAS replication, stroke lesion transcriptomic DGE, and final causal-gene assignments. The missing tenth ALT-major site has a **verified upstream MAF-QC exclusion lineage** (not a verified colocalization result).
