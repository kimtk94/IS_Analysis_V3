# IS G0022 — Multi-tissue eQTL Catalogue audit and publication evidence gate
**Research date:** 2026-10-10 KST
**Branch:** `research/is-broad-discovery-20261009`
**Source provenance:** original eQTL Catalogue indexed GTEx gene-level GRCh38 nominal QTL tables, plus previously verified GIGASTROKE EAS AIS GRCh37 GWAS.

## Why we expanded to four tissues
The full-locus G0022 AIS SuSiE-RSS signal (4,649 ref-QC SNPs; 1000G EAS n=504 LD; study-level case/control approximate N_eff≈70,474) yielded **one exploratory 95% credible set of four variants**. Single-window analyses had previously produced three credible sets, but full-locus modeling reduces that to one. All are exploratory, not validated causal signals.

The first GTEx Artery Aorta nominal QTL panel could not assay any of the four high-posterior GWAS variants. We extended actual raw QTL extraction and allele harmonization to artery coronary, artery tibial and brain cortex. No canonical data or the 2,225-gene IS candidate universe was overwritten.

## Data provenance and real audit numbers

Source `eQTL-Catalogue/eQTL-Catalogue-resources` tabix dataset manifest: `QTS000015` (GTEx), quantification `ge`, GRCh38, `chr12:111650000-112370000`. Source dataset sample counts below are metadata study sample sizes, not verified per-SNP genetic effective N.

| Dataset | Tissue | Sample size (catalogue) | Regional raw association rows | Exact GRCh37 EAS AIS gene–SNP allele matches |
|---|---|---:|---:|---:|
| QTD000131 | Artery Aorta | 387 | 27,504 | 3,858 |
| QTD000136 | Artery Coronary | 213 | 27,046 | 3,858 |
| QTD000141 | Artery Tibial | 584 | 34,490 | 3,858 |
| QTD000171 | Brain Cortex | 205 | 21,699 | 2,796 |
| **Total** | **4 tissues** | — | **110,739** | **14,370** |

These are *association rows / gene–SNP comparisons*, not unique genome-wide independent SNP signals. They were retrieved as individual, small targeted tabix region extracts, not bulk datasets (per eQTL Catalogue server access guidance). Every source file was gzip-tested and SHA256-stamped in `G0022_MULTITISSUE_QTL_SOURCE_MANIFEST.tsv`.

Methods:
- Official UCSC hg38ToHg19 chain and `pyliftover` with one-to-one, GRCh38 to GRCh37 zero-based interval conversion, alternative effect allele tracked explicitly;
- GTEx QTL Catalogue `ref`/ `alt` is genomic REF/ALT, and `alt` is the effect allele (verified against official `tabix/Columns.md`);
- GIGASTROKE EAS AIS `alt_effect_beta` matches the ALT-coded GWAS association estimate;
- exact `chr:pos:REF:ALT` join; palindromic SNPs excluded; ancestry allele frequencies not naively treated as identical;
- 8 genes audited: ALDH2, ACAD10, NAA25, HECTD4, BRAP, PTPN11, IFT81, ATP2A2; the remaining candidate gene universe was not restricted or discarded.

For the six genes with any common QTL/GWAS SNPs: 643 shared SNPs per gene in each artery sample, and 466 per gene in brain cortex; IFT81 and ATP2A2 lack matched SNPs in this 720kb extract (NOT_TESTED).

## Important discovery: *none* of the 4 GWAS credible-set variants is in the QTL regional input

We downloaded the official `hg19ToHg38.over.chain.gz` from UCSC and performed GRCh37→GRCh38 liftover *and reciprocal mapping back* for every GWAS CS SNP:

| GWAS CS GRCh37 | GRCh38 position (chr12) | 1KG EAS ALT AF | EAS GWAS ALT EAF |
|---|---:|---:|---:|
| 12:112168009:G:A | 111730205 | 0.1756 | 0.2319 |
| 12:112241766:G:A | 111803962 | 0.1736 | 0.2311 |
| 12:112468206:C:T | 112030402 | 0.1716 | 0.2208 |
| 12:112736118:A:G | 112298314 | 0.1607 | 0.2253 |

All four GRCh38 positions are within the region explicitly queried in every GTEx dataset. We searched **all original nominal QTL records, across all genes and alleles**, for each position (4 variants × 4 tissues = **16 source-position tests**).

**16/16 tests = NOT_TESTED_VARIANT_POSITION_ABSENT.**

This is not a strand or allele-label harmonization issue; the positions themselves are absent from these exact eQTL Catalogue regional source extractions. It does *not* prove GTEx donors do not carry the variants, and it does not establish biological absence of an eQTL; variant filtering, genotype availability, study ancestry and assay coverage are unresolved.

## Exploratory coloc.abf numerical diagnostics (DO NOT interpret as shared causal effect)

To verify software and data contracts, processed 32 prespecified gene–tissue comparisons. 24 comparisons had enough common SNPs to run numerical `coloc.abf`; eight had no shared QTL/GWAS SNPs. All 24 numerical posteriors are publication-gated as **BLOCKED: absent main GWS/CS variants**. Not a single validated independent causal gene.

For traceability only, the largest diagnostic PP.H4 by tissue:
- Artery Aorta: ALDH2 PP.H4 = 0.3814.
- Artery Coronary: HECTD4 PP.H4 = 0.0719.
- Artery Tibial: ALDH2 PP.H4 = 0.5615.
- Brain Cortex: NAA25 PP.H4 = 0.3381.

None qualifies as positive colocalization. The PP.H4 posterior is sensitive to the **selected variant universe** and coloc's one-causal-signal assumption, and 720kb QTL window excludes the critical GWS/CS variants. Passing a numerical posterior threshold would not resolve that defect.

## Fail-closed publication gate

`G0022_MULTITISSUE_QTL_PUBLICATION_GATE_SUMMARY.json` records:
- 4 target tissues, 32 gene–tissue analyses
- 110,739 original regional QTL association records
- 14,370 exact allele-matched gene–SNP records
- 24 diagnostic posterior computations
- **0 credible-set SNPs observed in these QTL extracts**
- **0 scientifically valid colocalizations**
- **0 validated causal genes**

Code actively throws if the source coverage count, expected tissue denominator or publication boundary unexpectedly changes. Re-assay or new credible set evidence must be audited before promotion.

## Reproduction

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# 1) Only if missing, download each official QTD000xxx whole 720kb REGION (one at a time):
tabix 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000136/QTD000136.all.tsv.gz' '12:111650000-112370000' | gzip -1 > /srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1/QTD000136_artery_coronary_GRCh38_chr12_111650000_112370000.tsv.gz
# Other files already downloaded as QTD000131, QTD000141, QTD000171.
python3 scripts/is/prepare_g0022_eqtl_catalogue_multitissue.py
R_LIBS_USER=/srv/is-analysis/.Rlib Rscript scripts/is/run_g0022_eqtl_multitissue_abf.R
python3 scripts/is/audit_g0022_qtl_cs_assay_coverage.py
python3 scripts/is/audit_g0022_multitissue_qtl_gate.py
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

Outputs:
- Server input archive: `/srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1/`
- Research evidence: `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/eqtl_catalogue_multitissue_v1/`
- Master: `G0022_MULTITISSUE_QTL_PUBLICATION_GATE_SUMMARY.json`
- Source coverage: `G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE.tsv`
- Separate diagnostic PP.H4: `G0022_GTEX_MULTITISSUE_ABF_DIAGNOSTIC.tsv`
- Per-source SHA256 and sample metadata: `G0022_MULTITISSUE_QTL_SOURCE_MANIFEST.tsv`

## Next research action

**Do not spend more compute repeating GTEx regions with the same underlying SNP coverage just to obtain a coloc number.** Prioritize finding additional molecular-QTL cohorts that actually assay the four critical AIS variants (ideally East Asian ancestry), or obtain a sufficiently rich cis-QTL summary-statistics set with complete variant coverage. Then perform true multi-signal molecular fine-mapping and `coloc.susie` under appropriate ancestry-matched LD and sensitivity tests.

In parallel, obtain vascular and neurovascular cell-state expression/ATAC enhancer-to-gene functional evidence for the *entire* 2,225-gene universe, without interpreting a missing tissue/assay record as a biological negative.

## Official references
- eQTL Catalogue source/access: https://www.ebi.ac.uk/eqtl/Data_access/
- eQTL Catalogue column definitions: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md
- UCSC official GRCh37→GRCh38 chain: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
- Original GIGASTROKE: https://www.nature.com/articles/s41586-022-05165-3
