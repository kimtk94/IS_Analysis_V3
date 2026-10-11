# IS EAS5 — public QTL ↔ GWAS exact-SNP GRCh37 allele harmonization audit

**2026-10-11 KST · `research/is-broad-discovery-20261009`**
**Primary result: VERIFIED INPUT VARIANT OVERLAP; NO NEW COLOCALIZATION OR CAUSAL GENE CLAIMS.**

## 1. Research question

Can the new 15 eQTL Catalogue nominal-QTL regional source files from five provisional East-Asian ischemic-stroke genomic groups be joined to their exact BBJ Japanese ischemic-stroke and GIGASTROKE EAS AIS GWAS SNPs, without assuming gene proximity is functional causality?

The original **2,225 positional genes**, **869 EAS region-linked genes**, **43 genes with archived BBJ–GTEx ABF**, and **646 historical ABF assays** remain unchanged. The five regions are **distance-based provisional genomic groups**, not 5 confirmed independent LD-defined causal loci.

## 2. Original data and computational science gates

- Original public nominal QTL (GRCh38, 19 columns): 5 groups × 3 datasets = **15/15 verified raw regional extracts**; **2,677,112 total gene–SNP regional records** including genes outside our positional list. Exactly **1,622,101** records belong to original five-group positional candidate gene IDs.
- QTL datasets: **GTEx v8 Artery Aorta QTD000131 (catalogue N=387)**, **GTEx Brain Cortex QTD000171 (N=205)**, **ROSMAP brain DLPFC QTD000434 (N=560)**. Source ancestry composition and genotype LD for those molecular-QTL donors were **not** independently confirmed.
- Correct **GRCh38 → GRCh37** conversion uses `/srv/is-analysis/data/reference/hg38ToHg19.over.chain`. Its UCSC chain header identifies target chromosome 1 length **248,956,422** (GRCh38) and query chromosome 1 length **249,250,621** (GRCh37). Reverse strand variants complement REF and ALT but do **not** reverse regression beta sign when the same ALT biological allele remains the effect allele.
- Verify single unique mapping, exact mapped **GRCh37 REF against `human_g1k_v37.fasta`**, biallelic SNP status, GWAS variant genomic ID, and GWAS effect/other allele orientation. Exclude unresolved or conflicting pairs instead of inferring swaps from frequency or unproven palindromic strand.
- The **official eQTL Catalogue documentation** explicitly states `ALT` **is the effect allele**, `AC` is its count, `AN` is total allele count: [official nominal eQTL column definitions](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md) and [EBI Data Access FAQ](https://www.ebi.ac.uk/eqtl/Data_access/). QTL beta is therefore interpreted **relative to the original source ALT**. Only SNPs with verified REF/ALT and coordinate transfer can be compared to GRCh37 ALT-oriented GWAS effect.
- GWAS processed/canonical inputs (originally from GRCh37 source): **BBJ 13,435,531 rows** and **EAS AIS 6,725,216 rows**; only the selected five ±500kb exploratory windows were indexed. In those windows, **29,353 BBJ** and **8,752 EAS AIS** source SNVs passed independent FASTA REF/effect allele/numeric QC. Source AIS GWAS metadata specifies **256,274 total participants** but does **not** break out case/control sample counts, which are required for independently calibrated case-control coloc setup.

## 3. Executed gene–SNP intersection

| Gate, record-level (gene × SNP × molecular-QTL dataset) | Records |
|---|---:|
| Original candidate gene QTL nominal source records | **1,622,101** |
| Biallelic, unique GRCh38→37 liftOver with GRCh37 FASTA REF concordance | **1,332,650** |
| Exact BBJ ischemic-stroke GWAS variant match | **805,037** |
| Exact GIGASTROKE EAS AIS GWAS variant match | **427,259** |
| Same mapped SNP in both BBJ and EAS AIS GWAS | **423,051** |
| Both GWAS, nonpalindromic A/T or C/G allele pair excluded | **376,334** |

**Important unit:** these are *gene–SNP–QTL-dataset records*, not independent SNPs or novel gene discovery counts. EBI can repeat records for multiple rsIDs, so source-level same-gene/same-SNP duplicates were checked explicitly. Within the matched derivative, **809,245 unique gene–SNP–dataset–region records** were available from one or both GWAS, **0 conflicting/duplicate rsID records** in this selected dataset after exact-variant verification. Actual duplicate behavior can differ in other QTL sources.

Among the **809,245** one-or-both GWAS matched records, **100,870** were palindromic; no automated direction inference or allele-frequency flip was applied for them.

Source-frequency heterogeneity warning: among SNP–gene records matching EAS AIS, **211,866** had absolute **ALT frequency difference ≥0.10** versus the QTL donor sample, and **93,024** had difference **≥0.20**. Those differences are not prima facie strand errors (GWAS and QTL REF/ALT were sequence-checked). The true explanation—ancestry composition, LD, genotyping QC, or residual source errors—requires study-specific validation. These are a **QC signal, not automatic exclusion or causal evidence**.

## 4. Five-region readiness partition

| East Asian provisional region | Original positional gene IDs | Genes with public QTL nominal rows + BBJ/EAS AIS matched SNPs | Newly source-eligible genes absent archived ABF | Both-GWAS exact **nonpalindromic** gene–SNP source records |
|---|---:|---:|---:|---:|
| `IS_XDATA_G0015` (T1) | 71 | 46 | **30** | **118,328** |
| `IS_XDATA_G0023` (T1) | 14 | 9 | **3** | **31,786** |
| `IS_XDATA_G0007` (T1) | 10 | 4 | **1** | **7,910** |
| `IS_XDATA_G0025` (T4) | 92 | 61 | **61** | **143,363** |
| `IS_XDATA_G0018` (T4) | 53 | 29 | **29** | **74,947** |
| **Five-group total** | **240** | **149** | **124** | **376,334** |

The **124 novel QTL-source genes** are not novel causal genes; they were previously *unscreened in the archived 646 ABF* and now have exact GWAS–QTL SNP overlap. No archived ABF input was reanalysed or modified. The original **826 EAS genes with no archived ABF** still number **826**, even though 124 now have new nominal QTL sources.

## 5. Crucial distinction: separate minimum P values vs SAME SNP association

A permissive preparation heuristic identified **36 gene–tissue source combinations**, all in `G0015`, with at least 100 overlapping nonpalindromic SNP rows plus a *minimum* AIS P <5e−8 and a *minimum* nominal QTL P <1e−4 anywhere among their overlap SNPs. However **those two minima need not be the same variant**.

An independent, strict *same exact physical SNP* check examined all **393** gene–tissue records with source/GWAS intersection:

- Original QTL `P ≤1e−4`, EAS AIS `P <5e−8`, exact matching GRCh37 SNP and REF/ALT, nonpalindromic requirement.
- **0/393 gene–tissue comparisons** contained a SNP passing *both* thresholds simultaneously.
- In G0015, **43 gene–tissue source comparisons** had some QTL-significant matching SNP, and **81** had some AIS-significant matching SNP, but **0** had an exact SNP passing both thresholds.

This does **not** falsify all colocalization: trait-specific tagging due to LD, variable power and single vs multiple causal variant models can produce different marginal maxima. Equally, simply observing two significant minima in a broad locus does **not** support a shared causal signal. The previous **36** comparisons are **EXPLORATORY SCREENING HEURISTIC ONLY** and must not be described as 36 causally supported eQTLs.

## 6. Science status and next action

| Gate | Status | Interpretation |
|---|---|---|
| Public nominal QTL full input availability for 5 pilot regions | **VERIFIED** | 15 SHA256-verified archived raw extracts |
| GRCh38 QTL SNP → GRCh37 FASTA identity / GWAS exact allele pair | **VERIFIED for included rows** | 376,334 nonpalindromic gene–SNP source records in both GWAS |
| Official QTL ALT-effect beta orientation and AC/AN ALT frequency | **DOCUMENTED AND APPLIED to verified rows** | No ambiguous automatic SNP flip |
| EAS AIS case fraction and per-dataset genotype/QTL population LD | **NOT VERIFIED** | Case-control coloc and multi-signal assumptions unresolved |
| QTL donor ancestry, expression measurement power, full cis SNP missingness and study independence | **NOT VERIFIED** | Strong ALT frequency heterogeneity warrants source audit |
| Same-physical-SNP strict joint marginal P threshold | **0/393** | Descriptive threshold only; not proof of no colocalization |
| New coloc.abf, coloc.susie, MR, causal gene promotion | **NOT RUN** | No new disease mechanism claims |

**Next defensible analysis:** recover study source/sample phenotype details and actual QTL cis tested-universe and molecular ancestry/LD, independently document AIS case/control fraction. Assess local LD mismatch and sample independence, then consider baseline `p12=1e−5` **predeclared exploratory** single-signal coloc on fully harmonized SNP sets with separate sensitivity. A high prior alone must never be used to promote a negative outcome to shared causal signal. Existing G0022/ALDH2 separate mechanistic branches remain independent and **BLOCKED** where the relevant validation gate is unresolved.

## 7. Reproducible artifacts, no source overwrite

Original input scans and QC:
- `scripts/is/audit_is_eas5_public_qtl_gwas_snp_overlap.py`
- `/srv/is-analysis/results/is/audits/is_eas5_qtl_gwas_snp_ref37_20261011_v2/`

Official ALT alleles/AC–AN QC and de-dup:
- `scripts/is/audit_is_eas5_official_alt_effect_and_dedup.py`
- `/srv/is-analysis/results/is/audits/is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1/`

Same physical SNP threshold + allele frequency:
- `scripts/is/audit_is_eas5_same_variant_qtl_ais_p.py`
- `/srv/is-analysis/results/is/audits/is_eas5_same_variant_p_af_20261011_v1/`

Figure: [Readiness by source-tested gene and exact nonpalindromic BBJ/EAS AIS SNP overlap](../figures/is/IS_EAS5_QTL_GWAS_SNP_OVERLAP_READINESS.svg).

Small published derivative tables:
- `docs/data/is/IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv`
- `docs/data/is/IS_EAS5_SAME_SNP_ASSOCIATION_GENE_TISSUE_QC.tsv`
- `docs/data/is/IS_EAS5_SAME_SNP_ASSOCIATION_BY_REGION_QC.tsv`
- `docs/data/is/IS_EAS5_SOURCE_QTL_GWAS_SNP_OVERLAP_QC.tsv`
- `docs/data/is/IS_EAS5_ALT_BETA_PROVENANCE_MANIFEST.json`

**Safeguards:** initial invalid `v1` source join was deliberately **quarantined** (`INVALID_RUN_DO_NOT_USE.json`); the corrected `v2` audit uniquely uses the **hg38ToHg19** chain, validated by file-name and GRCh38/37 chromosome-size header tests. Raw QTL 15 regional slices, BBJ/GIGASTROKE source GWAS, original 646 archived coloc model inputs, and candidate lists remained unchanged. Large row-level files (16 MB gz matched, 269 MB official ALT-matched TSV) are **retained on the server only**, not copied to Git.
