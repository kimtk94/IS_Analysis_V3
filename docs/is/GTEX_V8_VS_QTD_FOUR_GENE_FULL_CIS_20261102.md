# IS P3A — original GTEx v8 vs EBI QTD000161 four-gene cis variant universe audit

**2026-11-02 | CI [#41452103778](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41452103778) SUCCESS, 12/12 regression tests PASS.**

## Scope and provenance

- Original-imported GTEx v8 eQTL, Brain Cerebellar Hemisphere, supplied by official EMBL-EBI FTP: `imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz`. Original source MD5 `c07a1731eb96aafbd7aed131cb97515a` is present in the public official file list, but the *complete 3.4GB source bytes were not downloaded or rehashed*.
- EBI eQTL Catalogue reprocessed GTEx study QTS000015, dataset QTD000161 `QTD000161.all.tsv.gz`; exact historical release-6 compressed bytes have not been independently pinned, even though the *gene-level tested variation denominators* match the older 2026-10-11 research audit.
- GRCh38 coordinate windows retained from the independent historical study: chr4:79,000,000–81,500,000 and chr10:101,500,000–106,500,000.
- Two official sources accessed by bounded HTTPS Tabix, filtering to FGF5 ENSG00000138675, C4orf22/CFAP299 ENSG00000197826, CALHM2 ENSG00000138172, NEURL/NEURL1 ENSG00000107954.
- Explicit REF and ALT matched under exact gene-ID + GRCh38 pos + allele sequence; includes A/C/G/T-coded SNPs, indels, MNVs (max 200bp) instead of silently dropping non-SNP rows. Unsupported alleles, invalid p/β/SE and duplicated rsid/source keys are counted and screened. A source with incompatible duplicate statistics or excessive invalid rows is HOLD; data transfer is bounded at 1 million total source records/window.
- No independent original GTEx v8 full-cis denominator from a permutation file; do NOT label exact original-cis completeness as approved solely from these two intervals.

## Actual full-reprocessed-cis comparison

| Gene | Original-imported variants in region | EBI reprocessed full cis variants | Exact common gene+position+REF/ALT | Original-only | EBI-only | Beta direction agreement | Median abs β diff | Median abs SE diff | Median abs log10(p) diff |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FGF5 | 6,831 | 5,511 | 5,510 | 1,321 | 1 | 3,965/5,510 (71.96%) | 0.11438 | 0.05835 | 0.38433 |
| C4orf22 | 6,821 | 5,490 | 5,490 | 1,331 | 0 | 4,315/5,490 (78.60%) | 0.06560 | 0.00372 | 0.18713 |
| CALHM2 | 5,678 | 4,868 | 4,868 | 810 | 0 | 3,824/4,868 (78.55%) | 0.07728 | 0.02386 | 0.24333 |
| NEURL | 5,633 | 4,820 | 4,819 | 814 | 1 | 3,065/4,819 (63.60%) | 0.08578 | 0.05308 | 0.30960 |
| **Total gene–variant pairs** | **24,963** | **20,689** | **20,687** | **4,276** | **2** | **15,169/20,687 (73.33%)** | — | — | — |

The 20,689 reprocessed-variant denominator reproduces the previous independent four-gene QTD reprocessing report **exactly**: FGF5 5,511, C4orf22 5,490, CALHM2 4,868, NEURL 4,820. Original source has 4,276 more pairs absent from reprocessed QTD, while two QTD variants (one FGF5 and one NEURL) are not matched in original. Their exact identity/cause has **not** yet been investigated; do not attribute to model errors.

**Unit of these counts is gene–variant association, not an independent SNP/participant.** Same locus can contribute two genes and LD makes p-values/signs correlated. Do not perform a binomial test of β sign rates using nominal counts as independent observations.

The previous 19,155 common SNP-only pairs were a strict subset; the *non-SNP* 1,532 additional exact matches account for the move to 20,687. The earlier EBI denominator HOLD was caused by treating all non-single-base alleles as invalid. Retaining valid explicit ACGT indels/MNVs corrects this under an auditable allele-length and encoding rule. No automatic strand complementation or β sign flipping is performed.

## Science gate

**PASS for the following bounded question only:** Can the same 20,689 QTD full-cis tested gene–variant records be traced to an original-imported GTEx-v8 model under exact GRCh38 genomic/allele keys? 20,687 pairwise matches; two documented exceptions; exhaustive count independently reconciled.

**HOLD:** historical inherited 2026-10-11 QTL source and exact QTD r6 bytes, source-specific original GTEx v8 *full* cis variant denominator, original GWAS per-variant N/case fraction, QTL effect allele model equivalence, suitable matched-genotype LD, variant-missingness impact, p12 sensitivity and any MR, SuSiE, ABF coloc or causal gene claims.

A full pairwise mapping is **not a replication**: the GTEx-v8 original and EBI QTD both reuse underlying GTEx donors. Do not describe as independent cohorts. The results show systematic differences in effect estimates that may arise from model/normalization/covariate differences, not automatically data corruption. Carry each single-source QTL as an independent sensitivity track.

## Reproducible evidence

- Pipeline `scripts/is/audit_gtex_v8_qtd_four_gene_fullcis.py`
- Negative controls `tests/test_is_gtex_v8_qtd_four_gene_fullcis.py`
- GitHub CI `.github/workflows/is-gtex-v8-vs-qtd-fourgene-fullcis.yml`
- Actual CI result, proven source run and per-gene counters: [research/is/evidence/gtex_v8_vs_qtd_four_gene_fullcis_20687_live_20261102.json](../../research/is/evidence/gtex_v8_vs_qtd_four_gene_fullcis_20687_live_20261102.json)

**Next most targeted work:** identify two QTD-only variant keys, determine if original source has the same genomic position with different REF/ALT representation, and independently authenticate source-specific cis tested variant boundaries and QTL/LD study model before allowing formal coloc. Keep draft PR unmerged.


## Actual two-unmatched-QTD-variant provenance screen — 2026-11-02

[GitHub Actions #41455461447](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41455461447) **SUCCESS**, 16/16 deterministic regression tests PASS. The runner re-read the original and reprocessed official GRCh38 indexed source files across the original four-gene cis windows. Total source/QTD/matched gene–variant counts are unchanged at **24,963 / 20,689 / 20,687**.

| Gene | QTD-only SNP (GRCh38 position REF:ALT) | β | SE | p-value | Original same-gene exact position | Original same-gene positions ±5bp |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| FGF5 | chr4:81,266,618 T:G | -0.102885 | 0.0912328 | 0.261125 | 0 | 0 |
| NEURL / NEURL1 | chr10:102,493,750 T:C | -0.0133939 | 0.0899905 | 0.881870 | 0 | 0 |

Both unmatched records are **single-nucleotide variants, not indel/MNV discrepancies**. There was no eligible same-gene original source variant with a matching exact position, reversed REF/ALT order at that position, or a nearby position within ±5 bp. This supports a **source tested-variant set difference**, not a simple coordinate/allele notation change. However, the precise upstream reason (minor allele frequency, variant QC, missing genotype call, donor inclusion, original file mapping) was **not** proven by this test. The QTD-only p-values are non-significant in this dataset; they should not be called negative causal results.

Of the original-only **4,276** gene–variant records, **3,959 are single-nucleotide variants** and **317 are explicit indel/MNV records**, counted separately in the live audit. As another separate source QC observation, the official original-imported source contained **41 target-gene records with p=1, β=0 and a nonfinite SE encoded as -nan** (6 chr4 + 35 chr10); these were excluded from the strict valid-effect comparison. Their positions do not match either of the two above QTD-only variants in the bounded diagnostic. Therefore these source-QC exclusions are not evidence of a simple repair for the two unmatched QTD records.

### Effect allele interpretation

The [official eQTL Catalogue nominal association column definitions](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md) state explicitly that `alt` is the **effect allele**, and `ref` is the GRCh38 reference allele. The [official GTEx v8 methods](https://gtexportal.org/home/methods) likewise define normalized effect size `NES` / slope with respect to **ALT** relative to REF. Thus the original/reprocessed comparison uses documented matching effect-allele definitions, and large differences in β sign among exact REF/ALT matched variants **cannot simply be explained by swapping effect allele labels in one source**. This is not proof that raw-model coefficients have the same normalization/covariate scale: inherited pipelines remain scientifically distinct and LD-correlated gene-variant rows do not supply independent replication.

### Remaining boundaries

1. Original GTEx v8 *independent full-cis tested-variant interval boundaries* and source-specific `n_variants` metadata are not yet independently corroborated by the original eGene/permute tables.
2. Exact original source full-file hash, historical EBI release-6 file bytes and inherited local GTEx-QTL source identity are unverified; verify before promoting reproductions to definitive.
3. FGF5 and NEURL source-specific variant filtering mechanisms are unverified; no correction should synthesize or impute the two missing original QTL associations.
4. GWAS outcome-specific N/case fraction, ancestry-matched LD, independent EAS outcome replication and formal MR/colocalization remain scientific HOLD.

Actual diagnostic JSON: [research/is/evidence/gtex_v8_vs_qtd_2_unmatched_source_diagnostics_20261102.json](../../research/is/evidence/gtex_v8_vs_qtd_2_unmatched_source_diagnostics_20261102.json). The source-audit code includes exact/nearby variation and bounded discarded-row audit with negative controls.


## Independent GTEx-v8 GENCODE v26 TSS geometry — 2026-11-02

Actual source job [#41456646711](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41456646711), **SUCCESS, 7/7 tests PASS**, downloaded the official **GENCODE v26 gene GTF** (37,607,127 gzip bytes, SHA-256 `e84157bb3ddfe6b1dd0b8816ec71b553074ed162ab0940ac48b80375f6cb18bd`), parsed four gene features and tested the TSS ±1Mb cis definition published by GTEx v8.

| Gene | GENCODE v26 strand | TSS, GRCh38 | TSS ±1Mb boundaries | QTD-only SNP | Distance beyond nominal cis |
| --- | --- | ---: | --- | --- | ---: |
| FGF5 | + | 4:80,266,599 | 79,266,599–81,266,599 | 4:81,266,618 T:G | **19bp above upper** |
| C4orf22 | + | 4:80,335,720 | 79,335,720–81,335,720 | none | not applicable |
| CALHM2 | - | 10:103,452,402 | 102,452,402–104,452,402 | none | not applicable |
| NEURL | + | 10:103,493,979 | 102,493,979–104,493,979 | 10:102,493,750 T:C | **229bp below lower** |

All four 2Mb GENCODE v26 gene-feature cis proxies are contained within the historical broad chr4/chr10 screening windows. **Both QTD-only SNPs fall immediately outside their respective nominal original GTEx v8 TSS±1Mb proxies**, providing an independently sourced, mechanistically coherent explanation for both primary-source QTL gaps. Given that neither SNP had the matching same-gene original row at ±5bp and both QTD sources have MAF >30%, boundary-based differences are more strongly supported than a rare-variant explanation. This is a comparison of the annotation rule; it does **not** prove the exact input genotype or full cohort exclusion filter.

Crucial distinction: GTEx v8 used a **collapsed GENCODE v26 gene model**, and we have not yet independently confirmed that its collapsed exon-model TSS matches the regular GENCODE v26 gene-feature TSS. The original GTEx [public reference file](https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf) was explicitly identified in [the GTEx Consortium reproducibility instructions](https://github.com/broadinstitute/gtex-v8), but a read-only public GET in CI [#41457128784](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41457128784) returned **HTTP 403**. No requester-pays download, payment or credentials were attempted. Therefore **actual GTEx collapsed model TSS remains HOLD**, even though the standard GENCODE TSS proxy and exactly observed source exclusions are compelling.

Official source and observed distances saved in
[research/is/evidence/gtex_v8_gencode26_tss_two_qtd_only_20261102.json](../../research/is/evidence/gtex_v8_gencode26_tss_two_qtd_only_20261102.json). Reproducible CI and tests:
`scripts/is/audit_gtex_v8_gencode26_cis_geometry.py`,
`tests/test_is_gtex_v8_gencode26_cis_geometry.py`,
`.github/workflows/is-gtex-v8-gencode26-tss-cis-audit.yml`.

Further progress without paid sources can use the [official GTEx Portal v2 API](https://gtexportal.org/api/v2/redoc) endpoint `/dataset/collapsedGeneModelExon` with `datasetId=gtex_v8`, plus an authenticated GENCODE v26 versioned gene ID. Explicit API-derived collapsed-exon TSS should be independently compared with the two positions before upgrading the original cis-boundary gate.

**Scientific boundary:** Two-source QTL gene–variant mappings and full QTD screened counts PASS; original FastQTL exact collapsed TSS, historical release byte identity, original complete per-gene test boundaries, GWAS case fraction, ancestry LD, MR/coloc remain HOLD.


## Exact GTEx v8 collapsed gene-model TSS independently confirmed — 2026-11-02

Source: official [GTEx Portal API v2 collapsedGeneModelExon](https://gtexportal.org/api/v2/redoc), **explicit datasetId=gtex_v8** (the current API default is gtex_v10), with versioned GENCODE v26 identifiers extracted from the official GENCODE v26 GTF. The official API requests were performed **sequentially**, not in parallel. Real [GitHub Actions #41460057633](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41460057633): **SUCCESS, 9/9 deterministic tests PASS**.

| Gene | GENCODE v26 versioned ID | GTEx v8 API collapsed exons | GTEx v8 collapsed-exon strand-aware TSS (GRCh38) | Difference to regular GENCODE v26 gene-feature TSS |
| --- | --- | ---: | ---: | ---: |
| FGF5 | ENSG00000138675.16 | 6 | chr4:80,266,599 | 0bp |
| C4orf22 / CFAP299 | ENSG00000197826.11 | 12 | chr4:80,335,720 | 0bp |
| CALHM2 | ENSG00000138172.10 | 5 | chr10:103,452,402 (−) | 0bp |
| NEURL / NEURL1 | ENSG00000107954.10 | 9 | chr10:103,493,979 | 0bp |

**Key independently checked boundary evidence**:

- FGF5 QTD-only `chr4:81,266,618:T:G` lies **19bp beyond** the v8 collapsed TSS+1Mb upper limit `chr4:81,266,599`.
- NEURL QTD-only `chr10:102,493,750:T:C` lies **229bp below** the v8 collapsed TSS−1Mb lower limit `chr10:102,493,979`.
- Neither has an original-imported GTEx v8 gene–SNP association at the same position or ±5bp in the actual original indexed data, independently verified by [#41455461447](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41455461447).
- GTEx v8 original models and GENCODE v26 standard gene features agree in the TSS for **all four genes**, removing the model-TSS-drift concern as a plausible explanation for the two outside-window positions.

This confirms the exact **GTEx API collapsed exon model TSS**, not the separately stored **original FastQTL tested cis gene BED**. The evidence is therefore much stronger for a **cis window definition difference** than for an unexplained REF/ALT swap or low MAF in reprocessed QTD, but source-genotype/QC differential mechanisms and original FastQTL BED provenance are still not directly authenticated.

Evidence frozen under
[research/is/evidence/gtex_v8_official_api_collapsed_tss_qtd_exclusions_20261102.json](../../research/is/evidence/gtex_v8_official_api_collapsed_tss_qtd_exclusions_20261102.json). Reproducible CI `.github/workflows/is-gtex-v8-api-collapsed-tss.yml`, tested script `scripts/is/audit_gtex_v8_api_collapsed_exons.py` and tests `tests/test_is_gtex_v8_api_collapsed_exons.py`. Old direct public cloud file URL returned 403; this separate official GTEx API provided the necessary reference geometry without requester-pays downloads.

**Result gate:** GTEx-v8 API collapsed TSS agreement: PASS; 2 QTD-only variant outside-boundary location: PASS; exact original FastQTL tested cis BED and all other MR/coloc scientific gates: HOLD.
