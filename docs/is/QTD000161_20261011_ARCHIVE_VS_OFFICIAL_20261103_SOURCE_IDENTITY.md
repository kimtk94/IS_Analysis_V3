# IS — 2026-10-11 archived QTD000161 source vs 2026-11-03 official data provenance

**2026-11-03 (KST) | Historic four-gene QTL source and reprocessing version continuity**

## 1. Recovered actual 2026-10-11 source from connected Drive

Original file: `IS_PRIORITY4_EBI_REPROCESSED_FULL_CIS_VARIANTS.tsv`, Drive ID `12g0HrusD3ONV2v9ormnCDMfLNGdXqibU`, byte size **2,928,061**, **20,689 gene-variant records**, exact SHA-256 **`a80d5a6903aea551d5ee6f8bebf12772aea411a7cd29437bd545d5d2452768a3`**. Independently recomputed this in the working runtime and compared to the saved 2026-10-11 `IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGA_MANIFEST.json` source SHA.

Its actual archived values comprise `gene_id`, `chromosome`, `position_GRCh38`, `REF_GRCh38`, `ALT_GRCh38`, `QTD000161_nominal_QTL_p`, `QTD000161_beta`, `QTD000161_se`, plus historic source lineage flags. The numerical comparison was frozen in two minimal tracked panels from this exact input TSV:
- chr4: **11,001** gene–variant rows, SHA-256 `e34ec9d1200ce68838e39174a0c54f550f1b8d99af7d4eac9c1904d6f2f2555e` at `research/is/evidence/is_priority4_archived_qtd000161_chr4_source20261011.tsv`.
- chr10: **9,688** gene–variant rows, SHA-256 `a493c7d8ef32bb142e45e7775d288226bede0a095d001dfc5ba499c30cee7779` at `research/is/evidence/is_priority4_archived_qtd000161_chr10_source20261011.tsv`.

The CI checks these two **exact frozen SHA-256 values before reading association rows**. It does not silently regenerate a reference from the current live EBI FTP.

## 2. Actual 20,689-row independent original-to-live QTD comparison

Source [GitHub Actions #41553486065](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41553486065): **SUCCESS, 10/10 tests PASS**. Downloaded *only* two predeclared chr4/chr10 ranges from the current official EBI QTD000161 gene expression nominal summary file via indexed HTTPS tabix; compared exact stable ENSG ID + GRCh38 chromosome, coordinate, REF and ALT followed by signed ALT-beta, SE, p-value to archived records.

| Gene | Original archived records | Exact shared keys with current EBI | beta, SE, p concordant within prespecified precision |
| --- | ---: | ---: | ---: |
| FGF5 | 5,511 | 5,511 | **5,511/5,511 PASS** |
| C4orf22 | 5,490 | 5,490 | **5,490/5,490 PASS** |
| CALHM2 | 4,868 | 4,868 | **4,868/4,868 PASS** |
| NEURL | 4,820 | 4,820 | **4,820/4,820 PASS** |
| **Total** | **20,689** | **20,689** | **20,689/20,689 PASS** |

No missing, extra or numerically discordant gene–variant rows; 0 bounded mismatch examples. Tolerances for rounding differences: beta and SE relative `5e-5`, absolute `1e-6`; when p>0, absolute difference of log10 p ≤`1e-4`. This is **numeric/source consistency in the full four-gene archived QTD set**, NOT byte identity of the whole multi-gene `QTD000161.all.tsv.gz`.

Real CI log captured in `research/is/evidence/qtd000161_archived20261011_vs_official_live20689_20261103.json`.

## 3. Archived QTD permutation gzip — exact source byte identity PASS

The original 2026-10-11 `IS_PRIORITY4_EBI_FULL_CIS_RECONSTRUCTION_MANIFEST.json` (Drive ID `1dFLd6tcfq8E0oK6Wk44abM12IrAdIInH`) records an independently saved compressed permutation source hash:

`QTD000161.permuted.tsv.gz` SHA-256 **`b5de9603b25b7f65f5532442bb5cc90f5d546421387363633967b540b4424bbe`**.

Real [GitHub Actions #41553961484](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41553961484) **SUCCESS, 9/9 tests PASS** retrieved the **entire original compressed file** from current official EBI FTP, **831,384 bytes**, recomputed the full SHA-256, and found **EXACT byte identity** with the stored historic SHA-256. Its permutation summary contained 23,287 gene records. Four historical gene-level permutation test counts and `p_beta` matched:

| Gene | n_variants in historical and current permutation source | historical = current p_beta |
| --- | ---: | ---: |
| FGF5 | 5,511 | 0.0278074 |
| C4orf22 | 5,490 | 0.000591983 |
| CALHM2 | 4,868 | 0.0270452 |
| NEURL | 4,820 | 0.40605 |

These `p_beta` values are permutation-based beta-distribution p-values. **They are not FDR-corrected eGene q-values nor evidence of PP.H4**. Source schema: [eQTL Catalogue official Columns.md](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md).

Actual evidence: `research/is/evidence/qtd000161_permuted_20261011_to_live_SHA256_20261103.json`.

## 4. Explicit scientific decision

| Question | Gate |
| --- | --- |
| Historical derived four-gene QTL TSV independently matches prior saved SHA | **PASS** |
| Archived gene+SNP QTL beta, SE, p agrees with current EBI QTD | **20,689/20,689 PASS** |
| Archived entire QTD permutation compressed source matches current gzip | **SHA256 identical PASS** |
| QTD nominal `.all.tsv.gz` entire historical gzip bytes confirmed | **HOLD** |
| Full EBI 2023 release6 source processing vs late-2026 release model version provenance certified | **HOLD** |
| Historical inherited `GTEx_V8` archival source identical to canonical original FastQTL model | **HOLD** |
| GWAS AIS per-variant N/case fraction, matched ancestry LD and actual MR/coloc | **HOLD** |

The exact archived QTD-derived panel and current EBI QTD nominal associations are reproducibly consistent at four genes. Permutation gzip is byte-identical to the historically cached one. **This strengthens archived input provenance but does not establish an unchanged study-wide release or automatically validate historical coloc posteriors.**

## 5. Reproducible GitHub-only implementation

- `scripts/is/audit_qtd000161_archived_20261011_vs_live.py`
- `tests/test_is_qtd000161_archived_vs_live.py`
- `.github/workflows/is-qtd000161-historic-archive-vs-live.yml`
- `scripts/is/audit_qtd000161_historic_permuted_sha.py`
- `tests/test_is_qtd000161_historic_permuted_sha.py`
- `.github/workflows/is-qtd000161-historic-permutation-sha.yml`

No SSH server, paid GCS, full GTEx original 3.4GB download, source overwritten or merged PR. Next meaningful gate: authenticate the inherited historical `GTEx_V8` variant-level source independently against official original GTEx v8 (rather than by its directory name); assess exact AIS N/case fraction assumptions and QTL matched LD before coloc.
