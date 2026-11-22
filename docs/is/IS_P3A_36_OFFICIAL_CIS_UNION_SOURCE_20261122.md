# IS P3A — Three original QTL sources, full 20-gene Ensembl 105 reference cis union (2026-11-22)

**MASTER DEGREE / ischemic stroke P3A.** 36 genuine \`region_id × QTL_dataset × gene_id\` assays, 20 unique genes, **one** chr10 G0015 GWAS locus. This is **official QTL source recovery**, not a formal colocalization result or genetic causality claim.

## Evidence and reproducibility

- [Official Ensembl 105 phenotype position/source MD5 and prior 2.05Mb coverage run #45027864903](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45027864903), **SUCCESS**, evidence \`P3A_36_Ensembl105_reference_cis_window_evidence_HOLD.json\` (the original Ensembl105 record is [Zenodo 4715946](https://zenodo.org/records/4715946), source MD5 \`b999b21dcffeecbfe0bf886a5a91997b\`).
- [Prior authentic QTL baseline, type-disaggregated 2.05Mb run #45027864947](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45027864947), **SUCCESS**, \`P3A_36_QTL_official_GRCh38_indexed_2Mb_context_HOLD.json\`.
- Official cohort/source release provenance QTD000131/171/434 and study metadata: [main CI #45024493082](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45024493082). The real original UCSC GRCh37→GRCh38 6-point chain and verified GRCh38 chromosome 10 FASTA came from [#45024963210](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45024963210).
- **New actual indexed-source three-QTL 3.52Mb union CI:** [#45028586923](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45028586923), **SUCCESS**, focused synthetic tests **13/13 PASS**, full repository tests **300 executed, 299 PASS, 1 SKIP**.
- Original official eQTL Catalogue methodology describes ±1Mb cis-window around the gene start and gene-expression QTLs in uniformly reprocessed RNA data: [2023 eQTL Catalogue article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538656/), [qtlmap pipeline](https://github.com/eQTL-Catalogue/qtlmap). A reference positional window alone **does not prove** which variants were tested for each individual cohort and molecular trait.

## Declared source windows and observational reproducibility

| Interval | GRCh38 chr10 1-based inclusive | Length |
|---|---|---:|
| Former shared region | 102,403,919–104,460,593 | 2,056,675 bp |
| **Full union of the 20 Ensembl105 reference gene cis windows** | **101,837,472–105,353,833** | **3,516,362 bp** |
| Added left interval | 101,837,472–102,403,918 | 566,447 bp |
| Added right interval | 104,460,594–105,353,833 | 893,240 bp |

Exactly **one (1) bounded, remote indexed GRCh38 range request per dataset** using pinned PySam/HTSlib 0.24.1; 9-second sequential throttle, no automatic retries, no download of the ~11GB of complete original QTL BGZF files. Each source's original 2.05Mb data and gene×variant valid counts must be a strict **subset** of the expanded query; missing previously observed source records block CI. Only aggregate counts/digests/source identities are retained in the uploaded Actions artifacts.

### Actual source readback per data set

| Official QTL source | Tissue (catalogue N, not approved gene N) | Old 2.05Mb candidate gene×variant rows | New 3.52Mb candidate gene×variant rows | New eligible nonpalindromic SNP×gene records | Source INDEL/multibase excluded by SNP-only policy | New region total including non-candidates |
|---|---|---:|---:|---:|---:|---:|
| QTD000131 | GTEx artery aorta (N=387) | 46,737 | **57,780** | **44,693** | **5,916** | **241,786** |
| QTD000171 | GTEx brain cortex (N=205) | 42,647 | **58,823** | **45,570** | **5,985** | **264,947** |
| QTD000434 | ROSMAP brain DLPFC (N=560) | 74,501 | **100,822** | **58,081** | **32,553** | **366,252** |
| **Total** | 36 gene×tissue records, not 36 independent loci | **163,885** | **217,425** | **148,344** | **44,454** | **872,985** |

**Net additional evidence after expanding the query:**
- **+53,540** original nominal QTL candidate gene×variant rows (**217,425 vs 163,885**), an increase of **32.7%**.
- **+37,166** unique within-gene nonpalindromic SNP×gene assay records (**148,344 vs 111,178**), an increase of **33.4%**.
- **33/36 gene×tissue groups** gained new original source associations; the three CALHM2 assays correctly stayed constant because their reference cis domain was already wholly contained in the former query.
- Among **all 36 groups**, **0 source QTL association records fell outside their respective Ensembl105 gene-centered ±1Mb reference windows** during this expanded query. This is stronger positional compatibility evidence with the original QTL mapping convention, but does not establish original cohort-specific phenotype metadata identity or enumerate every input genotype.
- For canonically formatted single-base non-palindromic SNPs, **0 reference FASTA mismatches**, **0 invalid finite beta/SE/MAF observations**, and **0 duplicate conflicting effects** were observed in the streamed reference region. **44,454 legitimate type=INDEL source rows** were excluded by the previously declared SNP-only mapping policy, not declared corrupt or biologically irrelevant.
- All **36/36** authentic original gene×tissue assay groups retained observed source records; the actual study groups and their 20 gene IDs were not changed.

These are **SNP×gene×tissue records, not 148,344 unique genomic SNPs**. A SNP may contribute to multiple genes and tissues. The original 3.52Mb source scan SHA256 checksums describe text returned by the **Tabix region iterator**, not complete publisher-file hash digests:
- \`QTD000131\`: \`451ba8ba15493986f431b1751aff9ea512ab3dbc04ef41e38f95217d49cf3f53\`
- \`QTD000171\`: \`5964ba14a274f59247eebf88602da21706fdaea6c6e5df57990a3157f7854953\`
- \`QTD000434\`: \`89e6950472b8563442747eafc63fdbe5e75aa27fb28e3f0fec6ede5fc114ffa2\`

## 36-key independent technical QC and science blockers

\`scripts/is/audit_p3a_official_cis_union_36.py\` rejects:
- Any changed 36 source identities or a missing original gene among the 20.
- Source-provenance/MD5 or ref-FASTA SHA drift from pinned authentic results, or changed original GRCh38 region/Ensembl105 reference union.
- Truncated original physical 19-column associations, off-chr10 records or improperly formatted beta/SE/MAF, REF vs official FASTA, SNP allele or variant keys.
- A new source region that loses any previously observed gene-assay association or valid nonpal SNP.
- Treating 20-gene reference ±1Mb positional coverage as proven full study-specific tested nominal cis variants, inferred per-test effective N, or causality.

**Scientific acceptance state remains unchanged, even when all source checks succeed:**
1. Original QTL **study-specific phenotype metadata and the exact source-tested full cis SNP sets** have not been independently accepted. Source QTL tissue N values 387/205/560 are catalogue dataset metadata, not necessarily per-gene effective N.
2. QTL GRCh38 vs original GIGASTROKE GRCh37 **all SNP effect-allele/other-allele beta/SE/EAF, strand, multiallelic and INDEL normalization and liftover have not been independently reconciled**. The original G0015 C/G lead remains strand-palindromic without independent orientation confirmation.
3. GIGASTROKE EAS AIS exact source GWAS **cases, controls and case fraction, cohort overlap (including BBJ), phenotype definition and valid ancestral matching** are not independently confirmed, and source genome-wide summary alone cannot be reported as independent replication.
4. \`coloc.abf\` eligibility **0/36**, calculations **0**. \`coloc.susie\` eligibility **0/36**, calculations **0**. Do not create H4/H3 results from a partial or unapproved QTL-GWAS merge.

**Next research gate:** Recover per-study QTL phenotype source identity and GIGASTROKE EAS AIS outcome effective sample sizes/cases/controls, perform **variantwise** GRCh37↔GRCh38 canonical reference allele liftover with per-assay complete SNP QTL–GWAS overlap and ancestral frequency checks; retain formal \`HOLD\` for any missing assumptions. No fabricated N_eff, imputed SE or "independent BBJ replication" labels.
