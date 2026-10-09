# IS G0022 molecular QTL: GTEx v8 DAP-G and eQTL Catalogue nominal summary statistics
**Research snapshot:** 2026-10-09 KST
**Branch:** `research/is-broad-discovery-20261009`
**Interpretation:** STRICTLY EXPLORATORY, no validated molecular colocalization/causal gene.

## Background and purpose

Whole-locus (≈4.05Mb) East-Asian AIS G0022 SuSiE-RSS now converges to **one purity-filtered 95% candidate credible set, 4 SNPs**, not three independent clumps. Source GWAS GIGASTROKE EAS AIS `GCST90104545`, aligned GRCh37, 4,649 genotype-QC SNPs, with 1000 Genomes EAS n=504 external LD; case/control effective N is approximate, not SNP-specific. This is **not validated fine-mapping**, and all pQTL/eQTL evidence must remain ancestry- and study-specific.

CS variants (GWAS SuSiE exploratory only):
- `12:112241766:G:A` PIP 0.37484 (GENCODE v19 gene body ALDH2).
- `12:112468206:C:T` PIP 0.29776 (NAA25).
- `12:112168009:G:A` PIP 0.27247 (ACAD10).
- `12:112736118:A:G` PIP 0.03778 (HECTD4).

PIP is not a biological causal-gene nomination.

## 1. GTEx v8 DAP-G fine-mapping signals — allele-level screen

Source: official GTEx Portal v2 API:
`https://gtexportal.org/api/v2/association/fineMapping` with `datasetId=gtex_v8`, `gencodeId=...`, and paginated capture for 8 prespecified genes: ACAD10, ALDH2, NAA25, HECTD4, BRAP, PTPN11, IFT81, ATP2A2.

GTEx uses **GRCh38**, original GIGASTROKE is **GRCh37**. Converted each *biallelic SNP* with UCSC `hg38ToHg19.over.chain` and `pyliftover 0.4.1`, zero-based coordinate conversion and reverse complement on minus-strand. Required unique one-to-one mapping and **exact chr:position:REF:ALT match**, did not infer ambiguous allele flips or transfer Ensembl gene-version numbers across GENCODE releases.

Real-world audit:

| Candidate gene | GTEx v8 DAP-G records | SNPs matched to G0022 AIS genotype |
|---|---:|---:|
| ACAD10 | 11 | 9 |
| ALDH2 | 768 | 521 |
| NAA25 | 377 | 305 |
| HECTD4 | 171 | 122 |
| BRAP | 0 | 0 |
| PTPN11 | 0 | 0 |
| IFT81 | 70 | 7 |
| ATP2A2 | 89 | 22 |
| **Total** | **1,486** | **986 records** |

Liftover outcomes: 1,355 unique-mapped source records; 131 non-biallelic/other-chromosome variants excluded; 369 mapped but not in the restricted genotype-QC SNP universe. Of 986 matched QTL DAP-G records, **0** directly overlap the four GWAS CS SNPs. Negative overlap is *not* negative colocalization: GTEx fine-mapped variant sets do not contain all assayed variants, populations differ, only 8 genes were scanned, and the GWAS CS is not yet validated.

Tissue-stratified screen: 36 gene-tissue pairs with matching records (4 arterial gene-tissue pairs and 7 brain gene-tissue pairs), still no direct 4-SNP CS overlap.

## 2. eQTL Catalogue — original all-nominal variants in GTEx Artery Aorta

Source: EMBL-EBI eQTL Catalogue, source dataset **QTD000131 (GTEx artery_aorta, gene-level, n=387)**, `GRCh38`, region **chr12:111650000–112370000**. Retrieved from official indexed BGZF TSV using one targeted tabix-region request (follow service rate-limit guidance; REST API is deprecated).

```bash
tabix 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000131/QTD000131.all.tsv.gz' '12:111650000-112370000' | gzip -1 > QTD000131_artery_aorta_GRCh38_chr12_111650000_112370000.tsv.gz
```

Original region extraction `QTD000131_artery_aorta_GRCh38_chr12_111650000_112370000.tsv.gz`, **614,326 bytes**, SHA256 `99749308ee1957d880c0928f70ef1fb167a3de772d1d387fca39881f9fcd2330`.

Parsing uses the documented eQTL Catalogue 19-column schema, including position/REF/ALT, gene_id, beta, SE, MAF, allele count and p-value. rsID duplicates collapsed **only if identical SNP-gene effects**; distinct alleles, ambiguous liftover and strand ambiguity are blocked. Separate ancestry-specific GWAS `alt_effect_beta` orientation preserved.

Actual processing:
- **27,504** regional QTL association rows;
- **7,524** records for the 8 targeted genes;
- **3,858** distinct matched gene-SNP records (643 each for 6 genes);
- **306** palindromic records excluded from this target slice;
- **2,604** QTL records without an exact eligible GWAS allele pair;
- **756** non-biallelic/unsupported type records excluded.

**Important coverage deficit:** among every target gene's 643 common SNPs, best AIS GWAS p-value is **3.399e-5**; the full-locus AIS GWS lead and other strongest GWAS signals are absent from this shared genotype/QTL SNP subset. The nominal extract is also limited to only 720kb around the top CS; it does not represent the complete cis-QTL universe. This is an explicit **NO_GWS_COMMON_SNP** evidence gate.

## 3. coloc.abf diagnostic (NOT scientific evidence of colocalization)

To verify the pipeline runs against real QTL and GWAS beta/SE data, applied `coloc::coloc.abf` to 6 gene-specific 643-SNP subsets. Hypotheses H0–H4 use priors p1=1e-4, p2=1e-4, p12=1e-5, **single causal signal assumption**; GWAS AIS case-control `N_total=256274`, `case_fraction=19032/256274` (not verified SNP-specific); GTEx artery aorta gene-level `N=387`. All variants required exact genome REF/ALT alignment.

| Gene | Common SNPs | Min nominal QTL p | Min AIS p | Numerical PP.H4 (diagnostic only) | Scientific gate |
|---|---:|---:|---:|---:|---|
| ALDH2 | 643 | 2.691e-5 | 3.399e-5 | **0.3814** | **BLOCKED** |
| HECTD4 | 643 | 0.01551 | 3.399e-5 | 0.0972 | **BLOCKED** |
| BRAP | 643 | 0.01914 | 3.399e-5 | 0.0650 | **BLOCKED** |
| NAA25 | 643 | 0.10188 | 3.399e-5 | 0.0591 | **BLOCKED** |
| ACAD10 | 643 | 0.06345 | 3.399e-5 | 0.0582 | **BLOCKED** |
| PTPN11 | 643 | 0.14593 | 3.399e-5 | 0.0554 | **BLOCKED** |
| IFT81 / ATP2A2 | 0 | NA | NA | NA | **NOT TESTED** |

**DO NOT promote PP.H4=0.3814 for ALDH2 as positive coloc.** Shared SNP set lacks GWS signal, the region is truncated, independent GTEx donor population and LD differ from Japanese/EAS meta GWAS, and coloc.abf single-causal assumption is particularly unsuitable for originally multi-clump G0022. This exercise is computational QC only. Validated molecular colocalizations **0**.

Previous `BBJ_IS_L003` GTEx coloc best H4: ALDH2=0.07904 and PTPN11=0.11315. These are separate BBJ studies and **not directly compared or pooled** with current AIS coloc diagnostic.

## 4. Reproduce without modifying canonical results

Worktree `/srv/is-analysis/worktrees/is-broad-discovery-20261009`.

```bash
python3 scripts/is/audit_g0022_gtex_v8_dapg.py --execute
python3 scripts/is/summarize_g0022_gtex_v8_tissues.py
python3 scripts/is/prepare_g0022_eqtl_catalogue_aorta_coloc.py
R_LIBS_USER=/srv/is-analysis/.Rlib Rscript scripts/is/run_g0022_eqtl_catalogue_aorta_coloc.R
python3 -m unittest discover -s tests -p 'test_is_*' -q
```

The first command accesses a public API and caches exact gene-specific JSON pages in the output folder. The QTL Catalogue input must be separately retrieved once by tabix as above. Both outputs are under `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/`, in separate `gtex_v8_dapg_liftover` and `eqtl_catalogue_aorta_v1` directories. `pyliftover` installed in isolated, untracked reference directory. The canonical 80 discovery groups and 2,225 Ensembl positional genes are unchanged.

## 5. Next scientific gate

The strongest GWAS GWS markers must be directly represented in the full nominal eQTL tested SNP universe **before** interpreting coloc. Evaluate aorta/coronary/brain cortical QTL and QTL test availability for index GWAS SNPs. Add full cis-region signal coverage, allele frequency concordance, matched ancestry LD, variant/sample N and QTL fine-mapping, then `coloc.susie` with independent signals. This may require alternate public QTL cohorts because original GTEx/GIGASTROKE East Asian allele coverage is incomplete. Never assume no overlap = absence of causality.

References:
- GTEx Portal API documentation: https://gtexportal.org/api/v2/redoc
- eQTL Catalogue access: https://www.ebi.ac.uk/eqtl/Data_access/
- eQTL Catalogue nominal columns: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md
- GIGASTROKE: https://www.nature.com/articles/s41586-022-05165-3
