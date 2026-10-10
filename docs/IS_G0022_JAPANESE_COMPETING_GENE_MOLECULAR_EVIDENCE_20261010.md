# IS G0022 — Four-variant Japanese immune-cell QTL / JCTF competing-gene audit

**Analysis:** 2026-10-10 KST · Branch: research/is-broad-discovery-20261009
**Scientific gate:** SOURCE_ASSOCIATIONS_VERIFIED, FULL_GWAS_QTL_COLOC_NOT_PERFORMED, MULTIVARIATE_MEDIATION_NOT_IDENTIFIED, VALIDATED_CAUSAL_GENE_COUNT=0.

## 1. Question, rationale and source-level mistake prevented

Prior G0022 effort correctly prioritized functional rs671/ALDH2, but focusing on ALDH2 alone can mask **other genes in the same very-high-LD East Asian stroke region**. ImmuNexUT's ALDH2 gene-centric page had displayed the four candidate SNPs, yet it could not reveal other gene targets of those SNPs. A new four-variant-centric query was made against the original ImmuNexUT **significant-eQTL browser pages**, then compared with existing source-audited JCTF Japan Omics variant-level eQTL/pQTL.

Source pages were fetched *once each*, stored unchanged with SHA256, and parsed only after checking exact rsID and GRCh38 position matching pre-verified GIGASTROKE EAS AIS G0022 GRCh37 four-SNP CS:
- rs11066015 ⇔ chr12:112168009 G>A (GRCh37) ⇔ chr12:111730205 (GRCh38)
- rs671 ⇔ chr12:112241766 G>A ⇔ chr12:111803962
- rs11066132 ⇔ chr12:112468206 C>T ⇔ chr12:112030402
- rs77768175 ⇔ chr12:112736118 A>G ⇔ chr12:112298314

## 2. ImmuNexUT genuine browser association matrix

Four small real variant pages contain **14 significant gene–variant associations** over four unique genes (4+4+3+3). No separate private genotype data downloaded. Browser statistics represent the most significant immune-cell association per gene–SNP and **not all 28 cell types**.

| Gene | Four-CS variant coverage | Representative cell | reported minimum browser p (same within gene) | Unharmonized browser beta | source interpretation |
|---|---:|---|---:|---:|---|
| **BRAP** | 4 of 4 | Neutrophils / Neu | **8.22219×10^-21** | -0.362395 | Significant marginal eQTL |
| **RPH3A** | 4 of 4 | Classical monocytes / CL_Mono | **5.53234×10^-10** | -0.335162 | Significant marginal eQTL |
| **ALDH2** | 4 of 4 | Neutrophils / Neu | **4.88802×10^-9** | -0.207695 | Significant marginal eQTL |
| **MYL2** | 2 of 4 | Plasmablasts | **7.56457×10^-9** | -0.30853 | Significant marginal eQTL at 2 markers; absence at other 2 is NOT biological negative |

**Important:** BRAP has a smaller nominal eQTL p than ALDH2 in ImmuNexUT, but the genes differ in assay, expression level, cell specificity, and variance. P-value order is **not causal-gene rank**. Four correlated markers do not represent four independent biological QTLs. The reported beta lacks individually validated effect-allele orientation for the AIS GWAS, so negative signs are not comparable with AIS estimates.

The JPT104 subgroup LD to rs671 is nearly perfect for the top three G0022 SNPs, and 1000G EAS n504 LD rank≤503 for 4,649 SNPs remains exploratory. The current candidate set also includes locus distal genes and the 2,225-gene discovery universe remains intact.

## 3. Direct cross-reference to the existing real JCTF Japan Omics QTL file

Source is the existing 72-row JCTF four-variant report, a significant/PIP-filtered release. This is an **orthogonal source/assay comparison**, not a meta-analysis or independent human stroke replication.

| Gene | ImmuNexUT 4-CS eQTL coverage | JCTF 4-CS eQTL | JCTF 4-CS Olink pQTL | Max JCTF SuSiE variant PIP among four | Main caveat |
|---|---|---|---|---|---|
| ALDH2 | 4; Neu | **4**; min p=3.51e-12 | **4**; min p=6.59e-49 | eQTL **0.354**, pQTL **0.140** | Coding rs671 may alter Olink assay recognition; transcript vs protein vs enzyme-mediated stroke unproven |
| BRAP | 4; Neu | **0 observed in filtered JCTF** | **4**; min p=2.61e-5 | pQTL **0.0767** | ImmuNexUT eQTL and JCTF pQTL different assays; no validated shared causal signal |
| RPH3A | 4; CL_Mono | **4**; min p=1.61e-21 | **0 observed in filtered JCTF** | JCTF eQTL **0.000** for the tested GWAS SNPs | Extremely small marginal eQTL p is **not fine-mapped causal membership** |
| MYL2 | 2; Plasmablasts | **0 observed in filtered JCTF** | **0 observed in filtered JCTF** | Not available | No observation can reflect JCTF release threshold or assay coverage, NOT biological absence |

All JCTF entries are from the *filtered* QTL release (p<0.05 OR PIP>0.001) and the given PIPs refer to QTL-side variant models, not to the AIS GWAS SuSiE. No posterior can be fabricated from two marginal p-values. JCTF and ImmuNexUT study donor overlap and source-specific imputation/calibration remain unaudited.

**New conclusion:** ALDH2, BRAP and RPH3A all have source-linked molecular associations in Japan-related datasets, through different molecular traits. BRAP and RPH3A must remain explicit G0022 mechanistic competitors and should not be suppressed by the rs671 coding annotation. This strengthens multi-gene locus annotation but does **not** make three causal genes.

## 4. Full-nominal archive and alternate QTL source access decision

- ImmuNexUT FAQ explicitly identifies the *all-nominal including non-significant* dataset as E-GEAD-420.
- DDBJ 2026-10-10 HTTP source inspection: E-GEAD-420.processed.zip exact size 43,507,433,191 bytes (official MD5 listed), one nested file nominal_eQTL_2.tar.gz, itself about 43.50 GB inside the ZIP; compressed outer ZIP entry is deflate method 8 and the inner TAR is gzip-compressed. ZIP central entry cannot individually index SNP/ALDH2 from the nested tarball. **No bulk download started**.
- Official E-GEAD-398 contains *conditional/significant* eQTLs, 564,882,709-byte ZIP including conditional_eQTL_FDR0.05.tar (uncompressed size 4,133,191,680). It is not the full nominal test universe, and would not unlock unbiased new coloc by itself.
- NBDC NHA000166 has 141 Japanese donors, but its methods explicitly focus on **miRNA cis-eQTL**; it is not an alternative ALDH2 *mRNA* full-cis input. No unnecessary download attempted.
- JCTF v2 and v3 remain thresholded; MAEEA remains restricted per prior server audit (live timeout). Large NARD2 reference genotype/signed LD not locally available.
- **Next data acquisition goal:** obtain an approved **gene/cell-specific export of complete nominal E-GEAD-420 eQTL** including effect allele, beta, SE, N, all SNPs for ALDH2/BRAP/RPH3A/MYL2 and complete variant assay coverage. Alternatively a BGZF+tabix/parquet re-export is needed before region-level source processing.

The decision to avoid 43.5GB bulk download is operational, not scientific rejection of E-GEAD-420.

## 5. Required identification before any publishable ALDH2 claim

1. GWAS G0022 full-locus EAS AIS: per-SNP N / available larger ancestry- and cohort-matched signed LD, LD uncertainty and sample overlap.
2. Molecular source: unfiltered tested SNP–gene universe in matching immune/neurovascular tissues, actual effect allele, sample-specific beta/se, SNP frequency/INFO and signed matched LD.
3. Compare at least **ALDH2, BRAP, RPH3A and other positional/distance genes** at equal coverage and gene-specific power. Compute conditional/signal-level QTL credible sets and valid coloc (coloc.susie or equivalent with appropriate LD priors).
4. Independently validate AIS effect in disjoint EAS stroke GWAS. No claim that JCTF pQTL/eQTL is independent stroke replication.
5. Separate ALDH2 enzyme loss, alcohol intake/tolerance, blood pressure and aldehyde toxicity paths; assess BRAP/RPH3A perturbation evidence. An Olink missense pQTL may be assay-binding artifact.
6. If tissue/variant source unavailable, declare *marginal immune-cell molecular evidence and scientific limitations* rather than promoting any causal gene.

## 6. Reproduction and data provenance

Git code:
- scripts/is/audit_g0022_immunexut_competing_genes.py
- tests/test_is_g0022_immunexut_competing_genes.py

Immutable HTML source archive:
- /srv/is-analysis/data/is/qtl/immunexut/G0022_v1/variant_pages_20261010/{rs11066015,rs671,rs11066132,rs77768175}.html
- Exact SHA256 for each file in output summary JSON (source-page CSRF tokens make fresh HTML hashes vary, so snapshots are not silently overwritten).

Reference JCTF TSV:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/jctf_japan_omics_variants_v1/G0022_JCTF_4SNP_JAPANESE_EQTL_PQTL_EVIDENCE.tsv

Analysis output folder:
- /srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/immunexut_competing_genes_v1/
- G0022_IMMUNEXUT_4CS_ALL_GENE_VARIANT_ASSOCIATIONS.tsv (14 rows)
- G0022_IMMUNEXUT_JCTF_COMPETING_GENE_MATRIX.tsv (four-gene ledger)
- G0022_IMMUNEXUT_JCTF_COMPETING_GENE_SUMMARY.json (source hashes and fail-closed science status)

## Official links
- ImmuNexUT FAQ and conditional browser scope: https://www.immunexut.org/faqs
- rs671 variant-centric significant associations: https://www.immunexut.org/eqtlSnpsLink?snp_id=rs671
- ImmuNexUT cell names: https://www.immunexut.org/abbreviations
- DDBJ E-GEAD-420 archive filelist: https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/E-GEAD-420.filelist.txt
- DDBJ conditional E-GEAD-398: https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-398/E-GEAD-398.filelist.txt
- NHA000166 (miRNA): https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000166

**Safety:** original GIGASTROKE/BBJ summary stats, old QTL outputs, prior SuSiE/coloc, unrelated branches and production services unchanged.
