# IS Broad Discovery: Phase 2 molecular evidence work plan

Date: 2026-10-09
Branch: `research/is-broad-discovery-20261009`

## Verified inventory (real outputs)
- 6 currently analyzed BBJ/EAS GWAS canonical datasets, across 5 GIGASTROKE phenotypes.
- 45 source provisional regions collapse to **30 coordinate-overlap components**, NOT 30 independent loci.
- GRCh37 GENCODE v19: 869 gene/region associations, including 400 coding and 469 other gene IDs. A separate retained NEURL1 anchor outside current discovery windows brings v2 evidence to 870 rows.
- Legacy locus linkage by matching BBJ GWS intervals and same chromosome: L001->G0007, L002->G0015, L003->G0022, L004->G0023.
- Existing coloc linkages require BOTH original BBJ locus and stable Ensembl gene ID: 43 gene/region rows. Max ABF H4 is descriptive and never evidence of a causal gene by itself.
- Human targeted atlas assessed 9 symbols, mouse stroke descriptive pseudobulk assessed 4 genes.
- **826** positional gene/region rows have no matched historic coloc, so new molecular QTL screening is required; lack of evidence is not a biological zero.
- 30 GIGASTROKE accessions in official GWAS Catalog. EAS canonical present locally for 5; the other 25 are not locally available as EAS canonical. Registry UNKNOWN ancestry remains unresolved until per-study source verified.
- Lead reference alleles: 3 exact matches, 1 REF/ALT reversed, 26 lead positions absent from current regional 1KG EAS PVAR. These matches do not establish LD readiness.

## Reproducible pipeline
Run from the isolated IS worktree:

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
python3 scripts/is/run_broad_phase2_evidence.py
python3 -m unittest discover -s tests -p 'test_is_*' -v
```

Generated in `/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/`:
- `IS_LOCUS_GENE_EVIDENCE_V2.tsv`: locus- and Ensembl-specific evidence; 870 records
- `IS_ALL_CANDIDATE_MOLECULAR_WORK_QUEUE.tsv`: 869 gene-region molecular tasks
- `IS_MOLECULAR_TASKS_BY_REGION.tsv`: 30 region-level coverage records
- `GIGASTROKE_FULL_STUDY_INVENTORY.tsv`: exact accessions, labels and availability
- `IS_BROAD_PHASE2_REPRODUCIBILITY.json`: consistency assertions and SHA-256 fingerprints

## Critical next scientific gates
1. Obtain ancestry-matched genotype reference for non-anchor regions; harmonize REF/ALT and effect allele, inspect palindromic SNPs and variant overlap.
2. Obtain and normalize complete GWAS summary statistics from missing ancestry / meta-analysis source groups, with proper cohort-overlap documentation and build/allele harmonization.
3. Perform independent-signal assignment and per-ancestry fine-mapping. Do not rely on a megabase distance rule as proof of LD independence.
4. Source tissue-specific genome-wide eQTL/sQTL/pQTL summary statistics (brain, cerebral vascular/endothelial, blood/plasma as context). Match stable Ensembl gene IDs and variant IDs across same build; add n, EAF, effect sizes, QC/provenance.
5. Perform locus-specific ABF and multivariate/multi-signal colocalization where appropriate; assess PP.H3 vs PP.H4, priors, and signal ambiguity.
6. Expand human IS tissue and vascular cell-type datasets for all candidates. Normal human atlas and mouse MCAO model are supporting layers, not direct causal human stroke validation.
7. Include enhancer-linked distal genes outside the current 500 kb positional universe; retain the original 9 anchors, never pre-rank them solely because earlier work was deeper.
8. Only after data coverage audit, build an evidence-aware ranking with missingness indicators, replication and subtype specificity.

## Prohibited shortcuts
- No automatic causal gene claims from positional proximity.
- No conversion of NOT_TESTED to lack of biological association.
- No transferring locus-specific coloc H4 to another locus with the same gene symbol.
- No claim that 5 EAS files represent all 30 GIGASTROKE studies.
- No global use of the existing four BBJ-specific PGEN LD files.
