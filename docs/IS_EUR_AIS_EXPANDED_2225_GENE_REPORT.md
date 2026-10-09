# IS ancestry-expanded candidate universe — verified EAS/Japanese + EUR AIS
2026-10-09, `research/is-broad-discovery-20261009`

## Real data provenance

- Source European AIS GIGASTROKE `GCST90104540`, GRCh37, ancestry EUR, official GWAS Catalog source-MD5: `1bfe4ae8a5042fb24cb0a562b05b0d2f` **verified byte-for-byte** after full download.
- Source sample size in official YAML: **1,296,908** (study-wide source metadata; not asserted as per-SNP effective N).
- Normalized EUR AIS summary statistics: **7,482,032** valid source variant records, preserving effect allele, other allele, beta, SE, P and EAF. All are SNP pair keys, **NOT** asserted genomic REF/ALT orientations.
- EUR AIS observed **1,054** GWS (P≤5e-8) and **1,540** suggestive-or-better (P≤1e-6) SNP records with MAF≥0.01 in the provisional discovery analysis.
- **50 EUR coordinate-proximity clusters** at 1Mb separation: 26 GWS-containing, 24 suggestive-only.
- **5** EUR clusters overlap positional intervals of the 30 previous BBJ/EAS components. The two source sets remain separate to protect ancestry-specific conclusions.
- Original discovery: 30 EAS/Japanese component records, 869 positional Ensembl gene IDs and gene-region rows.
- New manifest: **80 ancestry-specific study/interval records** (30 existing + 50 EUR), **not 80 independent loci**.

## GENCODE v19 GRCh37 annotation result

Broad discovery with gene-overlap ±500kb windows generated:
- **2,425** gene-region assignments: **869** original and **1,556** EUR-region associations.
- **2,225 unique gene IDs**, containing **994 protein-coding** and **1,231 noncoding/other**.
- 1,356 more distinct positional genes than original 869-gene study universe.
- Original nine legacy candidate genes retained. Eight overlap current positional windows; `NEURL1` stays as an explicitly preserved anchor not automatically treated as a negative biological result.
- Candidate evidence table: **2,426** rows = 2,425 gene-region pairs + one anchor-only row.
- **43** old BBJ/EAS locus-and-stable-Ensembl-ID-matched coloc rows; **0** old coloc probabilities transferred to EUR regions by a shared gene name.
- **2,382** positional gene-region pairs lacking matched prior eQTL coloc (NOT_TESTED).
- All **2,425** positional pairs still require proper sQTL and pQTL screening at expansion scope.

## Important scientific guardrails

These are *candidate universe* counts, not 2,225 mechanistically established stroke genes. Distance clusters are not LD-independent loci. A genome build GRCh37 allele-pair match and gene-window overlap is not colocalization. EUR genotype LD must be built separately; reusing 1000G EAS for EUR would be invalid for primary fine-mapping. GIGASTROKE studies share participants across meta-analysis and subtypes; replication claims need disjoint cohorts or explicit overlap treatment.

A gene lacking assayed QTL or stroke-cell evidence is `NOT_TESTED`, not `NEGATIVE`. Normal cell-state / mouse model expression is not a substitute for human disease mechanistic evidence.

## Reproduce

```bash
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
python3 scripts/is/run_ancestry_expansion_v2.py
python3 -m unittest discover -s tests -p 'test_is_*' -v
```

Real data:
`/srv/is-analysis/data/is/reference/gigastroke/additional/GCST90104540_buildGRCh37.tsv.gz`
`/srv/is-analysis/data/is/processed/gigastroke/broad_v1/GCST90104540_AIS_EUR_GRCh37.allele_pairs.tsv.gz`.

Derived artifacts:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v2_ancestry/`
including `IS_EUR_AIS_PROVISIONAL_REGIONS.tsv`, `IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv`, `IS_ALL_GENE_WINDOW_UNIVERSE.tsv`, `IS_ANCESTRY_GENE_EVIDENCE_V3.tsv`, and SHA-256 reproducibility manifest.

## Next priority

1. Genotype reference for EUR GWS loci, proper normalization and EAF/QC. Do not copy EAS fine-mapping LD into EUR.
2. Proper signal-level SuSiE for G0022 chr12 (AS 2 and AIS 3 1KG EAS clump indices only suggest possible multi-signal locus, not credible sets).
3. Genome-wide molecular QTL/coloc and vascular/brain human stroke cell-state coverage for all 2,425 gene-region associations.
4. Additional ancestry-specific GWAS only when justified by samples/phenotype and full source provenance. Separate independent replication from overlapping meta cohorts.
