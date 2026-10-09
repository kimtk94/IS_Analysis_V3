# Ischemic Stroke — Broad discovery evidence gates (2026-10-09)

## Confirmed analysis
- 6 EAS/Japanese GWAS files were scanned.
- 45 phenotype-specific proximity clusters, collapsed into 30 coordinate-overlap components.
- 15/45 source regions contain genome-wide significant variants; 7/30 merged components contain ≥1 such signal.
- Only 4/30 lead variant coordinates occur in existing BBJ-focused EAS regional PVARs.
- REF/ALT audit: 3 exact, 1 swapped, 26 unavailable.
- **No component has been declared statistically independent; no new gene has been claimed causal.**

## Added execution gates
1. `build_broad_discovery_v1.py`
2. `audit_broad_regions.py`
3. `merge_broad_intervals.py`
4. `audit_broad_ld_coverage.py`
5. `build_expanded_region_manifest.py`
6. `audit_lead_alleles.py`
7. `map_all_genes_gencode_v19.py` (needs validated GRCh37 GENCODE v19 GTF)

## Annotation integrity
Source: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
Do not read incomplete .part downloads as a reference. Check gzip integrity and record SHA256 before mapping.

Output `IS_ALL_GENE_WINDOW_UNIVERSE.tsv` includes protein-coding AND noncoding genes within ±500 kb windows. This is a **window-based gene universe**, not proof of causal involvement. Include distal gene links from molecular QTL, chromatin accessibility, enhancer-to-gene and chromatin interaction before final ranking.

## Next critical source gaps
- Broad EAS ancestry-matched phased reference for non-anchor loci (26/30 lead positions missing).
- All-ancestry GIGASTROKE and other GWAS for higher-powered discovery; EAS-only subset is not full GIGASTROKE.
- Genome-wide molecular QTL datasets, subtype-specific fine-mapping/colocalization.
- Preserve discovery vs independent replication cohorts and overlap metadata.

## Guardrails
- No new causal-gene claim based on proximity alone.
- Do not interpret no-expression as negative disease evidence.
- Subtype datasets may share samples, so multi-dataset overlap ≠ independent replication.
- No server cron or automatic massive VCF transfer unless explicitly budgeted and source license checked.
