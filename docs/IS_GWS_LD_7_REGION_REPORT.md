# IS expanded GWS LD reference — verified results (2026-10-09)

## Research scope

We retain all 30 coordinate-overlap discovery components across BBJ and GIGASTROKE EAS, including 7 components with a GWS association and 23 suggestive-only components.

**Regional genotype reference acquisition is now complete for all seven GWS components**, not all 30 discovery components.

## Actual 1KG EAS genotype extraction

Three new GWS components, 504 EAS samples, indexed GRCh37 1000 Genomes phase3 VCF:
- IS_XDATA_G0004 chr2: 28,343 biallelic SNPs; lead 2:26929282:G:T; reference allele ordering reversed.
- IS_XDATA_G0008 chr4: 29,983 biallelic SNPs; lead 4:111704199:G:T; reference ordering exact.
- IS_XDATA_G0021 chr12: 28,375 biallelic SNPs; lead 12:90019178:G:C; reference ordering reversed.

The four pre-existing BBJ anchor regions remain separate original genotype sources: G0007/L001, G0015/L002, G0022/L003, G0023/L004.

**Corrected ID trap:** 1000G regional VCF has ID='.'; the initial PLINK2 PGEN conversion propagated missing IDs, and a first r2 calculation with lead='.' was therefore INVALID. It was corrected by re-importing VCF to new `*.stable.pgen/.pvar/.psam` with `--set-all-var-ids '@:#:$r:$a'`. Never use the unnamed legacy PGENs or their old `G0004_DOSAGE.traw` artifact for LD analysis.

## Verified corrected new GWS lead-centered reference r²

| Group | r² ≥ 0.2 proxies | r² ≥ 0.8 proxies | maximum proxy r² |
| --- | ---: | ---: | ---: |
| G0004 chr2 | 101 | 15 | 0.98934269 |
| G0008 chr4 | 194 | 37 | 1.0 |
| G0021 chr12 | 53 | 3 | 0.98810649 |
| **Three new GWS total** | **348** | **55** | — |

This is reference-genotype correlation, not an independent locus test, fine-mapping credible set or causal variant identification. Selection used autosomal ±500kb lead windows for reference extraction, with source GWAS interval details recorded in the execution manifest.

## Source provenance and data lifecycle

New reference files:
`/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas/IS_XDATA_G00*.1KG_EAS.GRCh37.vcf.gz`
and
`/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas/IS_XDATA_G00*.stable.{pgen,pvar,psam}`.

Outputs:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/new_gws_ld_proxy/`
and `IS_BROAD_REFERENCE_READINESS.tsv` plus its summary JSON.

The 23 suggestive-only components are still missing local regional EAS genotype references. New molecular QTL colocalization, signal independence, reference-standardized effect harmonization, and fine-mapping have **not** been completed.

## Next scientific steps

1. Prioritize EUR AIS GIGASTROKE GCST90104540 using source-verified MD5; run `acquire_gigastroke_verified.py --accession GCST90104540 --execute` deliberately. Do not equate it with independent EAS replication.
2. Normalize with `normalize_gigastroke_verified.py`, preserving effect/other alleles. The GWAS raw file has no REF/ALT, so use genotype/reference checks before multi-ancestry meta-analysis.
3. For seven GWS components, harmonize GWS, reference and gene-QTL alleles, handle allele flips/palindromes, then fine-map independently.
4. Expand the LD extraction to 23 suggestive-only components as budget and scientific priority allow; never label them genome-wide significant.
5. The 869-gene positional universe is inclusive but not yet fully regulatory/distal-gene mapped.
