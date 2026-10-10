# IS priority-four EUR LD source plan — 2026-10-11

**Scope:** FGF5, C4orf22, CALHM2, NEURL. Do not narrow the parent 2,225-gene universe. The 4 targets represent two original BBJ IS loci, not four independent association regions.

| Locus | Genes | Existing matched source GRCh37 SNP coordinate range | Proposed GRCh37 regional extract (+100kb) | GRCh38 original QTL position range |
|---|---|---|---|---|
| BBJ_IS_L001, chr4 | FGF5, C4orf22 | 80,681,087–81,684,150 | **chr4:80,581,087–81,784,150** | 79,759,933–80,762,996 |
| BBJ_IS_L002, chr10 | CALHM2, NEURL | 104,839,152–105,959,690 | **chr10:104,739,152–106,059,690** | 103,079,395–104,199,932 |

Source: the original GTEx/BBJ exact allelic SNP input from 10 numerically replayed gene–tissue assays. Each coordinate range is the union of all relevant input SNPs. UCSC hg19ToHg38 direction was separately regression-checked against rs671, rs11066015 and chr10:104839152->103079395. Chain ref and alt labels are not a substitute for independent GRCh38 reference-FASTA allele validation.

The server has 54 previously generated 1000G EUR region panels from an unrelated metabolic-resilience project. Their PVAR reference is b37 (hs37d5) and none overlaps any of the 30 original IS assays after attempted genome-build alignment; original EUR panels cover chromosomes 1,2,4,5,7,9,11,15,16,17,20,22, whereas the 30 IS assays cover chromosomes 4,10,12,13. Region-specific scope, not lack of EUR public data, is the primary limitation.

**Public prospective source:** 1000 Genomes Phase3 20130502 GRCh37 integrated genotype VCF indexed by chromosome, official FTP release at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ . Requires source endpoint and tabix remote-index verification; a direct HTTPS HEAD probe of chr4 index returned HTTP 404 during this audit, so no actual download or index capability is asserted. See also official project reference https://www.internationalgenome.org/announcements/ .

Before future inference: obtain small region-specific EUR samples from an authorized, checksum-verified source; confirm exactly 503 EUR individuals using sample metadata; harmonize b37->b38 variant identity and REF/ALT with GRCh38 FASTA; compare SNP intersection, MAF and missingness; compute signed LD and its allele convention; retain genotype ancestry-sensitive distinction. **GTEx donor-cohort LD is still not verified**, therefore a 1000G EUR LD matrix is only a sensitivity proxy, not an in-study QTL substitute or publication-ready coloc.SuSiE validation.

Created read-only existing-source regional extraction manifest: `/srv/is-analysis/results/is/audits/is_30_eur_region_acquisition_manifest_20261011_v1/IS_PRIORITY4_EUR_GRCH37_REGIONS.json`. No restricted data, production runtime, or original QTL tables modified.
