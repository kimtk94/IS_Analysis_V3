# IS P3A — official UCSC hg19 (GRCh37) → hg38 (GRCh38) chr10 reference audit

**Date:** 2026-11-22. **Study:** MASTER DEGREE / ischemic stroke, not Sales. **Scientific status:** `HOLD`; formal `coloc.abf` eligible/executed **0/36**, `coloc.susie` eligible/executed **0/36**.

## Real original reference execution

[GitHub Actions #45024963210](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45024963210) completed **SUCCESS** for real official source downloads, MD5 verification, six-point liftOver and real reference FASTA checks. Synthetic negative tests 12/12 PASS; repository-wide regression **256 executed, 255 PASS, one skip (Rscript absent)**.

| Original reference | Original compressed bytes | UCSC published MD5 | Observed SHA-256 |
|---|---:|---|---|
| UCSC `hg19ToHg38.over.chain.gz` | 227,698 | `35887f73fe5e2231656504d1f6430900` | `5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d` |
| `hg19/chromosomes/chr10.fa.gz` | 43,026,389 | `508a1e86cc0f64257dd0c9775462eebf` | `899b53d83b2e880a1e37433fa2001fc14a93e79fb06a80719c3eecca6960a5c7` |
| `hg38/chromosomes/chr10.fa.gz` | 43,157,332 | `b7d478cb06a84400f74a8932cd5fad0d` | `4425ab3327d51791e1c2361ffb283bb732ca77cfa1d76b9118d0f8eae1a568c8` |

UCSC checksums are from the original public [hg19 liftOver MD5](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/md5sum.txt), [hg19 chr FASTA MD5](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/md5sum.txt) and [hg38 chr FASTA MD5](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/md5sum.txt).

### Exact predeclared original and converted chr10 coordinates

| Region / anchor | GRCh37 position | GRCh38 position | FASTA REF | Lifted strand |
|---|---:|---:|---|---|
| Context start | 104,163,676 | 102,403,919 | C → C | + |
| Original selected locus start | 104,663,676 | 102,903,919 | A → A | + |
| Original selected locus midpoint | 105,192,013 | 103,432,256 | C → C | + |
| Historic EAS GWAS G/C lead | 105,596,641 | 103,836,883 | **G → G** | + |
| Original selected locus end | 105,720,351 | 103,960,593 | REF matched | + |
| Context end | 106,220,351 | 104,460,593 | REF matched | + |

**Outer candidate retrieval window:** GRCh38 `chr10:102403919–104460593`. **Inner predeclared locus:** `chr10:102903919–103960593`. These are six **reference-anchor** mappings, not confirmation of every base/indel/individual SNP or a proven full gene-specific cis interval. The archived GWAS lead is **C/G palindromic**, so aligned reference G does **not** remove strand ambiguity in signed association comparisons by itself.

## Operational safeguards

- `scripts/is/audit_p3a_chr10_ucsc_liftover_fasta.py` downloads only public reference chr10 FASTAs and chain within an ephemeral runner directory, verifies source MD5/SHA, maps 1-based original points correctly from 0-based chain coordinates, and checks true FASTA REF bases.
- `tests/test_is_p3a_chr10_ucsc_liftover_fasta.py` rejects duplicate position, unmapped anchor, wrong chromosome, wrong strand, mismatch REF, wrong candidate groups and false `coloc` approval.
- `.github/workflows/is-p3a-ucsc-liftover-fasta.yml` uses a pinned official conversion implementation and uploads only source hashes, mapped anchors and 36-group HOLD metadata.
- Any subsequent QTL Tabix query is strictly a **shared 2.05Mb contextual survey**, not an exhaustive groupwise full cis QC. Only an original per-gene tested cis window, variant-wise liftover and GRCh37 official GWAS harmonization can progress to ABF acceptance.

## Required next science

1. One rate-limited region query from each official QTL Catalogue GRCh38 `QTD000131/171/434` BGZF+Tabix file. Verify actual `gene_id`, `variant` key, REF vs official UCSC hg38 FASTA, ALT effect allele, measured beta/SE, MAF, duplicate rsIDs and study-wide sample metadata. Individual per-gene effective N must still be proven.
2. Map all needed QTL/GWAS SNP coordinates individually across builds, with multiallelic/indel, duplicate, unmapped, reverse-complement and palindromic exclusions. Determine full true cis window for each gene/assay; a shared context scan is not a substitute.
3. Independently approve the exact ischemic stroke GWAS phenotype, EAS/AIS cases and controls, effective sample size approach and BBJ overlap, then analyze only scientifically approved groups using prior sensitivity. Do **not** compute causal H4 from historic SNP-only flags.

**Current scientific verdict: Coloc approval 0/36; all source and cohort N gates remain HOLD.**
