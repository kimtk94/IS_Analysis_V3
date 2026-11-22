# IS P3A — verified indexed original QTL Catalogue chr10 context for all 36 assays

**Checked 2026-11-22, MASTER DEGREE / ischemic stroke P3A only.** This is a verified *original context source survey*, not evidence of gene-wise full-cis completeness, GRCh37 GWAS matching, causal inference, or any approved formal Coloc calculation.

## Actual CI evidence

- Official UCSC hg19 → hg38 chr10 6-anchor reference and FASTA run: [GitHub Actions #45024963210](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45024963210), **SUCCESS**.
- Three **actual** eQTL Catalogue GRCh38 source region fetches: [#45025825454](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45025825454), **SUCCESS**; corrected-parser and further 13-test branch validation: [#45025844092](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45025844092), **SUCCESS**.
- New indexed-source tests 13/13 PASS and repository regression 269 executed / 268 PASS / 1 SKIP on the final branch source CI. The skipped test requires unavailable Rscript and is unrelated to the new QTL source.
- Source release-7 identities fixed by immutable eQTL Catalogue Git commit `07918a267afdc1c84b5444bea9bd902e7201728f`, with publisher URLs/head checks validated in [PR #308](https://github.com/kimtk94/IS_Analysis_V3/pull/308).
- Original archive study/assay IDs fixed by immutable previous scientific source commit `0bf646fd6f29322d17da878253f37def24c01c74`.
- Source queries: **one 2,056,675-bp outer context** `GRCh38 chr10:102403919–104460593` per tissue QTL; the genomic scope is based on six independently FASTA-checked GRCh37/hg38 reference anchors. The full genomic region was NOT lifted SNP by SNP. HTSlib/PySam indexed remote reader, sequential read and throttle; 3.47/3.50/4.04-GB whole original sources were not downloaded.
- Only count-level QC metadata, URLs, reference checksums and no raw per-SNP association rows are uploaded to Actions.

## Authentic full source region statistics (NOT distinct or independent SNPs)

| Source QTL | Tissue / catalogue donor N | Groups | Original 2.05Mb region gene×variant rows | Original 36-candidate source rows | QC-valid unique SNP×gene records | QC-valid **nonpalindromic** SNP×gene records | Rejected from single-letter/strict variant key validation | Palindromic or other non-SNP among key-valid |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `QTD000131` | GTEx artery aorta, catalogue N=387 | 10 | **181,828** | **46,737** | **41,779** | **36,141** | **4,958** | **5,638** |
| `QTD000171` | GTEx brain cortex, catalogue N=205 | 11 | **196,820** | **42,647** | **38,089** | **33,029** | **4,558** | **5,060** |
| `QTD000434` | ROSMAP brain DLPFC, catalogue N=560 | 15 | **267,103** | **74,501** | **49,282** | **42,008** | **25,219** | **7,274** |
| **Total** | Three tissue cohorts, not three independent discoveries | **36** | **645,751** | **163,885** | **129,150** | **111,178** | **34,735** | **17,972** |

For all three sources, among records meeting the single-nucleotide canonical `chr10_position_REF_ALT` key, tested REF vs official GRCh38 UCSC FASTA yielded **0 mismatches**, numeric signed `beta`, observed positive `se`, MAF and p validity **0 errors**, and count of conflicting duplicate variant/rsID entries **0**. The 34,735 failed strict single-letter/variant-key screen must be examined separately for legitimate multiallelic or insertion/deletion formatting before any inference; they are not 34,735 proven erroneous genotypes. A group may share a SNP with other tissues/genes, so **111,178 is not a unique genomic SNP count**. All 36 of the original QTL gene×tissue assays had nonzero valid SNP association records in the queried context; **none is independently approved for Coloc**.

Exact SHA-256 of streamed QTL regional TSV:
- QTD000131: `634d468798831fc0aae6b5a409dec495b948d3c84a80b8b1cbd2977e0e618603`
- QTD000171: `97ca7865a70257ccab0e5680e7c19f52cd870d5655846599409ce8174ab89edc`
- QTD000434: `5e7091e60937a0f0f27fa8490982ec548e00dbba0fb6d11a11581ab9cc3f327c`

These regional stream hashes are **checksums of returned indexed query text**, not publisher-verified digest of the complete `.all.tsv.gz` files and not evidence of complete gene cis overlap.

## Clinical and statistical decisions

- **Technical source consistency:** PASS for 3 actual source indexed chr10 contexts, candidate assays 36/36 with observed QTL rows.
- **Verified liftOver reference anchors:** PASS for 6 UCSC hg19/hg38 coordinates including old palindromic `10:105596641:G:C` → hg38 `10:103836883` with REF G on both. The original lead remains C/G palindromic; do not assert independently resolved signed effect orientation.
- **All 36 gene-specific full cis windows:** NOT independently recovered. A shared chr10 ~2Mb positional context is a bounded survey, not nominal full tested cis region per gene.
- **Per-SNP GRCh38 QTL ↔ GRCh37 original GWAS harmonization:** NOT executed. REF/ALT effect, multiallelic, strand, EAF and allele orientation still require per-SNP bidirectional build validation and original GWAS matching.
- **QTL effective gene/tissue sample N:** NOT approved. Release-7 catalogue donor N is source-wide metadata, not per-test observed effective sample size.
- **GIGASTROKE exact phenotype/Ncases/Ncontrols and BBJ inclusion:** NOT source-approved, and GIGASTROKE EAS is not established as independent of BBJ.
- **Formal `coloc.abf` eligible / executed:** **0 / 36**, **0** calculations. The standard `coloc.abf` does not require LD itself, but requires valid matched summary data and model parameters.
- **Formal `coloc.susie` eligible / executed:** **0 / 36**, **0**; later requires ancestry-matched, allele-signed, order-consistent LD and fine-mapping QC.

## Repo changes and follow-up

- `scripts/is/audit_p3a_chr10_ucsc_liftover_fasta.py` / `tests/test_is_p3a_chr10_ucsc_liftover_fasta.py` / `.github/workflows/is-p3a-ucsc-liftover-fasta.yml` / [UCSC evidence](IS_P3A_OFFICIAL_CHR10_LIFTOVER_FASTA_20261122.md).
- `scripts/is/audit_p3a_indexed_qtl_real_chr10.py` / `tests/test_is_p3a_indexed_qtl_real_chr10.py` / `.github/workflows/is-p3a-official-indexed-qtl-chr10.yml`.
- Only official reference and aggregated QTL source QC artifacts are uploaded; no private genotypes, patient data or complete original genomic source.

**Next science gate:** preserve the 36 key registry; resolve 34,735 excluded variant records by type, independently establish QTL source-specific *true full cis windows* and observed donor N, then perform per-SNP GRCh37→38 source/liftover and GWAS allele-SE matching with exact stroke GWAS Ncase/Ncontrol. Only source-verified eligible groups may enter `coloc.abf`, followed by prior sensitivity and assumptions/limitations. All current groups remain **HOLD**.
