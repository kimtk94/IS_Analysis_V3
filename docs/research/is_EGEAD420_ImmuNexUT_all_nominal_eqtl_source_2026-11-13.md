# GIGASTROKE IS regulatory follow-up: genuine East Asian full nominal cis-eQTL source (2026-11-13)

**Internal research provenance note / execution gate; no inference of stroke causality.**

## Decisive correction: the earlier 105-donor Japanese NHA000051 is significance filtered

Original public NHA000051 `hum0099.v1.eqtl.v1.zip` was verified as an actual publicly accessible ZIP with `eQTL_gene_level.tar.gz`, `eQTL_exon_level.tar.gz` and original methods DOCX. The real nested `eQTL_gene_level/README.txt`, read and verified in [GitHub Actions #43428442728](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43428442728), explicitly states:

- `raw_data/[CD4,CD8,B,NK,Mono,PB]/chr[1:22]_cis_eqtl_mapping.txt.gz` from MatrixEQTL, ±1 Mb each gene, **filtered to p ≤ 0.05**, header `SNP, gene, beta, t-stat, p-value`, beta alternative allele effect on normalized expression.
- `permutation/[CD4,CD8,B,NK,Mono,PB]_permutation_q_value_0.5.txt.gz`: **eGene permutation Q ≤ 0.5**, not all tested pairs. 105 healthy Japanese donors.
- Therefore this public gene-level archive is **NOT suitable for unfiltered region-wide coloc.abf or SuSiE-coloc by naïvely treating available SNP-gene pairs as all cis tested variants**. Its 1.8GB physical size and the average number of reported SNP/gene in later methods cannot overrule the actual explicit filtering statement. Keep candidate/power or replication *sensitivity* role only.

## Newly identified primary EAS eQTL candidate: ImmuNexUT E-GEAD-420

**Authoritative sources**:
- [Official NBDC E-GEAD-420 dataset](https://humandbs-production.ddbj.nig.ac.jp/en/dataset/E-GEAD-420): **Unrestricted-access nominal QTL, including nonsignificant associations**; n=416 overall, 79 healthy and 337 cases across 10 immune diseases. Reference genome GRCh38, QTLtools nominal pass within **TSS ±1 Mb**. The deposited study states ALT allele corresponds to effect-size slope; date 2021-06-09, study hum0214.
- [ImmuNexUT original FAQ](https://www.immunexut.org/faqs) explicitly distinguishes E-GEAD-398 (significant conditional eQTL) from **E-GEAD-420 nominal INCLUDING nonsignificant**, 28 immune-cell contexts. Browser FDR hits/top1,000 per gene are NOT representative of the raw all-nominal archive.
- [Publisher DDBJ directory](https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/) and [official filelist](https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/E-GEAD-420.filelist.txt): `E-GEAD-420.processed.zip` **43,507,433,191 bytes**, published MD5 `0231fc7ca80b4a604f31661d72ae62d8`, containing `nominal_eQTL_2.tar.gz` **43,495,814,577 original member bytes** with independent published MD5 `41544fb671177aeee7fd2d0634b90ae2`. Published MD5s are authoritative **manifest values, NOT independently validated whole-file checksum unless all 43.5GB downloaded**. `E-GEAD-420.idf.txt` official MD5 `0672e92775fd103fd4a848d653de1f1a`.

**Official study design**:
- Japanese immunologic cohort (predominantly EAS, but every eQTL subgroup ancestry still needs verification).
- 416 participants overall: 79 healthy, plus individuals with systemic lupus, myositis, systemic sclerosis, mixed connective tissue disease, Sjögren disease, rheumatoid arthritis, Behçet disease, adult-onset Still disease, ANCA vasculitis, Takayasu arteritis. **Do NOT treat 416 as independent healthy participants, or assume 416 analyzable in every cell type**.
- 28 immune cell contexts, QTLtools normalized expression with PEER factors; top 2 genotype PCs, sample collection phase, clinical diagnosis and sex as covariates. Tested variants within **±1 Mb TSS**; imputed GRCh38 WGS / REF ALT effects.
- Need participant numbers per cell context, matched disease groups and sample stability. This is peripheral immune regulatory data, not ischemic stroke brain eQTL or stroke-case causal response.

**Current executable probe**:
- `scripts/is/audit_EGEAD420_nominal_zip64.py`, `tests/is/test_EGEAD420_nominal_zip64.py`, `.github/workflows/is_EGEAD420_nominal_source_qc.yml`, [GitHub CI #43430780352](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43430780352).
- Downloads only tiny original DDBJ filelist/IDF and last few hundred KB ZIP64 footer and member's first ~2MB raw bytes. Validates small IDF MD5, ZIP64 directory/member metadata and nested nominal TAR first headers. Full publisher ZIP MD5 **NOT** verified unless entire file downloaded. At time of this note, this workflow was **QUEUED**. No claims of a PASS.
- For pilot locus-level use, next assess ZIP member compression method, whether the nominal TAR can be efficiently read/ranged by cell context and chromosome without downloading 43GB; if not, use a GitHub Actions cache/object store or bounded per-locus query source, never repeatedly download entire 43GB.

## EUR comparator and GWAS source controls

- EBI BLUEPRINT `QTD000021.all.tsv.gz` [actual source row QC #43428082712](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43428082712) **PASS**: original **300 nominal association rows** validated, **287 p>0.05**, beta/SE/MAF/GRCh38 REF/ALT and Ensembl gene IDs confirmed. Entire QTL variant universe and BLUEPRINT donor ancestry not established by a small prefix.
- Public GTEx v8 EUR cis-eQTL archive actual TAR source bytes [#43427462917](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43427462917) PASS; it may include only significant eQTLs, so do not use as a full-cis comparator unless member-level full nominal associations are verified. Official GTEx complete allpairs files can be large/requester pays.
- Original GIGASTROKE IS EUR GCST90104540 and EAS GCST90104545 whole-source, matching ancestry 1000G Phase3 BIM **22/22 chromosome coordinate and SNP pair QC PASS**, cross-evidence [#43033560457](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43033560457). **No beta-effect direction alignment, EAF resolution of palindromic SNPs, chromosome build liftover, population LD R matrix or colocalization yet**.
- GIGASTROKE GWAS is GRCh37, ImmuNexUT E-GEAD-420 is GRCh38, so use validated liftOver (including alternate allele/ref concordance), not position-only joins. Compare only identical variant sets after QC and report coverage, power differences and tissue match.
- No MAGMA, LDSC-SEG, SuSiE, coloc.abf or MR outputs have been generated by this research-source preflight.

## Explicit go/no-go gates

1. **P0 data existence**: ImmuNexUT official all-nominal association listing **DOCUMENTED**, archive public. Real raw member-byte/readability still pending.
2. **P1 source integrity**: verify official manifests for source; full 43.5GB archive MD5 **NOT_TESTED**. Tiny source ZIP64 ranges can verify directory format, but not full file contents.
3. **P2 row schema and allele provenance**: inspect actual nominal per-variant records from a genuine relevant cell type (e.g. CD14 monocytes), requiring GRCh38 variant, REF/ALT, beta, se or t-stat plus appropriate residual degrees of freedom, allele frequencies, sample N. **HOLD** until real association rows.
4. **P3 complete cis region**: check that nonsignificant variants exist per selected gene and cross-check expected tested variants / local 1000G after conversion; never infer from the data publisher's title alone. **HOLD**.
5. **P4 final colocalization**: only when harmonized GWAS/EAS eQTL + EUR comparator full nominal, LD matrix and sample size uncertainties are documented; model multi-signal conditional coloc if feasible. **BLOCKED**.
6. **Research originality**: 2026 MAEEA already establishes utility of ancestry-matched eQTL for multiple phenotypes; novelty must be **specific IS locus-to-gene transitions** robust to ancestry LD/MAF, cell types and case-control power, not the generic use of EAS eQTL.

**Decision:** Promote publicly documented E-GEAD-420 (not NHA000051) to **EAS full-nominal source pilot P0**, with genome-wide colocalization still BLOCKED.
