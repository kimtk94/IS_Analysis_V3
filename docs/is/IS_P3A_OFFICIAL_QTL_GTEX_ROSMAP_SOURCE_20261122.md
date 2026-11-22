# IS P3A — Official QTL Catalogue source discovery and genome build gate (2026-11-22)

**Scope:** MASTER DEGREE / ischemic stroke P3A, 36 actual gene×QTL groups, 20 genes, one chr10 locus. This is source provenance and download readiness only, **not new Coloc evidence**.

## Genuine official source evidence

Official eQTL Catalogue reference [eQTL Catalogue release 7 dataset metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/07918a267afdc1c84b5444bea9bd902e7201728f/data_tables/dataset_metadata_r7.tsv) and [Tabix paths](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/07918a267afdc1c84b5444bea9bd902e7201728f/tabix/tabix_ftp_paths.tsv) were both downloaded, Git blob SHA1 verified against **immutable commit** \`07918a267afdc1c84b5444bea9bd902e7201728f\`, and cross-checked against the real [36-group manifest](evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json).

| QTL dataset | Study | Tissue (official) | Catalogue sample N | P3A assay groups | Published paper PMID | Verified uncompressed schema |
| --- | --- | --- | ---: | ---: | --- | --- |
| \`QTD000131\` | GTEx (\`QTS000015\`) | artery (aorta) | **387** | **10** | 32913098 | \`ge\`, \`.all.tsv.gz\` |
| \`QTD000171\` | GTEx (\`QTS000015\`) | brain (cortex) | **205** | **11** | 32913098 | \`ge\`, \`.all.tsv.gz\` |
| \`QTD000434\` | ROSMAP (\`QTS000025\`) | brain (DLPFC) | **560** | **15** | 28869584 | \`ge\`, \`.all.tsv.gz\` |

The sample N values are explicitly **dataset-level metadata**. Do not fill per-gene effective N for an ABF model until source-specific sample and test metadata are audited. The QTL Catalogue file is a nominal gene-expression QTL summary, not the significant-only BH05 eGene table. Historical 25/11 server source-stat memberships remain individually unresolved.

## Original remote files (HEAD checked in actual GitHub Actions, HTTP 200)

| Dataset | Source compressed bytes | Tabix .tbi bytes | Official FTP path |
| --- | ---: | ---: | --- |
| QTD000131 | 3,473,236,781 | 2,174,813 | \`ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000131/QTD000131.all.tsv.gz\` |
| QTD000171 | 3,496,087,269 | 2,247,876 | \`ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000171/QTD000171.all.tsv.gz\` |
| QTD000434 | 4,039,836,192 | 2,139,304 | \`ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz\` |

Total source compressed size (all three): **11,009,160,242 bytes**. The HTTP HEAD check establishes remote existence and advertised size **only**: it does **not** download compressed full source or tabix indexes, or independently verify variant row count or file digest. This CI avoids full transfers of ~11 GB from EBI.

Original QTL Catalogue [column documentation](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/07918a267afdc1c84b5444bea9bd902e7201728f/tabix/Columns.md) establishes **GRCh38**, \`ref\`, \`alt\` as the effect allele, \`beta\`, observed \`se\`, \`maf\`, \`an\`, \`ac\`, chromosome/position, and molecular trait ID. It may repeat exactly the same variant-effect rows under multiple rsIDs; deduplicate by \`chromosome:position:ref:alt:molecular_trait_id\`, not rsID alone.

## Critical build and scientific HOLD

The legitimate current P3A GIGASTROKE original GWAS \`GCST90104540\` (EUR) and \`GCST90104545\` (EAS) files were previously examined as **GRCh37**. The official nominal molecular QTL datasets are **GRCh38**.

**Never join their \`chromosome:position:ref:alt\` keys before a verified build transformation.** Neither the legacy candidate region \`chr10:104663676–105720351\` in GRCh37 nor the palindromic lead \`10:105596641:G:C\` can be directly reused as a GRCh38 Tabix region without a checked coordinate conversion.

Required remaining evidence:
1. Pin official \`hg19ToHg38\` UCSC chain checksum, perform chr10 region boundary and anchor SNP liftOver, reject unmapped/multimapped coordinates, recheck REF/ALT against versioned GRCh38 reference FASTA. Reverse strand and palindromic alleles require explicit adjudication.
2. Only then use **one carefully scoped, rate-limited** indexed Tabix query per QTL dataset on GRCh38, rather than downloading all 11 GB. EBI explicitly warns that frequent incomplete Tabix requests can trigger anti-abuse protections; avoid aggressive retries.
3. For all 36 group keys, validate actual nominal \`gene_id\`, beta, observed SE, MAF/alleles and full cis coverage. Retain null/missing source groups and excluded duplicate/ambiguous variants in the QC manifest.
4. Obtain the **exact** GWAS phenotype/source-specific Ncases/Ncontrols/fraction, and verify original GWAS per-SNP beta/SE REF/ALT and full cis coverage. GTEx/ROSMAP N values do not fill this GWAS component.
5. Only when matched per-variant inputs, allele direction, sample ascertainment and predeclared full cis set are individually verified may \`coloc.abf\` be considered. The standard ABF does not require LD; \`coloc.susie\` requires signed ancestry-matched LD and adequate fine-mapping QC.

**Coloc.abf eligible 0/36; executed 0/36. SuSiE eligible 0/36; executed 0/36. No biological causality or independent replication established.**

## Reproducible GitHub Actions evidence

- [Source recovery, sample N and build preflight CI #45024290928](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/45024290928): SUCCESS, independent metadata + remote source/index HEAD. Synthetic 11/11 tests; repository-wide regression 244 tests executed, 243 passed and 1 skipped.
- \`scripts/is/audit_official_p3a_qtl_catalogue_provenance.py\` validates immutable metadata, cross-checks all three datasets, and emits exact 36 rows with \`HOLD\`.
- \`tests/test_is_official_p3a_qtl_catalogue_provenance.py\` rejects sample, tissue, study, build and original source drift.
- \`.github/workflows/is-p3a-official-qtl-source.yml\` runs everything on GitHub Actions. Artifacts contain only metadata, URLs, file size and explicit 36-group blockers; never full raw genomic data or patient-level information.

**Scientific next:** Build the official cross-build locus conversion and reference FASTA allele check, then fetch at most three appropriately scoped indexed QTL regions with complete source provenance and merge the resulting 36-key evidence matrix.
