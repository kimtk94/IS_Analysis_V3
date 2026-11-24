# IS P3A — BBJ 2021 public phenotype source cross-check (2026-11-02)

## Empirical CI execution

Run [#41430972588](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41430972588) completed SUCCESS. Seven offline regression tests passed. The bounded live network check found:

- **PheWeb.jp official BBJ phenotype page**: all 4 checks matched ("Ischemic stroke", 22,664 cases, 152,022 controls, Sakaue 2021 citation); source `https://pheweb.jp/pheno/IS`.
- **IEU OpenGWAS secondary cross-check**: HOLD. The legacy gwas.mrcieu.ac.uk URL redirects to the new JavaScript-rendered opengwas.io page. The non-authenticated HTML did not provide machine-readable dataset metadata in the Actions job; **the source has not contradicted the study metadata**.
- **Overall metadata crosswalk**: HOLD, because the two-source predicate was not satisfied. Successful CI execution is **not** scientific PASS.
- **EAS5 original compressed file SHA-256**: still missing from the old V1 comparison output. No way to prove unchanged source bytes from log arithmetic alone.

A successor-site redirect-aware parser and regression for `CLIENT_RENDERED_PAGE_NO_ACCESSION_METADATA` were added after the initial CI run. Check the later workflow logs before upgrading the parser verdict.

## Additional public source evidence (not data provenance authorization)

- Original official published [PheWeb phenotype](https://pheweb.jp/pheno/IS): 22,664 cases, 152,022 controls, 174,686 participants; Sakaue & Kanai et al. 2021.
- Public [IEU OpenGWAS dataset index](https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90018644/) previously indexed GCST90018644, East Asian, HG19/GRCh37, 22,664/152,022. The **current site is a JavaScript application** and metadata API calls typically require authentication, so the new CI does not scrape protected endpoints or request tokens.
- [OpenGWAS API documentation](https://api.opengwas.io/api/) says all endpoints except /status and /batches generally require JWT. This is not a valid reason to bypass authentication.
- [BIGA's public BBJ reprocessing notes](https://www.bigagwas.org/documentation) describe GRCh37 coordinates and Allele2 as the effect allele. This is **secondary harmonization documentation**, not the exact original BBJ source manifest or an independent FASTA/REF check.

## Research gate

A historical [25,441-variant BBJ/EAS5 comparison](BBJ2021_EAS5_SOURCE_COMPARISON_EVIDENCE_20261102.md) strongly suggests the archive was derived from Sakaue 2021. However, none of the above proves **per-variant REF/ALT harmonization**, source EAS5 SHA-256, a specific matched GWAS SE field, or release-level independent replication.

**Scientific status: HOLD**. Do not adopt BBJ N into coloc.abf or claim BBJ and GIGASTROKE EAS as independent cohorts. Next actual source-dependent gate is re-read original EAS5 .tsv.gz bytes from its authorized storage, record SHA256 and exact SNP/allele source mapping; leave that gate HOLD while storage is inaccessible. The CI-only public metadata work needs no SSH.
