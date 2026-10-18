# IS P3A official GWAS source identity – acceptance gate v1 (2026-10-18)

## Verified baseline

- [Server acceptance #189](https://github.com/kimtk94/IS_Analysis_V3/issues/189): full-cis archived 24-column TSV has QTL beta/SE and BBJ/EAS beta/p, but has **no BBJ/EAS N, Ncases/Ncontrols, sample case fraction or directly observed GWAS SE fields**.
- 36 ABF review groups; 25 groups with source FASTA flags and numeric summary-stat observations; 11 without. 68,612 source-qualified rows, not independently harmonized unique GWAS SNPs.
- Current genuine coloc.abf execution count: **0**.

## Candidate GWAS sources — not automatically matched to original EAS5

| Candidate | Candidate accession | Public study metadata | Status |
|---|---|---|---|
| BBJ ischemic stroke | `bbj-a-129` | ncase=17,671, ncontrol=192,383; HG19/GRCh37, East Asian | CANDIDATE_ONLY |
| GIGASTROKE East Asian AIS | `GCST90104545` | Stroke consortium AIS East Asian GWAS accession; precise source-specific Ncases/Ncontrols must be verified | CANDIDATE_ONLY |

Primary references:
- https://gwas.mrcieu.ac.uk/datasets/bbj-a-129/
- https://www.nature.com/articles/s41586-022-05165-3
- https://www.pgscatalog.org/score/PGS002725/ (GWAS Catalog accession references)
- https://www.ebi.ac.uk/gwas/studies/GCST90104545

**Important:** This list is not an adopted source manifest. GIGASTROKE combines phenotypes and ancestries; the EAS aggregate counts cannot be substituted without verifying the *exact AIS analysis*. Even a valid BBJ catalog number does not prove that the old in-house BBJ ALT effects were derived from that same release.

## Automatic provenance contract

The immutable research inputs stay in the existing server research directories. A separate optional review-only manifest `science_inputs/P3A_APPROVED_GWAS_SOURCES.tsv` is parsed by `scripts/is/harness_archived_evidence.py`.

Template: `docs/is/templates/P3A_APPROVED_GWAS_SOURCES.template.tsv`.

Each BBJ/EAS-AIS entry must document a source accession, approved absolute file path restricted to the server research tree, exact SHA-256 digest, GRCh37 build, case-control trait definition and counts, HTTPS reference, reviewer, date, and at least a claimed number of independent SNP matches. **A manifest with manually entered counts is never sufficient to mark the source as accepted:** the audit explicitly reports `MATCH_EVIDENCE_NOT_REPRODUCED_BY_THIS_AUDIT` until a separate independent SNP/allele/beta/p comparison is implemented and reproduced from the source files. The gate never publishes private source paths or individual records in GitHub Issues.

Missing registries yield `MISSING_APPROVED_REGISTRY`, `accepted_sources=[]`, `NO_GWAS_METADATA_ADOPTED`, and no change to coloc execution gates. Changed registry metadata triggers the next full audit.

## Next scientific engineering task

Implement the reproducible **genomic-variant matching** runner for separately downloaded, authoritative BBJ and GIGASTROKE AIS files: build/REF-ALT/strand check, at least an independent representative SNP panel, numeric beta/p matching with tolerances and orientation reversal checks. Require per-trait exact Ncases/Ncontrols, phenotype type and source-specific hashes. Then restrict P3A discovery to full-region QTL/GWAS overlap with a reviewed variant-exclusion report; conduct prior-sensitivity analysis only after all input provenance and clinical GWAS metadata gates pass.

No source genotype access, arbitrary remote command execution or live-analysis mutations are implied by this document.
