# IS: Ancestry-aware multi-signal readiness of 30 numerically replayed ABF assays

**Date:** 2026-10-11 KST. This is an existing-data-only audit, not an executed SuSiE-coloc or causal gene confirmation.

- Source: 30 numerically reproduced SNP-level ABF assays, representing 14 stable Ensembl genes.
- The 2,225-gene, 2,425 gene–region universe remains unchanged.
- Eight of the 30 original locus/gene/tissue records exactly join previously audited 1000 Genomes Phase 3 EAS504 *reference* LD alignment with all listed GWAS–QTL matched SNP IDs covered. A source-record count is not a sample count.
- **Zero** of 30 assays have GTEx donor-cohort signed QTL LD attested in the existing readiness files; **zero** are publication-ready for validated multi-signal SuSiE colocalization.
- BBJ AIS GWAS is Japanese/EAS, while GTEx molecular tissue cohorts do not supply equivalent disease- and QTL-cohort signed LD here. The EAS GWAS reference cannot safely be substituted as a GTEx QTL panel. Even when variant IDs agree, LD patterns can differ by ancestry, sample composition and quality filters.
- The 22 unlinked assays may have other LD assets elsewhere, but **this audit does not verify such assets**; they remain NO_ASSAY_SPECIFIC_SIGNED_LD_ATTESTATION rather than absence of all possible reference files.

Reproduce using `scripts/is/audit_is_30_joint_multisignal_readiness.py` with:
`IS_LEGACY_646_SNP_REPLAY_STATUS.tsv`,
`IS_2225_GENE_EVIDENCE_MATRIX.tsv`,
`IS_PRIORITY8_EAS_LD_READINESS.tsv`.

Outputs:
`/srv/is-analysis/results/is/audits/is_30_multisignal_ancestry_readiness_20261011_v1/IS_30_MULTISIGNAL_ANCESTRY_GATES.tsv`
and the JSON with original SHA256.

Remaining gates: study-relevant QTL signed LD, full cis SNP universe, independent in-cohort variant orientation, multiple-signal diagnostics, prior and tissue-selection sensitivity. No null biological finding or causal claim is inferred from unavailable data.
