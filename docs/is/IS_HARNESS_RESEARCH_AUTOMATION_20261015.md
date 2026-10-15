# IS research Harness operations — 2026-10-15

## Purpose
Serve reproducible evidence-readiness checks for the ischemic-stroke research program without overwriting source GWAS/QTL files. This is **not** a fresh MR, colocalization, or causal inference analysis.

## Allowed GitHub Issue requests
Create an Issue authored by \`kimtk94\` with exact label \`is-harness-queued\` and exact JSON body:

| Body | Meaning | Output |
|---|---|---|
| \`{"operation":"eas5-audit"}\` | Stream and validate EAS5 harmonized GWAS/QTL records | \`eas5_audit.json\` |
| \`{"operation":"candidate-2225-audit"}\` | Validate 2,225 gene frozen QTL readiness TSV | \`candidate_2225_audit.json\` |
| \`{"operation":"figure-manifest"}\` | Render two numerical QC figure specifications from previously validated audits | \`research_manifest.json\` |
| \`{"operation":"research-refresh"}\` | Run both audits, then produce figure specifications | All three files |

No arbitrary shell commands, paths, custom scripts, or user-supplied parameters are accepted. Additional preexisting operations (\`status\`, \`inventory\`, \`headers\`, \`tests\`, \`mkdir-output\`) are preserved.

## Fixed scientific inputs

- \`/srv/is-analysis/results/is/audits/is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1/IS_EAS5_OFFICIAL_QTL_ALT_GWAS_HARMONIZED_VARIANTS.tsv\`
- \`/srv/is-analysis/results/is/audits/is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1/IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv\`
- \`/srv/is-analysis/results/is/audits/is_eas5_qtl_gwas_snp_ref37_20261011_v2/IS_EAS5_GENE_TISSUE_GWAS_SNP_OVERLAP_QC.tsv\`
- The frozen 2,225-gene TSV in commit \`fdb1983c436d8f438a606c14c3322bb157349d77\`, path \`docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv\`

The script enforces **five** EAS regions, 2,225 unique Ensembl gene IDs, and the four expected pilot readiness category sizes (1356, 629, 149, 91). Any unexpected change is a validation error, not a reason to silently change a denominator.

## Generated artifacts

Stored atomically under \`/srv/is-analysis/results/is/github_actions_runs/evidence/\`:

- \`eas5_audit.json\`
- \`candidate_2225_audit.json\`
- \`research_manifest.json\`

Source files are never edited. Each audit includes source hash(es), provenance and a clearly noncausal interpretation. QC counts are **variant-gene-tissue record counts**, not automatically unique genomic loci.

On \`research-refresh\` completion, a sanitized summary is added to its GitHub Issue comment as \`IS_EVIDENCE_JSON_V1={...}\`; the Vercel IS page fetches it through \`/api/is-evidence\`. The browser receives only whitelisted aggregated counts and text; it does not receive server raw rows, API credentials, or paths.

## Tests and controlled deployment

- Offline tests: \`python3 -m unittest discover -s tests -p 'test_is_harness_research_jobs.py' -v\`
- Workflow: \`.github/workflows/is-harness-research-jobs-ci.yml\`
- Remote installation workflow: \`.github/workflows/is-harness-safe-worker-update.yml\` on approved \`main\` changes. It validates and transfers only two allowlisted Python files under \`~/is-harness/\` using protected SSH environment secrets. It does **not** change the checked-out research branch.
- Existing \`systemd --user\` timer invokes the worker every ~60s. The IS web polls new evidence/queue state every ~30s.

## Known limitations
- No real-time raw journal streaming or background science re-analysis.
- A single worker is required; label-based claim is not transactional across hosts.
- A failed scientific audit will create a failed Issue record, not fabricated findings.
- Issue START/END and duration are execution metadata only, not evidence of scientific validity.
- Persistent artifact JSON on the server is not a centralized research DB. The web source is the most recent successful sanitized GitHub Issue result.
- GitHub Issues and Vercel access must remain appropriately protected.
