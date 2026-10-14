# IS evidence audit and SSH continuation (2026-10-14)

## Scope and provenance
- Base branch: `feat/is-harness-dynamic-loop-20261013`
- Research root: `/srv/is-analysis/IS_Analysis_V3`
- This document is an execution plan, **not** a report of completed server validation.
- Preserve all local worktrees, raw data, and analysis outputs. Do not force-pull or run destructive cleanup.
- The older `docs/IS_MASTER.md` (2026-10-05) is a historical snapshot, not a complete 2026-10-14 state audit.
- Broadened discovery on `research/is-broad-discovery-20261009` must be reconciled before freezing a gene shortlist.

## SSH audit (read-only, run on the authorized server)
```bash
cd /srv/is-analysis/IS_Analysis_V3 || exit 1
printf '\n=== LOCAL ===\n'
git status --short --branch
git remote -v
git branch -vv
git log -8 --format='%h %ad %s' --date=short
printf '\n=== GITHUB SSH ===\n'
ssh -T -o BatchMode=yes -o ConnectTimeout=8 git@github.com 2>&1
printf '\n=== REMOTE REFS ===\n'
git ls-remote --heads origin 'main' 'feat/is-harness-dynamic-loop-20261013' 'research/is-broad-discovery-20261009' 'research/is-evidence-audit-20261014'
printf '\n=== KEY OUTPUT ROOTS ===\n'
for p in results/is results/IS data scripts/is tests; do
  if [ -e "$p" ]; then du -sh "$p"; else printf 'ABSENT %s\n' "$p"; fi
done
```

Do not interpret `ssh -T` exit code 1 with the GitHub greeting as an authentication failure; GitHub does not provide interactive shell access. Never print SSH private keys or environment tokens.

## Workstream A: candidate-universe reconciliation (P0)
- Compare historical four primary loci (FGF5, ALDH2, SH3PXD2A, COL4A2), the 2026-10-09 broadened discovery branch, and harness-generated candidates.
- Create a frozen, versioned `candidate_universe.tsv` containing locus ID, gene ID, genome build, GWAS source, ancestry, variant ID, evidence tier, inclusion reason, and input hashes.
- Record attrition by explicit step; do not equate non-detection in a single molecular QTL tissue with exclusion.
- Avoid cherry-picking by keeping the tested universe and all exploratory candidates.

## Workstream B: allele, LD, and fine-mapping audit (P0)
- Harmonize REF/ALT, effect allele, strand, rsID, genome build, positions, frequencies, and ambiguous palindromic alleles.
- Distinguish BBJ, EAS, JPT, EUR LD; label proxy-LD findings as sensitivity, not cohort-matched inference.
- Confirm conditional signals with locus-level diagnostics, regional plots, SuSiE convergence, credible set purity, and variant counts.
- Treat ALDH2 rs671 as coding/protein/metabolic candidate; avoid claiming bulk-eQTL mediation merely from expression localization.

## Workstream C: colocalization and MR (P1)
- Report locus/sample sizes, priors, PP.H0–H4, PP.H4/PP.H3, window selection and allele-matching counts.
- Track single-signal ABF independently from multi-signal SuSiE; absence of a credible set is unresolved, not a definitive null.
- For MR, record harmonization, instrument F statistic, IVW, weighted median, Egger where sufficiently powered, heterogeneity, pleiotropy and sensitivity analyses.
- A statistically significant MR association **alone** does not establish mediation or a shared causal variant.

## Workstream D: cell-type evidence (P1)
- Stratify brain, endothelial, pericyte, smooth muscle, macrophage, monocyte, and microglial datasets by sample source and disease state.
- For SH3PXD2A and COL4A2, prioritize disease-relevant vascular/immune regulatory annotation and ATAC peak-to-gene linkage.
- Separate expression detection, differential expression, cell-type localization, regulatory overlap, and genuine QTL colocalization.
- Record dataset versions, sample/cell counts, QC exclusions, donor counts, and replication.

## Workstream E: harness / reporting (P1)
- Each experiment must write config, input manifest and hashes, git SHA, seed, environment, test results, metrics, and explicit PASS/REVIEW/BLOCKED status.
- Keep REVIEW_REQUIRED when new metrics lack independent validation; success of a command is not scientific acceptance.
- Regenerate figures from frozen outputs, and publish only claims with evidence provenance.
- Cross-check manuscript figures against exact output paths and result checksums.

## Required final artifacts
1. `candidate_universe.tsv` and stepwise attrition audit.
2. `variant_harmonization_audit.tsv` with build/allele/LD consistency checks.
3. `locus_evidence_matrix.tsv` preserving negative and unresolved results.
4. Reproducible figure-generation manifest.
5. Thesis evidence summary tagged VERIFIED / EXPLORATORY / BLOCKED, with explicit limitations.

## Completion gates
- Local/remote commits compared without clobbering the worktree.
- All target inputs identified and checksum-verified.
- Candidate-universe scope reconciled with 2026-10-09 branch.
- Allele/LD diagnostics and coloc/MR scripts tested on representative fixtures.
- All automated findings distinguished from independently reviewed scientific claims.
- Server runs and publication only after these validations.

## Next operation
Read the server audit output and diff `feat/is-harness-dynamic-loop-20261013` against `research/is-broad-discovery-20261009`. Decide on merge/cherry-pick based on the file-level diff; never merge blindly.
