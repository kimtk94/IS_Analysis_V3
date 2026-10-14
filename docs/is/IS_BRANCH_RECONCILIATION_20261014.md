# IS branch reconciliation and next execution gates — 2026-10-14

## Verified GitHub source inspection
- Broad-discovery source branch: `research/is-broad-discovery-20261009`.
- Harness source branch: `feat/is-harness-dynamic-loop-20261013`.
- GitHub compare reports **diverged**, broad-discovery ahead 212 commits and behind 4 commits relative to harness (direction relative to comparison base). Do not assume fast-forward is possible.
- This is a GitHub/document audit, **not** a server execution or independent reanalysis.

## Existing documented universe (source: IS_2225_COMPLETE_UNIVERSE_WORKFLOW_READINESS_AUDIT_20261011.md)
- 2,225 distinct stable Ensembl IDs, 2,425 gene-region pairs.
- EAS/Japanese 869 pairs; EUR 1,556 pairs.
- Workflow readiness tiers: A 14; B 29; C 826; D 1,356. **These are data-readiness tiers, not causal ranking.**
- Historical 646 ABF assay records: 30 PASS and 616 MISSING_INPUT; **not** 30 unique genes.
- 40 candidate existing-data follow-ups outside G0022. Preserve the full universe.
- Historical stable anchors: FGF5, ALDH2, SH3PXD2A, COL4A2, COL4A1, C4orf22, CALHM2, NEURL1, INA; no automatic ranking bonus.

## Discovery gate
- Six EAS/Japanese GWAS files; 45 phenotype-specific proximity clusters, 30 coordinate-overlap components.
- Fifteen of 45 source regions have genome-wide-significant variants; seven of 30 merged components contain such a signal.
- Four of 30 lead coordinates overlap existing regional PVARs; allele audit 3 exact, 1 swapped, 26 unavailable.
- No component is yet asserted statistically independent; ±500kb gene-window annotation is not causal assignment.

## Immediate code review targets
1. `scripts/is/audit_existing_2225_evidence_priority_queue.py` and `tests/test_is_2225_workflow_priority_queue.py` — denominator and tier conservation.
2. `scripts/is/build_broad_discovery_v1.py`, `audit_broad_regions.py`, `merge_broad_intervals.py`, `audit_broad_ld_coverage.py`, `build_expanded_region_manifest.py`, `audit_lead_alleles.py`, `map_all_genes_gencode_v19.py`.
3. CI files `.github/workflows/is-broad-discovery-ci.yml` and `.github/workflows/is-evidence-guardrails-ci.yml`.
4. Harness experimental metric provenance, review gates and deterministic output manifests.

## Fail-closed acceptance tests
- Every row retains stable gene ID, region ID, ancestry, source, phenotype, evidence status, and explicit missingness.
- Assert 2,225 unique IDs and 2,425 pairs **only when using the frozen 2026-10-11 input snapshot**, not arbitrary future expansions.
- Assert tier totals = 2,225 and ancestry pair totals = 2,425 for that snapshot.
- Distinguish allele swapped from strand ambiguity; fail if 26 unavailable are silently counted as harmonized.
- Distinguish region proximity from LD independence and molecular colocalization.
- No promotion of NOT_ASSESSED/MISSING_INPUT to NEGATIVE or VERIFIED.
- No claims of replication without disjoint cohorts and documented sample overlap.

## SSH/Codex handoff
Codex SSH is reportedly configured, but this conversation's tools cannot execute commands in that Codex terminal. Run on the authenticated Codex environment:

```bash
git -C /srv/is-analysis/IS_Analysis_V3 status --short --branch
git -C /srv/is-analysis/IS_Analysis_V3 remote -v
ssh -T -o BatchMode=yes git@github.com
git -C /srv/is-analysis/IS_Analysis_V3 fetch origin
git -C /srv/is-analysis/IS_Analysis_V3 diff --stat origin/feat/is-harness-dynamic-loop-20261013..origin/research/is-broad-discovery-20261009
```

The path above is the research server path, not necessarily a Codex-local path. If Codex is on a different machine, locate its actual clone first. Do not reset, rebase, merge, or overwrite an uncommitted worktree until its status is recorded.

## Next commit boundary
Add source-backed regression tests and a versioned evidence matrix in the **broad-discovery code branch** after inspecting actual implementation. Cherry-pick only validated changes into the integration branch. Server analyses and thesis figures require separate execution evidence.
