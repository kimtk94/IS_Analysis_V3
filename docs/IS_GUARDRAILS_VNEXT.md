# IS pipeline vNext: additive guardrails

This branch adds a **non-destructive, opt-in** audit layer. Historical BBJ IS GWAS, locus estimates, ALDH2 freeze and Phase 9D outputs remain unchanged.

## Use

```bash
python3 workflow/is_evidence_guardrails.py --mode variant --input VARIANT_CANONICAL.tsv --output /tmp/is_variant_qc.tsv
python3 workflow/is_evidence_guardrails.py --mode ld --input LD_MANIFEST.tsv --output /tmp/is_ld_qc.tsv
python3 workflow/is_evidence_guardrails.py --mode coloc --input COLOC_EVIDENCE.tsv --output /tmp/is_coloc_qc.tsv
python3 workflow/is_evidence_guardrails.py --mode functional --input FUNCTIONAL_EVIDENCE.tsv --output /tmp/is_functional_qc.tsv
python3 -m unittest tests.test_is_evidence_guardrails -v
```

## Source schemas and key gates

- **Variant** TSV must specify `build, chr, pos, ref, alt, effect_allele, other_allele, beta, se, p, eaf, ancestry`. `variant_key` includes **genome build**; REF/ALT cannot be inferred from effect/other alleles. Indels, palindromic SNVs, unresolved builds or reference allele mismatches are blocked for review, not silently flipped.
- **LD** requires `gwas_ancestry, ld_ancestry, ld_cohort, ld_file_verified`. A positive ancestry match is **necessary, not sufficient** for cohort match. Report JPT/1000G-EAS as reference proxies, never as BBJ cohort-matched LD.
- **Coloc** requires `method, ld_match, harmonization_qc, pp_h4, pp_h3, signal_pair_valid, finemap_status`. Strong support is possible only after a validated **SuSiE signal pair** with H4 >= 0.8 and H4 > H3. ABF-only results are always provisional. Missing LD or QC evidence prevents a strong classification.
- **Functional** requires `gene, mechanism_branch, gwas_qc, coloc_class, ld_match, functional_context, functional_support`. Tiers are **ordinal audit labels**, not calibrated probabilities. Mechanism branch distinguishes regulatory, coding/protein and vascular/immune hypotheses. Tissue expression alone does not establish causal mediation.

## Research-specific limitations

- FGF5: ABF/SuSiE convergence is promising but must not be called definitively causal without appropriate QTL LD and validated harmonization.
- ALDH2 rs671: coding/protein/metabolic path; no mandatory bulk eQTL mediation claim. Preserve its frozen allele, effect and coordinate.
- SH3PXD2A: assess macrophage/monocyte/microglia plus endothelial contexts.
- COL4A2: assess smooth muscle, pericyte, fibroblast and endothelial contexts; do not merge COL4A1 attribution.
- Phase 10 sc/snRNA and vascular ATAC integration still requires verified datasets, data contracts and full execution. New guardrails do **not** assert that this analysis has been completed.
- All four gates are opt-in until verified end-to-end on representative BBJ, eQTL and pQTL datasets; never silently overwrite canonical outputs.
