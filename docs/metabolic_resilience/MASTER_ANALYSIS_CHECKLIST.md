# MASTER ANALYSIS CHECKLIST

## Core causal layer
- [ ] Stage3-C0R candidate-specific EUR LD reference PASS 8/8
- [ ] Stage3-C1 LD clumping PASS primary + stringent
- [ ] Stage3-C2 frozen independent IVs PASS
- [ ] Stage3-D harmonization and multi-SNP MR complete
- [ ] Multiple-testing families frozen and overlap metadata reviewed
- [ ] Continuous Steiger sensitivity reviewed; binary Steiger not overclaimed
- [ ] coloc ABF complete with p12 sensitivity
- [ ] Multi-signal coloc gate resolved
- [ ] Protein-altering/coding/splice artifact audit resolved

## Cross-ancestry layer
- [ ] CKB target/platform mapping reviewed
- [ ] KoGES/BBJ/TWB schemas reviewed
- [ ] Sample overlap explicit
- [ ] EAS-specific LD reference documented
- [ ] Concordance/non-significance/opposite/not-testable classes reported

## KoGES layer
- [ ] Feasibility GO/MODIFY decision
- [ ] GRS variant coverage complete or partial coverage explicitly reviewed
- [ ] Cox primary model
- [ ] Repeated MBI LMM
- [ ] PA interaction
- [ ] Missingness policy
- [ ] PH diagnostics and model diagnostics

## Interpretation / reproducibility
- [ ] Functional/tissue annotation complete
- [ ] No automatic candidate ranking from supportive annotation
- [ ] Environment snapshot
- [ ] Result SHA256
- [ ] Drive↔Git drift PASS
- [ ] MASTER_ANALYSIS_GATE PASS_CORE_READY
