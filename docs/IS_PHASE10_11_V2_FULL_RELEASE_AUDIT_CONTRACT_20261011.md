# IS Phase10/11 — full 646 common-sdY V2 independent result release gate

**2026-10-11 | Independent validator prepared and tested, but V2 full-646 Google Colab has NOT yet run. Do not promote native 8,956/8,956 or quantitative sdY 646/646 before verified artifacts exist.**

## Why V2 needs an independent final audit

The original G0 numeric replay 646/646, G1 internal REF/ALT input consistency 646/646, and scalar MIN-N stress 646/646 are already complete. The following are different scientific questions:

1. **Native BBJ GWAS identity G1b:** Does the full union of **8,956 unique original BBJ GWAS SNPs** (including **647 outside the prior cache30**) match an independently generated effect-allele/BETA/AF/GRCh37-to-GRCh38 chain reference? The server's native GWAS source was streamed over 13,435,541 original associations to produce 56,189 source-verified SNP fingerprints. The old **8,309 cache30** subset passed, but the **647 new** SNPs remain UNVERIFIED pending Colab V2.
2. **Observed SNP-specific N–assisted common-sdY G2b:** Does full original 646-gene-tissue use of variant-specific observed eQTL N **only in a common quantitative-trait sdY estimate** materially change single-signal coloc posteriors? The previous 30-source direct pilot yielded maximum |ΔPP.H4|≈0.0000462 and zero threshold crossings; the 646 full-source result remains unexecuted.

## New direct one-cell V2

https://colab.research.google.com/github/kimtk94/IS_Analysis_V3/blob/feat/is-p0-evidence-20261010/notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_NATIVE_8956_REPAIRED_ONECELL.ipynb

- Source original index/master SHA256 frozen; 646 original per-SNP GTEx sources mounted from Google Drive.
- Native Sakaue/Kanai GWAS plus hg19ToHg38 chain reference fingerprint compressed file SHA256 pinned **acbec8c5b56177f4a51f0a5e6a21d46a190746abe8a1296a252b6bdfdb809a2a**. Verified GitHub raw link HTTP 200. No raw native GWAS BETA/EAF packaged for redistribution.
- Exact prior 30 source files have independent frozen SHA256s, and must reproduce native-SNP subset size 8,309.
- Full source set must contain exactly 8,956 unique GWAS SNPs; 647 must be novel relative to that 30-source subset. All full SNP effect/allele/liftover fingerprints must pass; otherwise the notebook writes an audit failure ledger to a separate folder, refuses scientific release and stops before R fitting.
- **After native 8,956 passes only**, runs all 646 original selected gene–tissue coloc cases with observed per-SNP QTL N entering the *single common sdY estimator*. Original per-SNP beta/SE, MAF and p1=p2=1e-4/p12=1e-5 unchanged. No actual per-SNP regression SE correction or matched QTL genotype LD.
- New output folder: MyDrive/MASTER_DEGREE/IS_COLAB/results/IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V2; no overwrite of earlier 646 original/structural/MIN-N analyses.
- Output manifest schema updated to **IS_PHASE10_11_646_VARIANT_N_SDY_COLAB_V2**, disambiguating the repaired V2 from the previously failed original V1.

## Independent post-Colab audit prepared now

- Code: scripts/is/audit_is_646_variant_n_sdy_v2_results.py
- The checker refuses missing/uncommitted outputs or existing audit output directory, and *independently* recomputes:
  - Manifest vs frozen historical 646-index/master/source R script and per-output SHA256.
  - Source SNP SHA256 equality against archived G0 for **all 646 original locus/tissue/gene keys** and same file names/SNP counts.
  - All 646 historical H3/H4 posterior vectors (absolute tolerance ≤1e-8) in V2 baseline and integrity of every variant-N sdY-adjusted H4 delta.
  - Full G1 646 source structural output status and no skipped/failed input.
  - The **8,956-row unique SNP native validation ledger** against every matching reference fingerprint; subset partition **8,309 cached + 647 newly verified**, no duplicates or unexpected statuses.
  - Per-threshold H4 crossing counts, maximum/median actual |ΔH4|, by-locus original SNP denominator, reported eqtl_an=2n checks and both code/number provenance boundaries.
- Standalone audit output: IS_646_VARIANT_N_SDY_V2_INDEPENDENT_AUDIT.json in new immutable local audit directory. Never to be generated unless complete source artifacts are present.
- Five tests in tests/test_is_646_variant_n_v2_independent_audit.py confirm corrected V2 schema and reference SHA256, no overwrite, missing-data fail-close, previous real 30-cached N-sdY numerical base result matches G0, and no raw BBJ summary-stat columns in the reference.
- **When V2 finishes**, pull its V2 results through connected Google Drive and run the auditor. Do not mark G1b 8,956 or G2b 646 as complete based solely on the Colab printout without saved-ledger readback.

## Interpretation rules

- Completion of G1b native allele + liftover means source association-direction/coordinate provenance, not proof of GRCh38 FASTA reference bases, GTEx donor LD, tissue-specific causality or gene mechanism.
- Completion of G2b common-sdY comparison is one model-parameter sensitivity test; NOT sample-level N correction, multiple independent causal variants, prior calibration, or tissue selection multiplicity.
- Entire analysis pertains to historical **4 BBJ ischemic stroke GWAS loci / 646 original gene–tissue tests**. The expanded **80 provisional GWAS windows / 2,225 positional genes** remain separate.
- Causal gene/protein mechanism **NOT_ESTABLISHED**; G2 matched molecular QTL LD and G3 stroke-specific orthogonal functional validation remain PENDING.
