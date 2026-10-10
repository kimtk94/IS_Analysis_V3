# IS Phase10/11 — native BBJ GWAS full-646 SNP universe repair

**2026-10-11. Previous V1 assertion fixed in NEW V2 notebook, NOT YET RUN. No gene causal claim.**

## Actual error

The first full original-646 variant-N-common-sdY Colab found INPUTS_INDEXED=646, AVAILABLE=646, MISSING=0; original index/master SHA256 unchanged; GWAS_UNIQUE_VARIANTS_ALL_646=8956 (new full SNP-union SHA256=4d287d3df4430fd49cc5046c3232b4de3b33970ac8107b1dc040deba05025431). It then incorrectly asserted that the full SNP universe must equal the **8,309** distinct SNPs derived from only **30 cached gene–tissue files**.

This is a **validation scope/denominator bug**, not a failed genetic association or coloc posterior check. Full original 646 sources have **8,956 unique SNPs**, which is a **net 647 SNPs beyond the previous cache30 subset**. Do NOT claim original GWAS/chain verification for those 647 before direct primary source audit. The first V1 Colab aborted before running any new R per-SNP-N common-sdY analysis. The historical G0 646 numeric replay, G1 646 source structural QC, and MIN-N 646 sensitivity remain unchanged.

## Recovery: original BBJ GWAS ZIP + independent chain-derived portable reference

The original Japanese BBJ ischemic stroke GWAS native ZIP (13,435,541 autosomal summary-stat rows) was re-streamed on the local server. Four original genomic windows plus 1 Mb flanks were independently mapped with the original hg19ToHg38 chain. **56,189 distinct native variants** passed primary source effect allele Allele2=GRCh37 ALT, unique chain, and GRCh38 allele-aware match-key derivation; invalid/noncanonical native variants were explicitly excluded by code and enumerated in the manifest.

Only a compressed **160-bit SHA256 fingerprint** derived from original native locus/ID/lifted match-key/BETA/AF_Allele2 is included in the Git branch, **NOT raw native GWAS beta/EAF values**. Canonicalization uses fixed 12 decimal places to preserve rigorous ~1e-12 source checks despite native-to-derived GWAS beta differences of order 1e-16.

- Native ZIP SHA256: a8032286cf2206a02c746b138a263930dc7d371c57b94352820005f8a32c8d2b
- Liftover chain SHA256: 5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d
- Compressed 56,189 reference fingerprint file SHA256: acbec8c5b56177f4a51f0a5e6a21d46a190746abe8a1296a252b6bdfdb809a2a
- Frozen 30 local source TSV SHA256 manifest: data/is_reference/native_verified_original_30_input_SHA256.json
- Direct local cross-check: **8,309/8,309 cached distinct GWAS SNP fingerprint matches, zero mismatches**, independently verified.
- Additional **647 distinct full-646 SNPs: still PENDING source fingerprint comparison**. Colab V2 checks these directly.

The portable reference does not independently check GRCh38 FASTA reference bases, GTEx donor genotype LD, or association causality.

## V2 one-cell Colab (READY, NOT EXECUTED)

https://colab.research.google.com/github/kimtk94/IS_Analysis_V3/blob/feat/is-p0-evidence-20261010/notebooks/is/IS_PHASE10_11_646_VARIANT_N_SDY_NATIVE_8956_REPAIRED_ONECELL.ipynb

This script: freezes all 646 original input and index/master hashes; downloads only the published 56,189 SHA256 native GWAS-chain fingerprint reference and verifies its digest; cross-checks exactly the original 30 cached source file hashes; reconstructs the full 646-input 8,956-GWAS-SNP union and 8,309-cache subset (647 new); checks all original 8,956 source GWAS native fingerprint matches with zero conflict; generates a separate failure ledger and stops if any SNP fails; **only then** runs the actual 646 original input per-SNP-observed-N-assisted common-phenotype-sdY recalibration under coloc 5.2.3. The latter does NOT recalculate the SNP SE or use matched GTEx-QTL LD.

New result directory: MyDrive/MASTER_DEGREE/IS_COLAB/results/IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V2, separate from old full-646 numerical replay, G1 structural input and scalar-MIN stress results. Existing outputs are not overwritten.

**The original first V1 notebook (IS_PHASE10_11_646_VARIANT_N_SDY_ONECELL.ipynb) contains the bad 8309==full union gate and must NOT be rerun. It remains only as reproducible failure evidence.**

## Scientific gates

- G0 original four-locus 646/646 numerical coloc H0–H4: COMPLETE, unchanged.
- G1a original 646 eQTL source structural consistency: COMPLETE.
- G1b original BBJ source effect allele + independent chain: prior 8309 verified, **full-646 8956 still requires V2**, 647 new pending.
- G2a original 646 scalar-N extreme MIN/FIRST/MEDIAN stress: COMPLETE.
- G2b original SNP-specific observed-N-assisted common sdY: local 30 source inputs PASS; full 646 V2 prepared, not run.
- G2c appropriate GTEx QTL genotype cohort LD, true SNP effective-N SE model, multi-causal signal fine mapping, tissue multiplicity and G3 stroke-specific cell mechanism remain pending. Expanded 80 GWAS windows / 2,225 gene universe is separate.
