# IS BBJ and GIGASTROKE reciprocal lead p-values, GRCh37 allele orientation

**2026-10-11.** Scanned the full frozen BBJ ischemic stroke and GIGASTROKE EAS AIS GWAS summary tables for the four lead sites identified in the preceding unrestricted regional analysis. Resolved reversed variant ID allele order by **chromosome + base position + unordered allele pair** and checked GRCh37 `human_g1k_v37.fasta` reference base at each site. The GWAS reported effect allele was then oriented to the FASTA reference allele for both datasets.

| Position GRCh37 | BBJ p | GIGA AIS EAS p | BBJ beta (FASTA REF) | GIGA beta (FASTA REF) | REF-frequency difference | Important limitation |
|---|---:|---:|---:|---:|---:|---|
| 4:81182554 | 1.2119e-8 | 2.888e-8 | -.066237 | -.0673 | .000876 | nonpal, orientation verified |
| 4:81181072 | 1.85449e-8 | 2.583e-8 | -.065771 | -.0679 | .001019 | **A/T palindrome; strand origin unverified** |
| 10:105459834 | 8.7620e-10 | 7.792e-7 | +.067700 | +.0568 | .003486 | nonpal, orientation verified |
| 10:105599770 | 9.24608e-8 | 5.584e-11 | +.059812 | +.0764 | .001889 | nonpal, orientation verified |

The **FASTA genomic REF bases** at these coordinates are T, A, T, and T respectively. The raw GIGASTROKE variant IDs are `4:81182554:C:T`, `4:81181072:T:A`, `10:105459834:C:T`, `10:105599770:G:T`; BBJ genomic ordered pairs are the exact reverse at these four sites. Missing entries in the older literal-key cross-lookup table were due to ID ordering, **not absent observations**.

At the A/T palindrome, agreement of FASTA-REF-oriented effect estimates and allele frequencies is strong **supporting consistency but not independent proof of original DNA strand**. Mark this site **palindromic strand unresolved**. Do not aggregate or declare causal direction using this single variant without more source provenance or unambiguous haplotype anchors.

There are only **two loci**, not four distinct independent validations. Two lead positions at chr10 show weak EAS504 reference LD (r²≈.0074), but do not prove independence without study-calibrated conditioning. BBJ/GIGASTROKE phenotype differences, sample overlap and the EBI QTL donor-LD gaps remain. These reciprocal p-values are marginal, not conditionally independent replication or coloc.

Code `scripts/is/audit_is_priority4_unrestricted_lead_alleles.py`. Server results `/srv/is-analysis/results/is/audits/is_priority4_unrestricted_lead_alleles_20261011_v1/` (four-row TSV and JSON manifest). Preserve the original 2,225 IS genes, historic 646 coloc fits and EBI-BBJ 120 ABF runs.

**Next gate:** resolve chr4 palindromic strand using source study allele/FASTA alignment metadata and extract the chr4 pair from raw EAS504 genotype VCF without nonpal filter to measure lead–lead r². This is a technical LD audit, not a substitute for conditional GWAS.
