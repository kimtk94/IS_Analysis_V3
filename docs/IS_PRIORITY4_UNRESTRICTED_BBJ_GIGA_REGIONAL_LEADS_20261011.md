# IS BBJ vs GIGASTROKE EAS AIS — unrestricted regional GWAS lead audit

**2026-10-11 KST.** This deliberately removes previous QTD000161 eQTL-assay SNP-selection bias. All **13,435,531 BBJ IS GWAS rows** and **6,725,216 GIGASTROKE EAS AIS GWAS rows** were scanned. The positions are GRCh37. No molecular eQTL p-values, LD availability, or non-palindromic filter determined the GWAS lead.

| Region | BBJ interval variants | GIGASTROKE interval variants | BBJ lead and p | GIGASTROKE lead and p |
|---|---:|---:|---|---|
| chr4 80,087,872–82,356,697 | 8,514 | 4,371 | `4:81182554:T:C`, **1.2119e-8** | `4:81181072:T:A`, **2.583e-8** |
| chr10 104,113,502–106,353,104 | 11,979 | 4,266 | `10:105459834:T:C`, **8.762e-10** | `10:105599770:G:T`, **5.584e-11** |

The earlier shared-QTL SNP analysis found identical lead SNPs on chr4. This is **no longer true when all available GWAS SNPs are included**. The new chr4 lead is a **palindromic A/T variant** and is excluded from the earlier nonpalindromic reference export. This dataset's EAS504 full regional VCF contains the site `4:81181072:A:T` (`rs10857147`), in reversed REF/ALT notation compared with GIGASTROKE `4:81181072:T:A`. The unfiltered VCF genotype must be properly resolved before reporting signed r or between-study direction. Because it is A/T, strand can be ambiguous. **Signed chr4 lead–lead LD has NOT been authenticated and is withheld**.

The unrestricted chr10 GIGASTROKE SNP key `10:105599770:G:T` is the same genomic position and alleles as the earlier `10:105599770:T:G`, represented with **swapped allele order**. EAS504 REF/ALT-orientated genotypes for the latter and BBJ `10:105459834:T:C` show **r=+.085863, r²=.007372**, after ALT dosage harmonization, exactly as in the previous study. This is **reference-panel LD**, not a conditional GIGASTROKE or BBJ inference. With reverse allele order, the source's reported GIGASTROKE beta sign cannot be compared before reorientation, so we withhold a between-GWAS lead-direction claim.

**Important provenance restriction:** this identifies each study's lowest p-value *within the two fixed genomic intervals*, not the author's official GWAS-defined independent conditional association. The traits BBJ IS and GIGASTROKE EAS AIS have potentially different case definitions and unknown sample overlap. Multiple signals remain hypothetical until study-calibrated conditional tests or multi-signal colocalization are possible. Neither a locus-specific p-value threshold correction nor an independence claim is made.

## Reproducibility

Code `scripts/is/audit_is_priority4_unrestricted_gwas_leads.py`. Files at `/srv/is-analysis/results/is/audits/is_priority4_unrestricted_gwas_leads_20261011_v1/` include `IS_UNRESTRICTED_BBJ_GIGA_LEADS.tsv`, `IS_UNRESTRICTED_BBJ_GIGA_CROSS_LOOKUP.tsv`, and `IS_UNRESTRICTED_GWAS_LEADS_MANIFEST.json`.

The cross-GWAS lookup currently uses the **literal un-normalized variant ID**, so its missing cross-study p-value cells mean **ID representation not harmonized** and must not be interpreted as absence of the variant in that GWAS. Reciprocal allele orientation, palindromic status, reference genotype availability and GWAS sample overlap should be audited independently before next causal interpretation.

The 646 historical coloc posterior results, 120 new EBI–BBJ ABF sensitivity scenarios and original 2,225-gene IS discovery universe were all preserved.
