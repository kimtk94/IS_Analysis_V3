# IS chr4 palindromic lead: public 1000G reference LD validation

2026-10-11 KST. Variant site GRCh37 `4:81181072` (1000G Phase3 REF=A, ALT=T), rs10857147, is a palindromic A/T SNP. BBJ GWAS's lead is `4:81182554:T:C` (ALT=C). GIGASTROKE EAS AIS has a *regional* lead at `4:81181072:T:A`, written with reversed alleles. Both variants meet GWAS marginal p<5e-8 in both BBJ and GIGASTROKE.

We extracted exact per-sample genotype calls from previously acquired bounded 1000G Phase 3 VCF files at the palindromic site and merged by sample IID with existing validated GRCh37 ALT dosage of the BBJ lead. All genotypes were called and biallelic diploid (no imputation). This measures **genomic ALT_C versus ALT_T reference genotype LD**. It does not fully authenticate the GIGASTROKE-reported palindromic effect strand.

| Population | N | BBJ lead ALT=C freq | Palindrome ALT=T freq | r ALT_C vs ALT_T | r² |
|---|---:|---:|---:|---:|---:|
| EAS504 | 504 | .361111 | .353175 | **+.970575** | **.942015** |
| JPT104 (subset of EAS504) | 104 | .293269 | .298077 | **+.987348** | **.974855** |
| other EAS400 | 400 | .378750 | .367500 | +.967082 | .935249 |
| EUR503 | 503 | .267396 | .264414 | +.987046 | .974261 |

**Interpretation:** Very high pairwise LD (`r²>.94` in both entire EAS504 and JPT104) supports that the two different GWAS regional lead positions are likely tags of one shared high-LD haplotype **in these public reference samples**, not evidence by themselves of two independent conditional associations. The reported BBJ and GIGASTROKE effect allele pair/strand on the palindrome remains `UNRESOLVED` without independent verified strand/effect-allele-provenance. GWAS source reported BBJ risk effect T β=+.065771 and GIGASTROKE A β=−.0679. Naively opposite signs would suggest agreement when reoriented under a common forward strand, but that is explicitly **conditional on strand**, and not an attested result.

The per-sample public reference panel is **not** the BBJ or GIGASTROKE in-study LD and the small JPT N=104 is not a valid COJO reference. Strong reference-LD tagging supports shared regional association, **not** confirmation of causal genes FGF5 or C4orf22, nor SuSiE colocalization.

Code: `scripts/is/audit_is_chr4_palindromic_lead_ld.py`. Audit outputs: `/srv/is-analysis/results/is/audits/is_chr4_palindrome_lead_ld_20261011_v1/IS_CHR4_PALINDROME_REFERENCE_LD.tsv` and JSON, plus two transient public 1000G genotype extracts. Original 646 historical coloc outputs, 120 newer EBI–BBJ ABF scenarios and 2,225-gene broad IS universe were untouched.
