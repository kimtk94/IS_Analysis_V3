# BBJ versus GIGASTROKE AIS: unrestricted lead reciprocal alleles and p values

**2026-10-11 KST.** All four chr4/chr10 study-specific regional GWAS lead positions were cross-looked up in the **original full canonical BBJ ischemic-stroke GWAS** and **GIGASTROKE EAS acute ischemic stroke GWAS**, without restricting to the EBI QTL-tested SNP set. Exact GRCh37 chromosome/position and unordered A/C/G/T allele pairs were compared. Effect allele sign was oriented to BBJ genomic ALT. Original datasets had opposite REF/ALT representation at several positions, so literal variant-ID equality is **not** an absence of replication.

| GRCh37 site | BBJ p | GIGA EAS AIS p | BBJ effect allele / beta | GIGA effect allele / beta | Allele pair matched | Normalized direction |
|---|---:|---:|---|---|---|---|
| 4:81181072 | **1.854491e-8** | **2.583e-8** | T / +0.0657711 | A / −0.0679 | A/T palindrome | **STRAND_UNRESOLVED** |
| 4:81182554 | **1.211904e-8** | **2.888e-8** | C / +0.0662372 | T / −0.0673 | C/T | Same direction |
| 10:105459834 | **8.762024e-10** | **7.792e-7** | C / −0.0677001 | T / +0.0568 | C/T | Same direction |
| 10:105599770 | **9.246082e-8** | **5.584e-11** | G / −0.0598125 | T / +0.0764 | G/T | Same direction |

**chr4:** Both study lead positions meet genome-wide marginal GWAS p<5e-8 in both studies. The observed different minima are therefore not evidence by themselves of discordant signals. GIGASTROKE lead chr4:81181072 is palindromic; the two reported opposite effect alleles with opposite beta signs are consistent with concordance if they refer to the same physical strand, but strand authentication is needed to make this a trustworthy oriented signal. 1000G EAS regional VCF reports 4:81181072 REF=A ALT=T, whereas GIGASTROKE canonical ID uses T:A. Earlier nonpalindromic PLINK extraction omitted this variant, so **signed chr4 lead-to-lead LD is not yet validated**. Ancestry-matched frequency plus reference sequence checks should be the next gate.

**chr10:** Study leads are at 10:105459834 and 10:105599770. Their variants are present in both original GWAS and reverse allele-pair notation is resolved. At both lead positions, effect-allele oriented GWAS signs agree, but each study's preferred lead falls just outside marginal genome-wide significance in the *other* GWAS. EAS504 reference LD r² between their sites is ~0.00737. This may reflect different regional signals, but formal in-study conditional GWAS and molecular QTL LD are needed to establish independence or colocalization.

**Interpretation:** The 4 SNP locations are not 4 independent tests. Sample overlap, differing outcome definitions, allele-strand provenance at palindromic loci, and the LD structure preclude independent replication claims. Original 646 historic coloc posteriors, 120 newer EBI–BBJ ABF sensitivity scenarios, and the 2,225-gene IS discovery universe remain unchanged.

**Reproducibility:** Source `scripts/is/audit_is_priority4_unrestricted_cross_alleles.py`; read-only source GWAS `/srv/is-analysis/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz` and `/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz`. Output `/srv/is-analysis/results/is/audits/is_priority4_unrestricted_cross_alleles_20261011_v1/LEAD_POS_ALLELE_CROSS.tsv`, JSON companion.

**Next:** Genotype-check the chr4 palindrome against an anchored reference, quantify lead-to-lead signed/r², and review marginal versus conditional signal interpretation. Formal COJO remains gated by absence of a sufficiently large cohort-appropriate LD reference.
