# CKD TOP3 — Original Stanzick eGFRcrea four-SNP verification (2026-10-18)

**Primary source evidence:** [dedicated read-only GWAS workflow success 39124577872](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124577872). Original `metal_eGFR_meta_ea1.TBL.map.annot.gc.gz`; **8,844,847 rows** scanned. All four exact rsID and GRCh37 chr:pos hits found. This is independently read from original GWAS, not derived Stage2 files. Separate UKB-PPP original source verification [39122708873](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39122708873).

|Gene|rsID|UKB-PPP effect allele|pQTL A1FREQ|GWAS Allele1/Allele2|GWAS Freq1|GWAS Effect|GWAS StdErr|GWAS P.value|GWAS P.value.GC|1KG EUR minor freq|
|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|
|SDCCAG8|rs35646972|T|0.0549205|T/C|0.0984|0|0.0006|0.9951|0.9956|0.167992|
|SDCCAG8|rs6663222|C|0.0835498|T/C|0.8437|-0.001|0.0005|0.03995|0.0646|0.194831|
|INHBC|rs11834807|C|0.695821|T/C|0.3206|+0.0011|0.0003|0.0001957|0.0008064|0.409543|
|UMOD|rs11860566|T|0.00141033|A/T|0.8832|+0.0004|0.0009|0.6818|0.7122|0.135189|

**Important frequency convention:** UKB-PPP `A1FREQ` is frequency of pQTL `ALLELE1`, not necessarily minor allele. INHBC source A1FREQ=0.695821 corresponds to MAF=0.304179. For rs6663222, pQTL C frequency=0.0835498, GWAS C frequency=1-0.8437=0.1563. For rs11860566, pQTL T frequency=0.00141033, GWAS T frequency=1-0.8832=0.1168, reference LD T frequency=0.135189. These concordant allele identities rule out simple flip of strand/effect allele as sufficient explanation of the large frequency difference.

**GWAS effect direction note:** GWAS `Effect` corresponds to Allele1, whereas source pQTL effect corresponds to ALLELE1 (from original TAR). Normalize reference to pQTL effect allele before comparing signs of effects. For INHBC, the pQTL effect C is opposite GWAS Allele1 T, so aligned eGFR beta would be negative the original +0.0011. For UMOD pQTL T is opposite GWAS Allele1 A, so aligned beta would be -0.0004. For rs6663222 pQTL C is opposite GWAS Allele1 T, aligned beta +0.001. The original GWAS `Effect` values are rounded, so zero and sign interpretation should not exceed source precision.

**Source verification gate:** pQTL original frequency and allele ✅, GWAS original allele/effect and frequency ✅, individual pQTL vs LD reference frequency mismatch unresolved, LD–GWAS Z score consistency ⏳, UMOD SuSiE non-convergence source ⏳. None is grounds for claiming mediation or causal validation.

The old combined audit [39123850576](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39123850576) failed due to SSH Broken pipe in earlier repeated per-protein scans; dedicated GWAS audit avoids those redundant stages. No source files or LIVE settings modified.
