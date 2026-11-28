# Skin regional molecular evidence and novelty decision — 2026-11-28

Status: HOLD_FULL_REGIONAL_AND_ANCESTRY_VALIDATION. No new MR/coloc, no new full data transfer. Existing real run readback and literature/source audit only.

## Verified source results

- Sun et al. Nature Communications 2026-08-15 DOI 10.1038/s41467-026-76575-4, https://pmc.ncbi.nlm.nih.gov/articles/PMC13578292/ . Data availability points to Zenodo 19638349 and author code. Published 36 significant pGenes; not a guarantee of complete all-tested assay universe.
- Run 45614961676 job 139207862687: official pQTL file 8,085,897 bytes MD5 5ed66ce34705c8aeff36b68162c73b64, 108872 valid rows, 36 author-matched protein keys; zero GCC2 or frozen chr2 window rows. No null inference. Run 45614087522 job 139204871712: author significant eQTL table 24788 rows, 282 traits, zero GCC2/window rows; reported paper row denominator 20999 differs. Retain discrepancy, do not substitute. Job 139204872255 lists second Zenodo eQTL file 70,906,052 bytes MD5 caca80b5a56f3fd9ceeb3b03b9b4cc07.
- Existing PR #362 evidence top3_fixed_marker_isolate_readback_20261128.json, run 46124756300 job 140916801245: full released selected-target eQTL file 870843 rows, 282 targets, one exact match of 16 frozen markers: rs4833267 G:A, ENSG00000138685, source slope .0527222, SE .0630581, nominal p .404752; effect allele orientation not newly certified. No GCC2 marker match. pQTL 0/16. Released file is not proven complete for all tested genes. Do not re-download merely to repeat point inventory.
- Author repository tree at 74b17724f23f6de983a0b76320f95f85439a289a lists significant tables and restricted overlap-named allpairs files; names alone do not establish full all-tested regional data. Request GCC2, EDAR, LIMS1 regional data and assay eligibility from authors if needed; no request sent.
- Run 45613505767 job 139206323250, GTEx v8 GCC2: 1072 exact common SNVs / 1492 GWAS / 2849 QTL, 71.85% and 37.63% coverage, all three contexts fail fixed 80% both-source gate. Spearman neg-log-p sun .167605, protected .237225, fibroblast .158168; top10 overlap zero. Descriptive, not Coloc, no independent SNP assumption.

## Novelty correction and additional sources

Endo et al. 2018, https://pmc.ncbi.nlm.nih.gov/articles/PMC5997657/ , 11311 Japanese women, previously associated rs6542772 with excessive hairiness and annotated GCC2/LIMS1 eQTL; rs3827760 with eyebrow thickness and rs260643 with hair morphology. Thus rs6542772–GCC2/LIMS1 skin-related phenotype context is prior work. Hairiness and wrinkle roughness are distinct; no same-site replication or shared causal variant established. Preserve LIMS1 as an explicitly labeled prior-literature candidate for future review, not a replacement for the existing frozen family. A new candidate panel would require a separate approved protocol before unused outcomes.

EBI eQTL Catalogue https://www.ebi.ac.uk/eqtl/Studies/ lists TwinsUK QTS000029 (433 donors across tissues, not skin N), distinct study candidate for contextual sensitivity, not East Asian matched source. Source access and exact tissue dataset not yet qualified. Michigan skin eQTL tables https://csg.sph.umich.edu/junding/eQTL/TableDownload/ are filtered at p<1e-5; not complete nominal data. GTEx v10 or skin sQTL/apaQTL may be contextual alternatives only after source/ancestry/assay eligibility; no indiscriminate expansion.

## Decision

Keep existing six-marker external context, 16-marker regional family, 1 Mb GCC2 window and >=80% dual coverage threshold. Priority: (1) original Korean study phenotype model/scale, sampling covariance, signed dosage LD and aggregate conditional estimates; (2) Chinese full nominal regional statistics, explicit tested-gene universe, sample size, tissue site, genome build and effect-allele convention; (3) genuinely independent same-site wrinkle phenotype; (4) only then Coloc and conditional sensitivity. If required sources unavailable, report descriptive findings or deprioritize manuscript. Do not claim GCC2 is novel or independent of EDAR. No new threshold, schedule, main merge, controlled data access or contact.
