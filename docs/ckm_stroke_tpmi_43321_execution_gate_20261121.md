# TPMI 433.21 targeted SNP and GIGASTROKE gate — 2026-11-21

## Verified evidence (do not promote to scientific GO)

- Original 58,302,393-byte PLINK2 file, SHA256 `ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62`, 1,739,048 GWAS records.
- Previous original-source run [#44872892525](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44872892525): offline-negative-tests and full-original-genomewide-source jobs completed successfully; original aggregate artifact ID 14694881855 still unexpired when inspected.
- Independent OBS_CT audit run [#44873247155](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44873247155): artifact ID 14695131212 still unexpired when inspected.
- One source variant meets P <= 5e-8, but its identity and GRCh37 mapping are **not yet verified in an observed run**.
- PLINK2 OBS_CT = 204099 at every recorded variant, while phenotype page lists 9249 cases + 304660 controls = 313909. The case/control breakdown among the 204099 regression observations is UNKNOWN; no binary N_eff is approved.

## Current executable pipeline

Workflow: `.github/workflows/ckm-tpmi-43321-lead-grch37.yml` in branch `research/ckm-stroke-staged-actions-20261111` (PR #284, not merged).

1. PR and push: test source one-hit parsing, bidirectional assembly mapping, and GIGASTROKE lookup request offline. This is **code contract validation**, not primary-source result validation.
2. Explicit original-source execution (workflow_dispatch `run_original_source=true`, or a push whose commit message includes `[tpmi-lead-original]`): verify the two previous source-run pins, download original TPMI 58MB **once**, SHA-check, scan all 1,739,048 rows, save ONLY the unique genomewide-significant variant.
3. Using validated single SNP: query official Ensembl forward and reverse assembly maps and original/mapped reference bases; reject ambiguous mapping, invalid position, unequal REF, and unsupported indels.
4. Publish GRCh37 lookup **target only**; GIGASTROKE effect is not downloaded or inferred and scientific result remains HOLD.

## Proof needed to move to gate 3 PASS

- A new GitHub Actions run ID with SUCCESS for `extract-pinned-original-lead`, `bidirectional-build-and-allele-map` and `prepare-gigastroke-grch37-target`.
- Archived outputs `tpmi_43321_GRCh38_single_lead_public_GWAS_QC.json`, `tpmi_43321_GRCh37_Ensembl_bidirectional_map_HOLD.json` and `tpmi_43321_GRCh37_GIGASTROKE_target_only_HOLD.json`.
- Record exact chr/pos/REF/ALT/A1, source P and A1 logOR, Ensembl roundtrip, mapped REF and target position; do not infer replication from a coordinate alone.

## Proof needed for gate 4 scientific consideration

- Independently sourced GIGASTROKE *original* target-row chr/pos/alleles, effect allele, beta or logOR, SE and P; handle flips, complements, and strand-palindromic SNPs with frequency evidence.
- Document each cohort's ascertainment and participant overlap. EUR and EAS effects must not be called independent simply because ancestry labels differ.
- Obtain actual PLINK2 433.21 per-analysis cases and controls; do not equate `OBS_CT` to effective sample size or use 2025 PheWeb/2026 study sample counts interchangeably.
- No LDSC, MR, coloc or independent replication claim until each corresponding analysis-specific gate is met.

## Current disposition

- Technical full-file QC: **PASS** from prior pinned run.
- New one-SNP extraction, Ensembl liftover, GIGASTROKE matched effect: **UNVERIFIED, HOLD**.
- Scientific independent replication and effective-N: **HOLD**.
- New workflow run: **not confirmed** from available GitHub commit workflow-run/status lookups (which do not enumerate all events). Do not interpret an empty lookup as proof of zero runs.

## Additional targeted AIS source matching implementation

Code: `scripts/ckm_stroke/match_tpmi_43321_gigastroke_ais.py`, with offline cases `tests/test_ckm_tpmi_43321_gigastroke_ais.py`.

The main TPMI workflow now includes an **optional** manual input `run_gigastroke_ais=true`, guarded by `run_original_source=true`. It downloads the official 283,482,604-byte GCST90104540 GRCh37 gzip **once**, verifies publisher MD5 `1bfe4ae8a5042fb24cb0a562b05b0d2f` after scanning all rows, enforces exact mapped coordinates and uniqueness, and extracts only one target association and raw effect orientation.

- On no target locus, ambiguous allele pair, duplicates, wrong published hash, or invalid numeric input, the job returns HOLD or fails closed.
- On a valid source coordinate plus non-palindromic allele match, the single-locus **technical** gate can PASS; sample independence, valid binary N_eff, independent replication and LDSC remain HOLD.
- No successful run, target rsID, numeric original AIS beta or allele direction has been observed from this new code at the time this section was committed.
- A draft/feature-branch workflow may not be available through GitHub's manual Run workflow button until registered on the default branch. This is a deployment constraint, not proof the code passed its tests.

## Verified executed TPMI source and GRCh37 mapping — 2026-11-21

**Source-backed GitHub Actions run:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44878942971

Individual GitHub job API conclusions:
- `offline-allele-assembly-tests`: **SUCCESS** (job 136769756519).
- `extract-pinned-original-lead`: **SUCCESS** (job 136769807292), artifact 14697201895, full pinned original SHA verified, 1,739,048 records, one source P <= 5e-8.
- `bidirectional-build-and-allele-map`: **SUCCESS** (job 136770356497), artifact 14697022641, Ensembl forward/reverse and both assembly REF bases checked.
- `prepare-gigastroke-grch37-target`: **SUCCESS** (job 136770422991), artifact 14696224277.
- `source-verified-gigastroke-ais-one-locus`: **SKIPPED** (not yet harmonized with GIGASTROKE original).
- `report-ci-status`: **FAILED** (job 136770469771) with `gh: Resource not accessible by integration (HTTP 403)` when attempting a PR #284 comment, unrelated to original GWAS parsing or Ensembl mapping. The overall run is therefore failure even though the four scientific-data preprocessing jobs passed. The main workflow was corrected in follow-up commit `288c109a4cde7ec9772e658d8a72b75c68fd5dd4` to write only to `GITHUB_STEP_SUMMARY`.

Original pinned lead extracted from real official TPMI PLINK2 GWAS:
- `rs187759247`, source `chr19_15462568_T_C`; GRCh38 `19:15462568:T:C`; A1 effect allele C.
- Source A1 frequency `0.0212101`, OR `1.68467`, log(OR) `0.5215696989267228`, log(OR) SE `0.0662893`, Z `7.86808`, P `3.60115e-15`, OBS_CT `204099` (not effective N).
- Successfully mapped to **GRCh37 `19:15573379:T:C`**, effect A1 C, reverse liftover and source/mapped reference alleles verified.
- No SNP effect from GIGASTROKE AIS GCST90104540 has yet been checked at the locus, and independence/actual case/control N remain unverified. **Scientific HOLD** for replication, LDSC, MR and coloc.
- Updated narrow milestone: original source integrity and coordinate mapping TECHNICAL PASS; GIGASTROKE comparison/outcome causal interpretation HOLD.

## Executed GIGASTROKE AIS exact-locus check — 2026-11-21

- Executed run: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44880059442
- Offline tests, source extraction, Ensembl liftover and target generation: **SUCCESS**.
- The `source-verified-gigastroke-ais-one-locus` job (`136773836212`) completed with **HOLD_NO_EXACT_GIGASTROKE_COORDINATE** for `19:15573379:T:C`; the generated HOLD artifact was `14698280048`.
- This HOLD state is emitted by the source parser only after completing the compressed input scan and matching the pinned original file size and published MD5. No beta, standard error or risk-allele alignment can be claimed from an exact target SNP in this GWAS file.
- **The exact locus is absent in this specific GIGASTROKE AIS EUR release**. This is not evidence of no biological effect, and it is not a failure of coordinate mapping.
- Do not run exact-coordinate retries indefinitely. Before any proxy analysis, inspect available nearby variants in the source, define a defensible East Asian LD panel with declared build/population, validate proxy r² and allele phase, confirm allele frequency and independence, then harmonize the proxy effect. A nearby physical locus alone is not a proxy.
- AIS source-backed direct SNP replication: **HOLD_NO_EXACT_COORDINATE**; CKM/IS LDSC, MR and coloc: **HOLD**, as previously.

## East Asian LD proxy source evidence — 2026-11-21

- [GitHub Actions #44881325356](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44881325356) **SUCCESS**: Ensembl GRCh37 1000 Genomes phase 3 LD query on `rs187759247` returned zero `r²>=0.80` rsID proxies in each of CHB, CHS, JPT, CDX, KHV (all 5 API queries succeeded).
- An additional check was required because a 200/empty LD response alone does not demonstrate that the lead rsID is in the underlying LD panel.
- [GitHub Actions #44881838774](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44881838774) **SUCCESS**: independently queried Ensembl `variation/human/rs187759247?pops=1` for GRCh37 mappings and 1000 Genomes named EAS population allele frequencies, and repeated all 5 population LD queries. The script returned `HOLD_NO_HIGH_LD_PROXY_IN_QUERIED_REFERENCE`, `variant_presence: PASS_VARIATION_RECORD_READ`, each population `PASS_SOURCE_QUERIED`, and `candidates:0`. Its narrow classification requires a GRCh37 `chr19:15573379` mapping and at least one named EAS population AF entry. Artifact `14698122275`.
- Conclusion is limited to **no `r²>=0.80` proxy found by this Ensembl GRCh37 1000G phase-3 EAS reference query**; it does not rule out all proxy choices using alternative reference panels or EAS-specific outcome GWAS and cannot be called successful replication.
- GIGASTROKE GCST90104540 EUR AIS complete official source scan: no exact `19:15573379` row in two verified executions `44880059442` and `44880521973`. Direct AIS matched effect size absent; proxy evidence absent at cutoff; **scientific replication HOLD**.
- Next research design decision: seek alternative AIS outcome dataset with lead SNP coverage and ancestry/sample-independence support, or explicitly declare this locus not replicable against current GIGASTROKE AIS release. Do not reduce r² cutoff without pre-specified sensitivity analysis.

## Alternative stroke GWAS candidate access — 2026-11-21

**Scientific motivation:** TPMI `rs187759247` (`GRCh37:19:15573379:T:C`) is not present at the exact coordinate in source-MD5-verified GIGASTROKE EUR AIS `GCST90104540` (two separate Action runs) and no `r²>=0.80` proxy was returned by GRCh37 1000G phase-3 EAS query for CHB/CHS/JPT/CDX/KHV. Direct GIGASTROKE replication is HOLD.

**Candidate BBJ East Asian:** IEU OpenGWAS dataset `bbj-a-129`, ischemic stroke, HG19/GRCh37, 17,671 cases, 192,383 controls, 8,885,705 reported variants; these counts come from the external dataset catalog, not a per-SNP source scan. Public catalog URL: https://gwas.mrcieu.ac.uk/datasets/bbj-a-129/

**Candidate FinnGen:** Release R12 `I9_STR_EXH` (ischaemic stroke excluding haemorrhages); official public download manifest https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv ; R12 sumstats are GRCh38, necessitating a separately checked GRCh37->GRCh38 mapping and effect-allele harmonization.

**Executed GitHub Actions:**
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44882159836 — SUCCESS; BBJ catalog and FinnGen R12 manifest HTTP 200; BBJ brittle source-text numeric metadata matcher failed, FinnGen first 500k scan didn't find target name; no per-variant outcome effect evaluated.
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44882209478 — SUCCESS; BBJ and entire FinnGen 773,612-byte manifest accessible, but direct page-text identifier matcher was false and the exact string `I9_STR_EXH` was absent from this manifest. BBJ actual result file and exact SNP coverage **NOT_QUERIED**. These are honest metadata/access diagnostics; not a scientific negative replication.

**Next evidence gate:** Verify current BBJ full summary/VCF endpoint is publicly retrievable and has actual GRCh37 rs187759247/19:15573379 effect allele, beta, SE, p-value and source provenance; independently evaluate overlap with TPMI before calling replication. FinnGen release/endpoint path should be reconciled with the public manifest, not constructed from assumptions.

## BBJ original access and locus-specific outcome evidence — 2026-11-21

- [Action #44882542342](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44882542342) **SUCCESS**. Public PheWeb.jp downloads catalog disclosed official 2021 BBJ IS Sakaue/Kanai source `https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.IS.v1.zip` (22,664 cases, 152,022 controls). This is distinct from 2019 OpenGWAS `bbj-a-129` (17,671 cases, 192,383 controls, 8,885,705 variants); never conflate study samples or file provenance.
- [Action #44882606387](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44882606387) **SUCCESS**: public original BBJ ZIP returns HTTP 200 and size `1,611,973,884 bytes` and supports HTTP byte ranges. Official GRCh37 exact locus PheWeb API at `https://pheweb.jp/api/variant/19:15573379:T:C` returns HTTP 200. Variant page confirms `rs187759247`, genome locus `chr19:15573379:T:C`.
- [Action #44882662184](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44882662184) **SUCCESS**: authentic PheWeb SNP API returned `rs187759247` and three PheWAS phenotype rows, but **zero `phenocode == IS` (ischemic stroke) rows**. Original raw PheWeb API result therefore contains no stroke outcome beta/SE/p for this SNP. Artifact `14698846110`; status `HOLD_NO_BBJ_PHEWEB_IS_ASSOCIATION_FOR_TARGET`.
- A variant's presence in PheWeb for other traits must **not** be used as an ischemic stroke replication association. No 1.61 GB transfer is justified without an explicit need to cross-check PheWeb completeness; the API evidence is clearly recorded.
- OpenGWAS `bbj-a-129` remains **NOT_QUERIED** at the per-variant level because official automated OpenGWAS API/VCF access requires a JWT; the 2021 PheWeb result **does not imply** that the distinct 2019 BBJ GWAS lacks the rsID.
- Scientific replication remains **HOLD**: TPMI `rs187759247` reliably extracted; GIGASTROKE EUR AIS exact coordinate missing (2 full-file scans); Ensembl GRCh37 1000G EAS r²>=0.8 proxy missing across CHB/CHS/JPT/CDX/KHV; 2021 BBJ PheWeb AIS association missing. Any future GO requires actual authentic effect summary in a declared non-overlapping outcome GWAS with valid allele orientation.

## 2019 BBJ authentic cerebral-infarction summary source — 2026-11-21

- [Action #44883013507](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44883013507) **SUCCESS** in technical runner: published JENGER source/home endpoints timed out under both HEAD and GET; NBDC `hum0014/v32` website accessible (HTTP 200), but that version's generic page is not itself proof of 2019 AIS source coverage.
- Official NBDC `hum0014-v17` release (2019-10-08) contains individual unrestricted GWAS disease datasets. [Action #44884266698](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44884266698) queried **all 40/40** published dataset pages and uniquely identified `NHA000096`, `Cerebral infarction I639`, linked to original ZIP `hum0014.v17.CI.v1.zip`. Source dataset page: https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000096 .
- Published `hum0014_v17_README.txt`: effect allele is **alternative** `Allele2`; coordinates **hg19/GRCh37**; SAIGE v0.29.4.2; `*.auto.rsq07.mac10.txt.gz` for the both-sex autosomal analysis, `Rsq>=0.7`, `MAC>=10`; original data include `CHR POS SNPID Allele1 Allele2 BETA SE p.value AF_Allele2 Rsq MAC`.
- [Action #44884351810](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44884351810) **SUCCESS**: original NBDC ZIP size **2,252,105,310 bytes**, remote byte-range and ZIP central-directory integrity PASS. Five verified members: a directory; X chromosome main `IS.chrx.rsq07.mac10.txt.gz`; X chromosome sex-stratified; autosomal sex-stratified; and main autosomal `IS.auto.rsq07.mac10.txt.gz`. Main auto member compressed within outer ZIP **906,186,722 bytes**, original inner gzip member size **906,046,078 bytes**. No allele effect inferred from ZIP metadata.
- [Action #44884446188](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44884446188) **started** source original `IS.auto.rsq07.mac10.txt.gz` bounded HTTP range extraction and verifies original ZIP member CRC32, full gzip stream, GRCh37 `19:15573379` REF T ALT C, and source raw beta/SE/p. Scientific HOLD pending job log and cohort-independence confirmation.
- These facts are independent of older 2021 PheWeb 22,664-case AIS release and **must not be combined as though sample sources or summary datasets were interchangeable**. 

## BBJ 2019 primary cerebral-infarction source scan completed — 2026-11-21

- **Real source-backed Actions evidence:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44884446188 (attempt 2), job `136788491844`, **SUCCESS**. Corrected a first-attempt parser bug (mistaken literal escaped whitespace regex, not a biological no-match); research code patch `59176384e717c4f7a3d94c50b7409a6beadda469`, rerun using GitHub connector's direct `rerun_workflow_job` operation.
- **Verbatim original header:** `CHR POS SNPID Allele1 Allele2 AC_Allele2 AF_Allele2 N BETA SE Tstat p.value p.value.NA Is.SPA.converge varT varTstar AF.Cases AF.Controls Rsq MAC`. Source is official NBDC `NHA000096` `hum0014.v17.CI.v1.zip` -> `IS.auto.rsq07.mac10.txt.gz`, GRCh37/hg19.
- **Scanned precisely 8,678,632 source association rows** after full compressed member/ZIP validation. Target mapped GRCh37 `19:15573379:T:C` matched **zero coordinate rows**. Result `HOLD_NO_EXACT_BBJ_2019_CI_COORDINATE`, explicit `matches:0`, `association:null`, `error:null`, `scientific_status:HOLD`. Attempt-2 artifact ID `14699471967`. The successful job is a *technical completion*; scientific independent replication has not been demonstrated.
- Source filtering: `IS.auto.rsq07.mac10.txt.gz` is original published autosomal summary association after `Rsq>=0.7` and `MAC>=10` filtering, not an unfiltered complete catalog. A missing row can reflect QC/MAF/imputation/coverage criteria, not a null genetic association.
- Independently, MD5-verified GIGASTROKE EUR AIS `GCST90104540` exact GRCh37 coordinate not present in two full source scans; PheWeb.jp 2021 BBJ exact PheWeb variant returned three other phenotypes but `IS` 0 rows; GRCh37 1000G phase-3 five EAS populations have no validated high-LD proxy with `r²>=0.80`.
- **Decision:** the single lead TPMI `rs187759247` has no confirmed independent AIS outcome effect in the three evaluated public sources. Keep AIS replication, MR, coloc and LDSC **HOLD**; do not mix unrelated phenotypes, switch to physical-nearest-neighbor proxy, or claim null effect. Choose a new well-powered outcome GWAS with demonstrated variant coverage and no cohort overlap, or move to a genome-wide/multi-locus scientific design that is not contingent on this one rare SNP. If OpenGWAS `bbj-a-129` is reattempted, prove its exact source identity and SNP inclusion separately, not by assuming NBDC equivalence.

## Genome-wide outcome pivot: GIGASTROKE 30 accessions and five EAS stroke outcomes (2026-11-21)

The single TPMI lead rs187759247 is not a sufficient standalone outcome-replication design. The next workstream prioritises **ancestry- and subtype-specific original GWAS coverage**, matching alleles, SNP intersection and valid original metadata, with LDSC/MR/coloc scientific HOLD pending those gates.

**Source-backed GitHub Actions:**
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44888185252 attempt 1 **SUCCESS**, 30/30 original publisher GRCh37 sumstats returned HTTP HEAD PASS, but the old GWAS Catalog REST endpoint gave HTTP 410 Gone for all 30. The original REST API v1 was retired in August 2026, not proof that metadata was missing.
- Same Action **attempt 2 SUCCESS** using official REST API v2; 25/30 study metadata HTTP PASS, 5 transient 429 Too Many Requests, while 30/30 original files remained available. Technical registry artifact ID 14701061482. This registry does not itself verify any SNP effects or sample independence.
- https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44888444640 **SUCCESS**, direct v2 phenotype EFO ontology and exactly described East Asian case/control N were validated for **5/5** studies. Artifact ID 14701515076.

| EAS discovery outcome | GWAS Catalog accession | Official EFO trait label | Cases | Controls |
|---|---|---|---:|---:|
| AS (any stroke) | GCST90104544 | stroke disorder | 27,413 | 237,242 |
| AIS (any ischemic stroke) | GCST90104545 | Ischemic stroke | 19,032 | 237,242 |
| CES (cardioembolic) | GCST90104546 | cardioembolic stroke | 926 | 237,242 |
| LAS (large artery) | GCST90104547 | large artery stroke | 1,735 | 237,242 |
| SVS (small vessel) | GCST90104548 | small vessel stroke | 5,532 | 237,242 |

The study official AIS metadata separately lists replication: 8,304 EAS cases and 174,855 EAS controls; do not pool discovery and replication uncritically.

**Priority outcome: GCST90104545 EAS AIS.** Smaller subtypes are secondary (especially CES n=926), and all have shared EAS controls. Cross-ancestry GIGASTROKE includes EAS samples, so its results are not independent of EAS-specific GWAS. An EAS discovery outcome is not automatically independent of the original TPMI cohort or suitable for MR.

- Actual full original EAS AIS source QC has been started: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44888595169 . Until its result is inspected, original genomewide variant coverage, gzip integrity, chromosome distribution and rs187759247 exact presence remain NOT VERIFIED.
- No approved per-SNP binary sample effective N, official cohort overlap clearance, TPMI genomewide GRCh38-to-GRCh37 harmonization or signed harmonized IV has yet been verified. **Scientific status: HOLD.**

## Official EAS AIS whole-file QC and exact TPMI phenotype clarification (2026-11-21)

### Whole-file EAS AIS original verified
- GitHub Actions https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44888595169 **SUCCESS**; job 136799230332; artifact 14701455627.
- Public original: EBI GWAS Catalog \`GCST90104545_buildGRCh37.tsv.gz\`. Full gzip stream was consumed and gzip CRC/integrity checked, the official URL and Content-Length matched the actual downloaded bytes, and a reproducible **computed** full-source SHA-256 was recorded: \`c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3\`. A separately published checksum was **not** independently matched; do not claim it was.
- Exactly **6,725,216 variants** scanned in the source; **143,669 chromosome 19 records**; **0** rows at the verified TPMI GRCh37 coordinate \`19:15573379\`. Log technical status \`PASS_SOURCE_VERIFIED_NO_EXACT_LOCUS\`. This does not imply absence of biological effect.
- Distinguish full genomewide GWAS usability from *single rare TPMI SNP* replication failure. EAS AIS outcome access and full source schema are verified; original effect of TPMI rare lead absent; genomewide GRCh38↔GRCh37 TPMI overlap and complete allele harmonization have **not** been performed.

### Phenotype identity is not ICD-9-CM identity
- GitHub Actions https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44888783440 **SUCCESS**; job 136799805352; artifact 14700697258.
- Actual TPMI PheWeb source page \`https://pheweb.ibms.sinica.edu.tw/pheno/433.21\` and corresponding \`pheno-filter\` were separately checked via live HTTP 200.
- **TPMI Phecode 433.21: "Cerebral artery occlusion, with cerebral infarction"**, **9,249 cases**, **304,660 controls**, circulatory-system category. These figures are direct phenotype-page metadata, not \`OBS_CT\` effective N.
- **ICD-9-CM code 433.21 is different**: "Occlusion and stenosis of vertebral artery with cerebral infarction." The TPMI identifier here is Phecode, not ICD-9-CM; do **not** use the ICD-9 definition to re-label the PheWeb trait.
- GIGASTROKE EAS AIS is "any ischemic stroke", 19,032 discovery cases and 237,242 controls. TPMI Phecode 433.21 and GIGASTROKE AIS are related but **not the same operational case definition**, and neither is automatically equivalent to TOAST CES/LAS/SVS.
- The source \`OBS_CT=204099\` remains inconsistent with phenotype-page case+control total \`313,909\`; the source effective sample size is unapproved. Confirm model- and SNP-specific observation counts, covariates and Phecode inclusion/exclusion definitions before any LDSC, LDSC-SEG, MAGMA, MR or coloc claim.

### Next required gates
1. Select thesis-level scientifically matched AIS phenotype(s), compare case ascertainment and sources/overlap, and document the different outcome definitions.
2. Harmonize all available TPMI GRCh38 summary variants to verified EAS AIS GRCh37 (or validate an official GRCh38 outcome file), record allele/strand QC and overlap counts, without assuming that a rare lead absent from the EAS outcome is replicated.
3. Explicitly approve phenotype-specific and SNP-specific effective N / LDSC preconditions **only when supported by source data**, otherwise HOLD.

## Genomewide allele harmonisation and independent N_evidence gates — 2026-11-21

This stage replaces the scientifically unsupported "single rare lead as required independent replication" approach with a measured whole-genome SNP overlap audit. It **does not** automatically allow LDSC/MAGMA/MR/coloc.

1. **Independent sample-size evidence run**: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889438772 , **SUCCESS**, artifact \`14701617112\`. Official phenotype source and GWAS Catalog v2 metadata both verified: TPMI Phecode 433.21 \`9,249 cases / 304,660 controls = 313,909\` portal samples, but the pinned original PLINK2 \`OBS_CT = 204,099\`. **Difference 109,810**, case/control breakdown among model observations unknown. GIGASTROKE EAS AIS discovery \`19,032 cases / 237,242 controls\` are study-level metadata, not per-SNP N_eff. Metadata technical PASS; N_eff/MR/coloc/LDSC scientific **HOLD**. No N_eff imputed.
2. **Original genomewide intersection run**: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889274048 , initially started on the official two sources: TPMI GRCh38 original \`1,739,048\` associations (pinned compressed SHA-256 \`ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62\`) and GIGASTROKE EAS AIS \`GCST90104545\` GRCh37 \`6,725,216\` associations (independently recomputed source SHA-256 \`c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3\`). Both source hashes must match before any overlap statistics are accepted.
3. Use publisher UCSC \`hg38ToHg19.over.chain.gz\`, MD5 \`ff3031d93792f4cbb86af44055efd903\` checked against the UCSC MD5 manifest, and require the previously independently verified TPMI chr19 lead to map precisely \`GRCh38:19:15462568:T:C -> GRCh37:19:15573379:T:C\`. Map only unique autosomal point coordinates, reject multimap/transchrom and invalid alleles. Never equate source A1 with ALT without checking.
4. Only aggregate outcome effect alleles matching the mapped TPMI ref/alt alleles; reject unresolved A/T and C/G palindromic SNPs; record swapped effect sign and AF discrepancy buckets. Upload **aggregate-only**, not a joint per-SNP effects dataset. Outcome variants remain unsuitable for causal inference until case definitions, sample independence, correct N and full harmonized summary QC are approved.

**Completion criterion** for this stage: original source hashes, official chain checksum and expected known locus check all pass, whole-genome row counts match, actual coordinate and allele-matched overlap counts are observed in an Actions artifact. Code existence or workflow success alone is insufficient.

## Whole-genome TPMI–GIGASTROKE EAS AIS allele overlap: observed result (2026-11-21)

**Verified Actions execution:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889274048 (job \`136801311574\`, **SUCCESS**; aggregate artifact \`14701104342\`). The original-source download checks, UCSC hg38ToHg19 chain MD5 gate, independent known lead mapping and offline allelic-orientation tests all completed before aggregate results were emitted.

**Observed stdout, not projected coverage:**
- TPMI original GRCh38 associations: **1,739,048**.
- GIGASTROKE EAS AIS \`GCST90104545\` original GRCh37 associations: **6,725,216**.
- TPMI mapped unique GRCh37 genomic positions: **1,722,752**.
- Intersecting genomic-coordinate association rows: **1,507,267**.
- Non-palindromic direct allele-matched intersecting rows: **1,507,261**.
- **6** coordinate rows with allele mismatch were rejected.
- Frequency deviation bins across matched rows, comparing TPMI A1_FREQ with harmonized EAS effect allele frequency:
  - absolute delta <=0.05: **1,087,668**
  - >0.05 through <=0.10: **328,074**
  - >0.10 through <=0.20: **90,338**
  - >0.20: **1,181**.
  - Thus **419,593/1,507,261 (~27.84%)** have absolute delta >0.05 and **91,519/1,507,261 (~6.07%)** have absolute delta >0.10. These are QC indicators and not an a priori reason for genome-wide exclusion across different ancestry samples.
- Script technical status: \`PASS_ORIGINAL_BUILD_AND_ALLELE_AGGREGATE_QC\`; source error \`null\`; scientific status **HOLD**. These counts are public summary-statistics aggregate only; no per-SNP merged association file was published.

**Interpretation boundary:** The 1,507,261 figure is a *technical allele-matched intersection*, not an approved LDSC/MAGMA/MR/coloc analysis set. Effect-size/SE/p-value eligibility per SNP, ambiguous allele resolution, ancestry-specific LD reference choice, QC of high-frequency-discrepancy subsets, input coverage beyond the matched SNPs, and cohort overlap remain separate gates. One rare TPMI lead \`rs187759247\` is still **absent** at the exact GRCh37 locus in EAS AIS and cannot be called independently replicated.

**N and phenotype HOLD (previous independent Actions):** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889438772 confirmed TPMI portal n=313,909 (9,249 cases and 304,660 controls) versus original PLINK2 \`OBS_CT=204,099\`, a discrepancy of **109,810**. The 204,099 model observations have no proven case/control breakdown; no valid TPMI per-SNP binary \`N_eff\` has been established. TPMI Phecode \`433.21\` ("Cerebral artery occlusion, with cerebral infarction") is not operationally identical to GIGASTROKE's broad "Any ischemic stroke".

**Next scientific gates:**
1. Reconcile actual association-level N, effect-size/SE/p-value valid count, the high-AF-discrepancy subsets, and EAS LD reference suitability. The frequency comparisons should be checked on valid, unambiguous, matched SNPs and may reflect ancestry differences.
2. Verify sample-overlap independence and an approvable case/control N_eff from genuine source metadata; preserve LDSC/MAGMA/MR/coloc HOLD until all relevant gates pass.

## Genome-wide effect-statistic eligibility audit (additional verification)

- The independently verified aggregate baseline is **1,507,261 allele-matched autosomal SNP pairs**, from Actions [#44889274048](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889274048) attempt 1.
- The matcher was extended at research commit `4c918199229ba1fea39d7db6837f6fe80469fff8` to retain, in memory only, TPMI original OR, LOG(OR)_SE, P, and OBS_CT, plus EAS AIS original beta, SE, P and OR for numeric QC of **the same allele-aligned pairs**. Bad/NA/non-finite/negative SE or invalid P cannot pass the numeric gate; EAS reported OR-vs-beta discrepancy is counted. No per-variant effects are exported.
- Negative unit tests added in `6e97fffb9313a13d46ec062b2f9c8a776f503146`.
- The existing pinned-original workflow [#44889274048](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889274048) was directly rerun using the GitHub rerun-job API (attempt 2). **Do not report a new validated eligible-pair count before inspecting that attempt's terminal scientific audit log.**
- Any numerical PASS is a necessary but insufficient prerequisite for LDSC, MAGMA, MR, coloc. The relevant LD reference, QC on allele-frequency discrepancies, diagnosis equivalence and case/control observations remain unresolved.
- Cohort independence is not established merely because published cohort names differ. TPMI is a Taiwan-based clinical cohort (Nature 2025), while GIGASTROKE uses several EAS source cohorts (Nature 2022). The cohort contribution and participant-level duplication must be independently substantiated before a cross-dataset replication claim. **No independence PASS by assumption.**

## SNP-level numeric validity — confirmed second Actions attempt (2026-11-21)

- [GitHub Actions #44889274048](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44889274048) **run attempt 2 SUCCESS**, job `136805213014`; artifact `14701694336`. Negative/offline tests **14/14 PASS**.
- Whole-origin counts reproduced: TPMI `1,739,048` summary associations; EAS AIS `6,725,216`; unique mapped TPMI positions `1,722,752`; original coordinate intersection `1,507,267`; **allele-matched `1,507,261`** and **6 allele mismatches**.
- All `1,507,261` matched variants passed numerical validity checks for TPMI OR, SE, P, OBS_CT and EAS AIS beta, SE, P, reported OR. Log `NUMERIC_PASS_FINITE_RANGE=1507261` and `NUMERIC_ELIGIBLE_MATCHES=1507261`, with **no** recorded `EAS_LOGOR_BETA_DIFF_GT_0.005` among matched records. This is a **finite/range consistency gate**, not phenotype-independent association, heterogeneity, INFO/MAF, LDSC harmonization, or p-value model validation.
- Original EAS-vs-TPMI A1 frequency discrepancy: `419,593` matched SNPs >0.05; `91,519` >0.10; `1,181` >0.20. Different East Asian cohorts may genuinely differ in population frequencies, and any unexplained outliers must remain QC flags rather than be corrected or discarded without ancestry-aware evidence.
- TPMI effective sample size remains **not approved** (PheWeb `313,909` portal individuals vs PLINK2 `204,099` additive observations). Do not manufacture binary cases/controls among those observations.
- **Cohort overlap evidence**: TPMI 2025 Nature cohort study (https://www.nature.com/articles/s41586-025-09680-x) describes Academia Sinica with 16 Taiwan medical centres; GIGASTROKE 2022 Nature study (https://www.nature.com/articles/s41586-022-05165-3) combines ancestry-specific stroke studies. These source descriptions establish different public study descriptions, **not verified participant-level non-overlap**. It is therefore wrong to assert that TPMI and GIGASTROKE EAS AIS are proven independent from study names or publication years.
- Scientific decision: proceed only to source-level cohort-inclusion / participant-independence evidence review and phenotype-compatible preconditions. Keep formal independent replication, LDSC, MR, coloc at **HOLD** pending effective N and independence resolution.

## Cohort independence evidence gate — executed (2026-11-21)

- [Action #44893035771](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44893035771) **SUCCESS** (job `136812756194`, artifact `14703275181`). Its source-provenance gate explicitly requires a documented original-cohort inclusion list and custodian-validated participant nonoverlap before interpreting EAS AIS as an independent replication.
- Observed summary: `cohort_list_documented=false`, `custodian_proven_no_overlap=false`, `status=HOLD_INDEPENDENCE_NOT_PROVEN`, `LDSC_MR_coloc_enabled=false`.
- **Interpretation:** The absence of documented evidence is not proof of actual participant overlap, nor proof of no overlap. GIGASTROKE discovery/replication samples must remain distinct. No independent replication, LDSC, MR or colocalization is approved based solely on this gate.
- Previously verified numeric source QC: action `44889274048` attempt 2, 1,507,261 matched variants valid under finite/range checks; 419,593 AF difference >0.05; 91,519 >0.10; 1,181 >0.20. The TPMI effective N mismatch (313,909 portal vs 204,099 regression observations) remains unresolved.

## Primary full-text source audit: TPMI/TWB overlap versus GIGASTROKE discovery (2026-11-21)

**Actual CI evidence:** [GitHub Actions #44895009897](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44895009897), **SUCCESS**, job `136818676105`, artifact `14703820461`, source extraction code `scripts/ckm_stroke/audit_gigastroke_tpmi_cohort_literature.py`. Both primary papers were retrieved as Europe PMC full-text XML and parsed with explicit source excerpts.

1. [GIGASTROKE 2022, Mishra et al., Nature](https://www.nature.com/articles/s41586-022-05165-3) / [primary full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9524349/fullTextXML): study text explicitly says the East Asian iPGS **model was trained on BBJ** and was additionally **validated** using Taiwan Biobank (TWB) participants with **1,399 prevalent AIS cases and 86,283 controls**. These TWB participants **must not** be classified as part of the `GCST90104545` EAS AIS GWAS discovery set solely on this text.
2. [Taiwan early-onset ischemic stroke paper 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13239414/) / [primary full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13239414/fullTextXML) explicitly states that **participant overlap between TPMI and TWB is approximately 1%**. This is an approximate cohort-level statement, not the outcome-specific `433.21` participant overlap or a validated `GIGASTROKE AIS discovery` overlap count.
3. **Interpreting together:** there is an actual published approximate TPMI↔TWB overlap, and GIGASTROKE used TWB in **iPGS validation**. Therefore treating all GIGASTROKE-derived datasets as inherently independent of TPMI would be wrong. Conversely, the TWB PGS validation association does not establish overlap with the GIGASTROKE discovery-only `GCST90104545` summary statistics.
4. We have **not established a complete, study-specific GIGASTROKE EAS AIS discovery cohort list**, their individual overlap with TPMI, or a custodian-confirmed overlap count. Keep independent replication and cross-sample analyses that require independence at **HOLD**. Before a non-overlap claim, evaluate the original GIGASTROKE supplementary study-population methods, specifically whether TPMI or any shared contributing GWAS study belongs to the EAS AIS discovery sample rather than independent follow-up or PGS validation.

**Unresolved regression N:** TPMI PheWeb `433.21` reports **9,249 cases, 304,660 controls (313,909 total)** but the verified PLINK2 summary has `OBS_CT=204,099`. Per [PLINK2 source format documentation](https://www.cog-genomics.org/plink/2.0/formats), `OBS_CT` denotes samples in regression, *not* the number of cases, controls, or binary `N_eff`. The **109,810 difference** does not identify which samples were removed or the case/control distribution. No inferential case/control split was computed or approved. `LDSC/MR/coloc` remain **HOLD**.

## Original GIGASTROKE EAS AIS discovery cohorts and sample count VERIFIED — 2026-11-21

This source-backed finding supersedes earlier "complete EAS AIS discovery cohort roster not verified" statements, but does **not** establish subject-level independence.

- [Action #44899518851](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44899518851) SUCCESS: GIGASTROKE primary XML identifies original supplement 41586_2022_5165_MOESM4_ESM.xlsx (Supplementary Tables 1–56).
- [Action #44899624455](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44899624455) SUCCESS, attempt 2: downloaded original [Springer publisher XLSX](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05165-3/MediaObjects/41586_2022_5165_MOESM4_ESM.xlsx), verified that Supplementary Table 1, EAS rows 50–65, columns K (AIS cases) and P (shared controls) contain the actual cohort contributions.
- [Action #44899774706](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44899774706) SUCCESS: independent CI totals matched original Table 1 and GWAS Catalog. Artifact ID 14705337599. Scientific HOLD for individual overlap and N_eff.

| EAS AIS discovery cohort | AIS cases (K) | Shared controls (P) |
|---|---:|---:|
| BIOBANK_JAPAN | 17,651 | 159,610 |
| HISAYAMA_FSR | 1,113 | 901 |
| INTERSTROKE_ASN | 222 | 698 |
| KCPS2_ASN | 46 | 12,366 |
| **Case-contributing subtotal** | **19,032** | **173,575** |
| CKBB_ASN ten regional study cohorts (AIS cases NA) | — | 63,667 |
| **Verified EAS AIS discovery total** | **19,032** | **237,242** |

EAS row 60 HISAYAMA_ASN is incident-only and explicitly not in the global combined incident/prevalent meta-analysis; do not include its cases in EAS AIS discovery. CKBB regional cohorts provide shared controls and any-stroke cases, but their AIS case entries are NA, **not zero-risk associations**. Original source row TOTAL EAS (row 65) confirms exact totals.

**No TPMI or Taiwan Biobank (TWB) study appears within this original EAS AIS discovery cohort roster.** The GIGASTROKE primary paper's use of TWB for iPGS validation is distinct from the discovery GWAS. The 2026 Taiwan paper's approximate 1% TPMI–TWB participant overlap **does not imply** that the GIGASTROKE discovery sample uses TWB or that TPMI participants overlap with BBJ/Hisayama/INTERSTROKE/KCPS2.

**Two-level independence outcome:** Source-cohort roster PASS with no named TPMI/TWB discovery contribution; participant-level nonoverlap remains **HOLD** pending individual overlap/custodian evidence. Do not claim verified zero overlap.

**Other unresolved scientific gates:** TPMI Phecode 433.21 portal cases 9,249 + controls 304,660 = 313,909, but PLINK2 association OBS_CT 204,099 (109,810 difference). No source-supported case/control regression N_eff. Operational Phecode "Cerebral artery occlusion, with cerebral infarction" is related to, but not identical to, GIGASTROKE broad AIS. Keep LDSC/MR/coloc at HOLD.

Earlier independent allele/quantitative summary QC PASS: 1,507,261 SNP pairs. This remains a technical overlap, not an approved causal or independent replication result.

## Primary GIGASTROKE AIS cohort gate rerun with recovered source table — 2026-11-21

- Scientific cohort gate [Action #44893035771](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44893035771) attempt 2 **SUCCESS**, job 136833365591; output artifact 14706045108.
- Unlike prior attempt 1 (which lacked a cohort roster), the latest gate **re-parses the actual publisher Supplementary Table 1**, imports its 14 EAS discovery cohort rows and verifies the reported 19,032 AIS cases / 237,242 shared controls before judging the independence evidence.
- **Observed:** cohort_list_documented=true, custodian_proven_no_overlap=false, status=HOLD_INDEPENDENCE_NOT_PROVEN, LDSC_MR_coloc_enabled=false.
- Source roster does not list TPMI or TWB in AIS discovery; that is evidence against direct inclusion of their named studies, **not** verified zero participant overlap. TPMI observed regression cases/controls remain unavailable, so effective N is not approved.

## TPMI Phecode 433.21: original SAIGE vs PLINK2 model-specific N resolved in part (2026-11-21)

**New source-verified discovery: The TPMI portal case/control counts belong to its SAIGE source in a genuine sampled portion; they are NOT evidence of the PLINK2 regression case/control breakdown.**

1. [Action #44906623458](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44906623458) SUCCESS; primary literature (Europe PMC XML) was directly parsed for 2025 publications [Chen et al., Nature 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12675292/) and [Yang et al., Nature 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12675286/). Chen et al. explicitly describes GWAS sample n=363,447 using SAIGE mixed-effect regression, an unrelated GWAS subset n=248,754 using PLINK2, plus 1:10 age-sex matching for imbalanced phenotypes with case/control ratio below 1:20. **This makes model-specific sampling a source-supported explanation category, not a demonstrated accounting for the precise difference of 109,810.**
2. [Action #44906755519](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44906755519) SUCCESS; inspected bounded initial 1MiB compressed bytes from *both original PheWeb phenotype 433.21 downloads*:
   - `https://pheweb.ibms.sinica.edu.tw/download/gwas_plink2/433.21` -> official PLINK2 `433.21.glm.logistic.gz`, source length **58,302,393 bytes**, 17-column header including `OBS_CT` and **no N_case/N_ctrl columns**; earlier full-original verified `OBS_CT=204,099`.
   - `https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21` -> official SAIGE `433.21.saige.txt.gz`, source length **816,992,833 bytes**, 20-column header including **`AF_case`, `AF_ctrl`, `N_case`, `N_ctrl`, `QC`**. This is a *distinct* summary-statistics analysis.
3. [Action #44906864247](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44906864247) SUCCESS; verified HTTP 206 for first 1MiB of SAIGE original gzip and extracted exactly **425 complete variant records**. Across the **425/425** valid rows, `N_case=9,249`, `N_ctrl=304,660`, sum **313,909**, identical to source PheWeb phenotype portal counts. Distinct N_case/N_ctrl pair in this **sampled prefix: one**. **This is not full-SAIGE-file N validation**; the remaining original GWAS rows might differ, and full gzip CRC was not examined.
4. **Scientific gate interpretation:** SAIGE has a directly documented *per-variant* case/control N field in its original summary statistics, and its initial 425 rows match portal N exactly. PLINK2 has only `OBS_CT=204,099`; the actual PLINK2 regression case/control split remains **UNVERIFIED**. Never copy SAIGE case/control N onto PLINK2 data. Do not use PLINK2 `OBS_CT=204,099` alone as effective N; do not claim the exact exclusion/matching mechanism is proven. If SAIGE is the desired GWAS for LDSC, its entire source, per-variant N cases and controls, allele QC, build/LD reference, and effect statistics must be verified separately.
5. **No downstream scientific GO:** TPMI cohort-level independent AIS replication, LDSC, MR, coloc remain HOLD pending model-specific N and assumptions. The previously validated **1,507,261 TPMI-PLINK2 x GIGASTROKE EAS AIS SNP pairs** must **not** silently be re-labelled as SAIGE pairs.

Progress: original 109,810 discrepancy is now explicitly classified as **between two different GWAS analysis sample bases**, not corrected by assuming that SAIGE N equals PLINK2 N.

## SAIGE full original source case/control audit initiated — November 21, 2026

- The primary TPMI PheWeb source clearly exposes *two different* GWAS methods: PLINK2 `/download/gwas_plink2/433.21` (58,302,393 compressed bytes, `OBS_CT` only), and SAIGE `/download/saige_gwas/433.21` (816,992,833 compressed bytes, `N_case` and `N_ctrl` per GWAS association).
- Prior source-bounded [Action #44906864247](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44906864247) read 425/425 first source SAIGE SNPs at `N_case=9249,N_ctrl=304660`, matching PheWeb phenotype page 313,909 total. **It did not verify all SAIGE original SNP rows or its gzip trailer.**
- [Action #44907119885](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44907119885), [source commit](https://github.com/kimtk94/IS_Analysis_V3/commit/2066520360c5230b586f70b5616931537ddeba92) and [workflow commit](https://github.com/kimtk94/IS_Analysis_V3/commit/72292e779d9465007f7d256eed65aba9f8717e4a): stream entire 816,992,833-byte official SAIGE source, verify full gzip integrity and downloaded size, report source SHA256, count all variants, N_case / N_ctrl pair distribution, missing N cases, beta/SE/p validity, and per-chromosome counts. Output **aggregate JSON only**, original per-SNP effects are not uploaded. **Do not mark full-source completion until the terminal run/artefact is checked.**
- The prior `PLINK2` study GWAS set used for the 1,507,261-variant harmonization is a different analysis sample and different summary statistical model. Do **not** reuse SAIGE sample N, effect estimates, or QC results as if they were PLINK2.
- A *future independent SAIGE-vs-GIGASTROKE* GRCh38↔GRCh37 chromosome liftover and effect-allele matching run requires separate original SAIGE source hash and correct coordinate build verification. EAS AIS original `GCST90104545` QC is already available but does not validate SAIGE harmonization by itself.
- Even if all SAIGE N_case/N_ctrl values are present, phenotype definitional non-equivalence (TPMI Phecode 433.21 vs broad AIS), participant-level independence, and LD reference/assumption checks remain. **LDSC/MR/coloc scientific HOLD** until model-specific preconditions are met.

## Independent SAIGE x GIGASTROKE EAS AIS original-format preflight (2026-11-21)

- [Action #44907344742](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44907344742) **SUCCESS**; job `136855638152`, artifact `14708412035`. The original TPMI `433.21` SAIGE and `GCST90104545` EAS AIS gzip headers and provider names were independently checked with strict bounded HTTP range requests. Output `PASS_BOTH_ORIGINAL_METHOD_SCHEMAS`.
- Official source SAIGE header contains `CHR POS MarkerID Allele1 Allele2 ... BETA SE ... N_case N_ctrl QC`; EAS AIS contains `chromosome base_pair_location effect_allele_frequency beta standard_error p_value odds_ratio ci_lower ci_upper effect_allele other_allele`. Case/control N is explicitly part of **SAIGE**, not the separate PLINK2 summary.
- This is an **original source/schema feasibility PASS**, not a matched SNP analysis. `harmonization_executed=false`, `SAIGE_build_verified=false` in this run, `SAIGE_allele2_effect_orientation_verified=false` and `scientific_status=HOLD`. The previous matched TPMI PLINK2 EAS AIS count (1,507,261) is not a SAIGE matched count.
- Full SAIGE per-SNP N case/control audit [Action #44907119885](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44907119885) was already launched using the actual 816,992,833-byte source, but must be separately checked for terminal success and full-source result before any claim about N_case/N_ctrl consistency or subsequent harmonization readiness.
