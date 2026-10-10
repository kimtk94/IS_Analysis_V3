# IS external East Asian replication readiness — TPMI / CKB
Date: 2026-10-11 KST. **Reference preflight only, not GWAS replication.**

## Motivation and exposure/outcome independence

Koyanagi 2024 Japanese alcohol GWAS and Japanese BBJ IS share a cohort-origin risk. GIGASTROKE EAS AIS includes BBJ, and therefore must not be labelled independently replicated. The currently available European-ancestry AIS GWAS is cross-ancestry context, not independent Japanese/EAS replication. Retain **ALDH2 BLOCKED** and **ADH1B EXPLORATORY**.

## TPMI independent East Asian outcome candidate

- Public PheWeb: https://pheweb.ibms.sinica.edu.tw/
- **Phecode 433.21 — Cerebral artery occlusion, with cerebral infarction**, with reported **9,249 cases and 304,660 controls** (based on local HTML probe snapshot).
- **GRCh38** (verified from locally saved TPMI PheWeb source JavaScript `window.model.hg_build_number=38`).
- Source HTML advertises a direct `Summary Statistics` download at `/download/saige_gwas/433.21`. Earlier server probe returned HTTP 200, later 429; a **single non-retried HEAD on 2026-10-11 KST returned HTTP 403**. No rate-limit workaround attempted. **No TPMI GWAS statistics downloaded or analyzed yet**.
- 433.21 is a fairly close cerebral infarction phenotype, but not necessarily identical in adjudication to BBJ/GIGASTROKE AIS. Evaluate endpoint semantics before quantitative meta-analysis.
- A public independent Taiwan dataset does not establish nonoverlap with BBJ for every participant without checking cohort details, but there is no automatic BBJ reuse assumption.

## GRCh37 → GRCh38 coordinate and reference audit

A first Ensembl Variation + Sequence REST pass verified 3/11, with HTTP 500 errors on 8 requests. It did **not** imply biological discordance. A second independent workflow used:
1. existing **UCSC hg19ToHg38 chain** at `/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz` (locally vendored pyliftover);
2. UCSC public **hg38 reference sequence** to establish the actual REF nucleotide;
3. Ensembl GRCh38 variation rsID mapping to cross-check identity and original allele set.

All **11/11 coordinates and REF/ALT pairs were verified with local chain plus UCSC reference**; **9/11 Ensembl rsID cross-checks completed successfully in this second execution**. rs671 and rs1229984 independently showed Ensembl HTTP errors during the second pass, despite rs671 having passed the first Ensembl sequence/variation attempt and rs1229984 mapping having been observed in a separate ad hoc variation API response. **Do not report 11/11 Ensembl rsID verification in one unified execution**. Follow-up can retry these 2 (respectfully) or use the chain+sequence verified map.

| Gene | rsID | GRCh37 REF:ALT | GRCh38 REF:ALT | UCSC/chain |
|---|---|---|---|---|
| ALDH2 | rs184590119 | 12:111900807:T:C | 12:111463003:T:C | PASS |
| ALDH2 | rs671 | 12:112241766:G:A | 12:111803962:G:A | PASS |
| ADH1B | rs1229984 | 4:100239319:T:C | 4:99318162:T:C | PASS |
| ADH1B | rs34144181 | 4:100293600:C:T | 4:99372443:C:T | PASS |
| ADH1B | rs12502498 | 4:100297551:T:C | 4:99376394:T:C | PASS |
| ADH1B | rs9997653 | 4:100299453:C:A | 4:99378296:C:A | PASS |
| ADH1B | rs1826906 | 4:100301048:T:C | 4:99379891:T:C | PASS |
| ADH1B | rs1442485 | 4:100306411:C:T | 4:99385254:C:T | PASS |
| ADH1B | rs12502290 | 4:100316666:G:A | 4:99395509:G:A | PASS |
| ADH1B | rs2584461 | 4:100321573:T:C | 4:99400416:T:C | PASS |
| ADH1B | rs17028965 | 4:100322106:C:T | 4:99400949:C:T | PASS |

Artifacts (under existing experimental `susie_rss_finite_ref_sandbox_v1`):
- `tpmi_grch38_mapping_v1/TPMI_11_SNPS_GRCH38_MAPPING.tsv` and manifest: first-pass Ensembl, **3/11**.
- **`tpmi_grch38_crossbuild_v1/TPMI_11_SNP_GRCH38_CHAIN_UCSC_ENS_AUDIT.tsv` and `TPMI_GRCH38_CROSSBUILD_AUDIT_MANIFEST.json`**: reference-validated **11/11**, Ensembl same-pass **9/11**, 2 transient API errors.

Source code and offline regression tests:
- `scripts/is/build_tpmi_grch38_mapping.py`, `tests/test_tpmi_grch38_mapping.py` (5 tests).
- `scripts/is/audit_tpmi_grch38_chain_ucsc.py`, `tests/test_tpmi_grch38_chain_ucsc.py` (3 tests).
- All **8 TPMI reference mapping regression tests PASSED locally on server**.

## China Kadoorie Biobank (CKB): access blocker is encryption, not absence of GWAS

- Source: https://pheweb.ckbiobank.org/pheno/i63
- CKB `i63`: locally downloaded encrypted ZIP `/srv/is-analysis/data/is/east_asia/china/ckb/CKB_i63_IS.zip` (169,496,188 bytes) contains flagged encrypted `i63.tsv` (uncompressed 544,277,491 bytes). No decryption password supplied.
- CKB PheWeb i63 currently describes **14,302 cases / 67,954 controls**; variation effects for this PheWeb are for reference/alternative alleles in GRCh38 (PheWeb indicates shift to GRCh38 in March 2025), confirm file schema after authorized decryption.
- Official download-key instructions: https://pheweb.ckbiobank.org/about — email `ckbaccess@ndph.ox.ac.uk` with subject `CKB summary statistics request` and researcher **name, institution, institutional address and institutional email**. Keys are not to be redistributed; decrypted results may only be shared with the immediate research group under their policy.
- **Important administrative distinction:** Public PheWeb *summary-statistics download key* request is distinct from more comprehensive **individual-level CKB data access**, which can require formal agreement and fees. Do not conflate these processes.
- The publicly documented CKB i63 cohort is a potential non-BBJ Chinese IS GWAS candidate; confirm exact phenotype, ancestry structure, variants, sample-overlap provenance and study design before calling it an independent replication.

## Research guardrails

- Do **not** infer zero effect from unreported variants such as ALDH2 rs671 missing in EUR AIS.
- Do **not** infer mediation from a significant ALDH2 alcohol association and significant BBJ/EAS AIS association.
- BBJ vs GIGASTROKE EAS cannot be treated as two independent replications without leave-BBJ-out/cohort-wise statistics.
- Fine-mapping credibility remains constrained by 1000G proxy LD and SuSiE RSS reliability warnings. Singleton CS PIP near 1 is not causal validation.
- Keep the broad IS candidate universe **2,225 genes** without narrowing due to this side analysis.

## Next steps

1. Once TPMI permits access, download phenotype 433.21 summary statistics using normal permitted access; verify phenotype, summary-stat effect allele, genome build and exact 11 mapped variants. Keep original file untouched and output a separate harmonized/verified table.
2. Request CKB PheWeb summary-stat decryption key through official published process (user-authorized institutional contact); confirm original file SHA256 and format after decryption.
3. Compare signs/p-values with BBJ/GIGASTROKE cautiously, test cohort independence and instrument pleiotropy, and maintain `BLOCKED/EXPLORATORY` states until matched-LD conditional reliability is resolved.
