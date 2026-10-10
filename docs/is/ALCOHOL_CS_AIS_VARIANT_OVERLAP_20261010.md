# Alcohol fine-mapped variants vs. AIS GWAS — audited evidence (2026-10-10)

## Scientific status

**DESCRIPTIVE ASSOCIATION LOOKUP ONLY. Not MR, not colocalization, not mediation, not independent replication.**

- **Official genome assembly:** Koyanagi 2024 J-CGE unstratified daily alcohol intake GWAS (`n=154,570`) explicitly lists **hg19 (GRCh37)** on the Japanese National Cancer Center / J-CGE distribution page: https://epi.ncc.go.jp/cgi-bin/cms/public/index.cgi/jcge/en/download/index
- GIGASTROKE EAS AIS accession GCST90104545 metadata lists **GRCh37** and `n=256,274`; local source metadata: `/srv/is-analysis/data/is/reference/gigastroke/metadata_phase8/GCST90104545-meta.yaml`.
- The original GIGASTROKE AIS source file checksum **MD5=904e12e33ec3c2eeefe4dc93bec4995c** matches that official local source metadata.
- 11/11 deduplicated SNPs independently matched the local GRCh37 FASTA REF; VEP GRCh37 colocated rsID plus coordinate-and-allele-set mapping was 11/11. The rsID audit is not independent dbSNP validation.
- **Cohort reuse risk confirmed:** Koyanagi 2024 includes 134,993 BioBank Japan individuals among six Japanese cohorts (https://pmc.ncbi.nlm.nih.gov/articles/PMC10816704/). GIGASTROKE explicitly includes BBJ in EAS GWAS (https://www.nature.com/articles/s41586-022-05165-3). BBJ ischemic-stroke GWAS and GIGASTROKE EAS therefore cannot be called independent replication; exact per-participant overlap remains unknown.
- 1000G EAS LD (504 individuals) is a population proxy, **not cohort-matched**. SuSiE RSS R reliability/sensitivity warnings persist. **ALDH2 fine-mapping remains BLOCKED; ADH1B remains EXPLORATORY**. GWAS statistical significance does not promote fine-mapping status.

## REF/ALT canonical identifier finding

The GWAS Catalog source of GCST90104545 provides **effect_allele** and **other_allele**, not declared FASTA REF/ALT fields. The historical canonical `variant_id` uses `chr:pos:other_allele:effect_allele` in observed records. That ID may equal GRCh37 `CHR:POS:REF:ALT` *only by coincidence* when `other_allele` is the reference.

- At 4:100239319 (rs1229984), GRCh37 FASTA `T:C`; historical GIGASTROKE canonical ID `C:T`. The original GWAS effect allele T was the FASTA REF. The source beta was -0.0124; properly aligned beta for ALT=C becomes **+0.0124**.
- Four of eleven queried variants have such reversed historical identifiers; they are now validated via **raw source + metadata MD5 + EA/OA + GRCh37 FASTA + beta/SE/p/EAF** in a separate derivative.
- Across a deterministic genome-wide thinned sample of **1,345 reference-resolved SNVs**, **583 (43.3%)** had historical canonical ID's first allele unequal to GRCh37 FASTA REF; **762** matched. **All 1,345 sampled IDs** followed `other:effect` ordering. This is a sampled QC rate, **not** a verified full-genome fraction.
- Source and canonical summary-statistical fields beta, SE, p, EAF, EA, OA all agreed for the ten observed query variants. Thus this is an **identifier semantics issue**, not evidence of source effect-statistic corruption.
- **No canonical source files were overwritten.** All new artifacts are in segregated experimental output folders. The helper `scripts/is/gigastroke_fasta_variant_contract.py` enforces FASTA-anchored REF/ALT and rejects unsafe biallelic/reference assumptions.

## Source-audited AIS p values (11 original alcohol credible-set variants)

| Gene | rsID | BBJ Japanese IS p | GIGASTROKE EAS AIS p | GIGASTROKE REF:ALT status |
|---|---|---:|---:|---|
| ALDH2 | rs671 | 1.99193648e-18 | 8.426e-18 | Source/FASTA exact |
| ALDH2 | rs184590119 | MISSING | MISSING | Missing |
| ADH1B | rs1229984 | 0.517197 | 0.3659 | **Repaired in derivative**, source/FASTA-verified |
| ADH1B | rs34144181 | 0.754749 | 0.5823 | Source/FASTA exact |
| ADH1B | rs12502498 | 0.577723 | 0.1665 | **Repaired in derivative**, source/FASTA-verified |
| ADH1B | rs9997653 | 0.586082 | 0.1656 | Source/FASTA exact |
| ADH1B | rs1826906 | 0.591041 | 0.1618 | **Repaired in derivative**, source/FASTA-verified |
| ADH1B | rs1442485 | 0.612576 | 0.1721 | Source/FASTA exact |
| ADH1B | rs12502290 | 0.587270 | 0.1850 | Source/FASTA exact |
| ADH1B | rs2584461 | 0.897024 | 0.5094 | **Repaired in derivative**, source/FASTA-verified |
| ADH1B | rs17028965 | 0.968565 | 0.5746 | Source/FASTA exact |

**22 total phenotype–variant pairs:** BBJ exact 10, BBJ missing 1; GIGASTROKE source+FASTA exact 6, canonical ID corrected in derivative 4, GIGASTROKE missing 1. No single ADH1B variant tested shows genome-wide significance for IS/AIS here. This does not exclude alternative regional variants, distinct phenotypes, mediation or pleiotropy; it simply describes these comparisons.

### Shared ALT beta

- ALDH2 rs671 `12:112241766:G:A`: alcohol ALT(A) beta **-1.3215** (log2[g/day+1]), BBJ IS beta **-0.1097493721** (p=1.99e-18), GIGASTROKE EAS AIS beta **-0.1514** (p=8.426e-18). Association directions concordant; **not independent replication**, no mediation inference.
- ADH1B rs1229984 `4:100239319:T:C`: alcohol ALT(C) beta **+0.1680**; BBJ IS ALT beta **+0.0083996654** (p=0.517), GIGASTROKE EAS ALT beta **+0.0124** (p=0.3659); GIGA canonical REF/ALT order repaired as above.
- ADH1B rs2584461 `4:100321573:T:C`: alcohol ALT(C) beta **-0.1048**, BBJ IS ALT beta **+0.00216995** (p=0.897), GIGASTROKE ALT beta **+0.0116** (p=0.5094). Neither outcome association is significant; avoid direction interpretation.
- Source alcohol p=0 for extreme ALDH2 signals is numerical underflow/representation, not a true p=0. Preserve original beta and SE.

## Executed artifacts (experimentally segregated)

Output root:
`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/susie_rss_finite_ref_sandbox_v1/`

- `ais_overlap_v2/ALCOHOL_CS_AIS_GWAS_DIRECT_OVERLAP.tsv` — first-pass BBJ/GIGASTROKE comparison.
- `gigastroke_refalt_provenance_v1/GIGASTROKE_AIS_11_SNP_FASTA_PROVENANCE.tsv` — raw/source/FASTA repair audit, source MD5 checked.
- `gigastroke_refalt_provenance_v1/GIGASTROKE_AIS_REFALT_AUDIT_MANIFEST.json` — 1,345-reference-SNV sampled QC and provenance.
- **`ais_verified_association_v1/ALCOHOL_CS_AIS_VERIFIED_ASSOCIATION_MATRIX.tsv`** — preferred corrected 22-row phenotype–SNP outcome table.
- `ais_verified_association_v1/ALCOHOL_CS_AIS_VERIFIED_MATRIX_MANIFEST.json` — diagnostic status + known BBJ sample-overlap limitations.

Code:
`scripts/is/audit_gigastroke_ais_refalt_provenance.py`,
`scripts/is/gigastroke_fasta_variant_contract.py`,
`scripts/is/build_alcohol_ais_verified_variant_matrix.py`.
All computations are read-only on input/canonical material; effects are oriented to FASTA-validated ALT.

## Required next work

1. Audit **all** downstream GIGASTROKE analyses for incorrect interpretation of historical `other:effect` ID as FASTA `REF:ALT`; only reprocess scientifically affected **derived** outputs, do not overwrite originals.
2. Seek truly independent non-BBJ AIS results, leave-BBJ-out GIGASTROKE summary estimates, or cohort-wise analysis where obtainable; ascertain overlap counts where possible.
3. Validate any new MR instruments against correlated pleiotropy, instrument independence and sample overlap before mediation claims; require suitable LD and locus-level conditioning/coloc for causal fine-mapping statements.
4. Keep ALDH2 **BLOCKED**, ADH1B **EXPLORATORY** until reliability warnings and conditional/LD issues resolved.
