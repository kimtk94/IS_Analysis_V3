# IS research — independent-cohort constraints, EUR AIS transferability and ALDH2 conditional LD audit

Date: 2026-10-10 (KST). Status: **source-audited exploratory associations only; no mediation, MR, fine-mapping, or independent EAS replication established.**

## Question

Do ADH1B/ALDH2 SNPs from Koyanagi 2024 alcohol-consumption credible sets show ischemic-stroke association outside an EAS data set sharing BBJ subjects? Can BBJ ischemic-stroke conditioning on ALDH2 rs671 establish a genuinely independent second signal?

## Source/case definitions and sample-overlap risk

- **Koyanagi 2024:** unstratified Japanese daily alcohol intake, n=154,570; Japanese population meta-analysis includes BBJ among constituent cohorts. Official J-CGE GWAS download catalogue declares **hg19/GRCh37**; study [PMC10816704](https://pmc.ncbi.nlm.nih.gov/articles/PMC10816704/), source [J-CGE data page](https://epi.ncc.go.jp/cgi-bin/cms/public/index.cgi/jcge/en/download/index).
- **BBJ ischemic stroke:** Japanese outcome; strong sample-origin overlap risk with Japanese alcohol-exposure GWAS.
- **GIGASTROKE EAS AIS (GCST90104545):** n=256,274; includes constituent BBJ cohorts. Cannot treat this as an independent replication of BBJ IS. Source metadata \`data/is/reference/gigastroke/metadata_phase8/GCST90104545-meta.yaml\`.
- **GIGASTROKE EUR AIS (GCST90104540):** European ancestry (as declared by local metadata), total n=1,296,908. Appropriate for a **cross-ancestry transferability check**, not an independent same-ancestry Japanese or EAS replication. It may have very different LD and MAF distributions, and its summary statistics do not document individual participant overlap.
- **CKB i63 IS:** downloaded \`/srv/is-analysis/data/is/east_asia/china/ckb/CKB_i63_IS.zip\`, 169,496,188 bytes; one \`i63.tsv\` file with ZIP encryption flag enabled. **No plaintext summary statistics or decryption password available**; do not claim CKB replication.
- **TPMI:** inspected public phenotype/probe files; no readily accessible AIS GWAS summary statistics in current server inventory.

## EUR AIS direct study: 11 curated alcohol credible-set SNPs

- Input: \`/srv/is-analysis/data/is/processed/gigastroke/broad_v1/GCST90104540_AIS_EUR_GRCh37.allele_pairs.tsv.gz\`
- The existing QC records the original-source MD5 \`1bfe4ae8a5042fb24cb0a562b05b0d2f\`, n=1,296,908, GRCh37 outcome; this is **recorded upstream QC**, not a fresh original EUR file re-hash in this analysis.
- Important: **EUR \`variant_pair_id\` may be lexicographically ordered allele pair** (observed nine matched records); it is **not REF:ALT**. The EAS historical canonical ID often uses \`other:effect\`, which is a different naming rule. Never transfer either ID naming convention to the other without separate reference checks.
- For every candidate the reference allele was fetched from the local GRCh37 FASTA and checked. Outcome betas and effect allele frequencies were oriented to the same verified reference ALT as the Japanese alcohol GWAS. Strict matching uses chromosome + position + allele set; no silent complement or harmonization by rsID alone.
- Nine of 11 SNPs were observed in EUR, **all 9 ADH1B**. Two ALDH2 SNPs (**rs671 and rs184590119**) were **not reported** in this EUR input.
- None of the nine EUR associations is genome-wide significant. The smallest EUR p-value among these candidates is **0.07948** at ADH1B **rs12502498**.
- **rs1229984:** Japanese alcohol \`beta_ALT=+0.168\`, Japanese alcohol \`ALT_EAF=0.2437\`; GIGASTROKE EUR AIS \`beta_ALT=-0.0102\`, \`p=0.7592\`, EUR \`ALT_EAF=0.9567\` (MAF 0.0433). The extreme frequency contrast underscores ancestry dependence; direction at nonsignificant AIS p values is **not reliable causal evidence**.
- **rs2584461:** EUR AIS \`beta_ALT=+0.0050\`, \`p=0.6471\`, \`ALT_EAF=0.8776\`; note low statistical strength.
- \`rs671\` is missing in this reported EUR source; **missing does not mean no effect**. Do not use EUR as ALDH2 replication.
- Nine tested EUR variants have reported MAF >=0.01. This does not imply high statistical power for an AIS effect or a comparable alcohol exposure model.

| Gene / rsID | Japanese BBJ IS p | GIGASTROKE EAS AIS p | GIGASTROKE EUR AIS p | EUR status |
|---|---:|---:|---:|---|
| ALDH2 rs671 | 1.99193648e-18 | 8.426e-18 | N/A | Not reported |
| ALDH2 rs184590119 | N/A | N/A | N/A | Not reported |
| ADH1B rs1229984 | 0.517197 | 0.3659 | 0.7592 | Matched, descriptive only |
| ADH1B rs34144181 | 0.754749 | 0.5823 | 0.5130 | Matched |
| ADH1B rs12502498 | 0.577723 | 0.1665 | 0.07948 | Matched |
| ADH1B rs9997653 | 0.586082 | 0.1656 | 0.09120 | Matched |
| ADH1B rs1826906 | 0.591041 | 0.1618 | 0.09907 | Matched |
| ADH1B rs1442485 | 0.612576 | 0.1721 | 0.08685 | Matched |
| ADH1B rs12502290 | 0.587270 | 0.1850 | 0.1017 | Matched |
| ADH1B rs2584461 | 0.897024 | 0.5094 | 0.6471 | Matched |
| ADH1B rs17028965 | 0.968565 | 0.5746 | 0.7112 | Matched |

**Interpretation:**
1. rs671 is convincingly associated with AIS in the BBJ/EAS analyses examined, **not independently replicated here**, and this does not demonstrate that alcohol mediates the association.
2. ADH1B alcohol instruments show strong exposure associations, but these *11 candidate SNP lookups* do not show genome-wide AIS significance across the examined outcomes. This is not a locus-wide negative result and cannot rule out other independent signals, low power, nonlinearity, interaction, or distinct disease definitions.
3. No causal effect ratio or MR estimate was computed; sample-overlap, horizontal pleiotropy and ancestry-specific LD assumptions remain unresolved.

## ALDH2 rs671 BBJ ischemic-stroke conditional proxy-LD sensitivity

**Important trait distinction:** The Phase8F conditional association table is **BBJ IS GWAS**, not the Koyanagi alcohol GWAS. This table cannot validate a second alcohol-exposure SuSiE CS.

Input: \`/srv/is-analysis/results/is/stage4_cross_eas/phase8f_l003_jpt_ld/L003_JPT_CONDITIONAL_VARIANTS.tsv\`, 3,033 SNPs.

- **JPT104 proxy LD:** 0 genome-wide significant residual variants after conditioning on rs671; minimum p \`1.37932409e-4\` at rs3782886, \`r²=0.946227\` with rs671.
- **EAS504 proxy LD:** 2 genome-wide significant residuals, **both very high LD with rs671**:
  - rs78069066 \`p=2.45783170e-15\`, \`r²=0.9871548\`;
  - rs3782886 \`p=8.36729717e-14\`, \`r²=0.9806765\`.
- At \`r² <=0.2\`, **no residual genome-wide significant variants** from either proxy panel. Threshold \`0.2\` is a predeclared in-house diagnostic heuristic, **not a universally sufficient test of independence**.
- This is external-LD **single-lead approximate conditioning**, not adequate cohort-matched Japanese conditional fine-mapping. It cannot prove zero secondary causal signals.
- The *different* Koyanagi alcohol SuSiE \`FINITE_REF_504\` and \`FINITE_REF_504_EB_MISMATCH\` results generated 2 singleton credible sets for ALDH2, but **all four total finite-reference SuSiE models emitted LD reliability and sensitivity warnings**; the ALDH2 alcohol GWAS/LD mismatch summary \`s=0.0908385487\` is high under this pipeline's in-house diagnostic gate. Therefore two apparent alcohol CSs must not be interpreted as two independent causal effects.

### Scientific gates

- ALDH2 alcohol fine-mapping: **BLOCKED**.
- ADH1B alcohol fine-mapping: **EXPLORATORY**.
- Validated independent ALDH2 secondary AIS signal: **NONE established** (not proof of absence).
- Validated causal link alcohol consumption → AIS: **NONE established**.
- Independent Japanese/EAS AIS cohort replication: **NOT AVAILABLE** from presently readable local data.
- Preserve the **2,225-gene original broad IS candidate universe**; this side-analysis does not shrink it.

## Generated files

Output root: \`/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/susie_rss_finite_ref_sandbox_v1/\`

- \`ais_eur_transfer_v1/ALCOHOL_CS_GIGASTROKE_EUR_AIS_TRANSFER.tsv\`
- \`ais_eur_transfer_v1/ALCOHOL_CS_GIGASTROKE_EUR_TRANSFER_MANIFEST.json\`
- \`aldh2_conditional_science_gate_v1/ALDH2_BBJ_IS_AND_ALCOHOL_CONDITIONAL_SCIENCE_GATE.json\`

Code:
- \`scripts/is/audit_alcohol_ais_eur_transfer.py\` — EUR source allele-pair and reference ALT QC.
- \`scripts/is/audit_aldh2_conditional_ld_gate.py\` — external-LD conditional residual and scientific gates.
- \`tests/test_alcohol_ais_eur_transfer.py\`, \`tests/test_aldh2_conditional_ld_gate.py\` — offline regressions.
- \`.github/workflows/is-alcohol-ais-science-gate-ci.yml\` — checks six allele/association/provenance scripts; CI workflow created, **GitHub Actions run result must be checked separately**.
- Previous: \`docs/is/ALCOHOL_CS_AIS_VARIANT_OVERLAP_20261010.md\`.

## Prioritized unblock conditions

1. Obtain a genuinely non-BBJ **Japanese or broader EAS AIS** association resource where individual cohort provenance and sample overlap can be assessed; investigate CKB ZIP authorized decryption/access and alternative individual-cohort public GWAS.
2. Obtain adequate **cohort-matched Japanese genotype LD** for rs671/ALDH2 locus (or matched in-sample conditional association summary statistics). Repeat both alcohol and AIS conditional analyses before interpreting multiple signals.
3. If MR is considered, use appropriate independent exposure instruments, matched/population-consistent allele reference, assess overlap and horizontal pleiotropy, and check whether models can distinguish alcohol-mediated and direct metabolic effects.
4. Review existing GIGASTROKE downstream consumers of \`variant_pair_id\`: EAS historical IDs may encode other:effect, EUR allele-pair IDs appear lexicographically sorted; in both cases **REF:ALT requires FASTA anchoring**.
