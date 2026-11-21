# SBP & T2D permission-gated input handoff for Orca

This is a **user-controlled** follow-up to the P1 source selection, and does not try to bypass any distribution restrictions.

## Only user action needed before validation

1. **Codex/Orca**: In the currently open Codex tab, select **Sign in with ChatGPT** and complete the login in the browser. Never send a password or auth token in chat.
2. **SBP**: Prefer authorized download of the actual Evangelou 2018 original summary data file `UKB-ICBPmeta750k_SBPsummaryResults.txt.gz`. The historical NIH GRASP link currently redirects. Alternatively, access OpenGWAS `ieu-b-38` under the provider's account/JWT rules, download the VCF.gz to a local secure directory, and preserve original metadata. <https://opengwas.io/datasets/ieu-b-38>
3. **T2D**: Review and accept the DIAGRAM/T2DGGI download terms on the official consortium page for **2024 EUR ancestry**. Obtain the exact ZIP and preserve it in a secure directory outside Git. <https://www.diagram-consortium.org/downloads.html>
4. Do not send full GWAS data or any credentials to chat. The local input verifier accepts paths directly, reads file metadata only and generates a compact JSON evidence file.

## Prepared CLI: authorization-safe QC

```powershell
# Run inside the Orca IS_Analysis_V3 worktree, after downloading data into a private folder.
py -3 scripts/ckm_stroke/audit_authorized_sbp_t2d_file.py --trait SBP --study ICBP_Evangelou2018_original_GRASP --file "D:\\GWAS\\SBP_original.gz" --out artifacts/ckm_stroke/sbp_supplied_qc.json
py -3 scripts/ckm_stroke/audit_authorized_sbp_t2d_file.py --trait T2D --study T2DGGI_Suzuki2024_EUR --file "D:\\GWAS\\T2DGGI_2024_EUR.zip" --out artifacts/ckm_stroke/t2d_supplied_qc.json
```

These are **example paths**, not verified file locations. The tool will not download or extract source files or store absolute paths. It checks archive layout, lists ZIP contents safely, computes SHA256/MD5 over source bytes, and tests T2D 2024 EUR ZIP against published MD5 `46f13bf073c5560d9b068b623229ff01` when that exact release is supplied.

**PASS_ARCHIVE_MD5_ONLY** is not LDSC-ready. SBP has no independently pinned publisher checksum: its result stays **HOLD_NO_PUBLISHER_CHECKSUM** even if gzip headers and SHA256 pass.

## Separate scientific go/no-go controls

- Confirm population ancestry, variant-level N, allele direction, build (GRCh37 if possible), effect scale and file version. Do not confuse SBP in mmHg, BMI beta, T2D log-odds.
- Rule out duplicated rsIDs/multiallelic mappings; prepare genome build and strand harmonization.
- Verify 1000 Genomes EUR LD reference/HapMap3 marker coverage and appropriate exclusions (MHC, ambiguous palindromic SNPs).
- Investigate UK Biobank overlap in GIANT, ICBP, GIGASTROKE and T2D sources, and quantify implications for genetic correlation uncertainty.
- Only after complete readiness release can `tasks.json` P2 LDSC commands be enabled. Current configuration remains deliberately empty and scientific HOLD.
