# Ocular P0 data access verification (2026-10-30)

## Status vocabulary
- `LISTED`: GEO/Mendeley advertises a downloadable file; never treat as transfer success.
- `TRANSFER_PASS`: CI downloaded real bytes, SHA256 computed, table/gzip header checked.
- `PARTIAL_G0_PASS`: series SOFT metadata plus at least one selected small data file downloaded; this is **not** full analysis readiness.
- `HOLD`: failed transfer, incomplete phenotype, insufficient donor pairing or restricted study.
- `HOLD_CONTROLLED`: dbGaP approval required.
- `GO_MOLECULAR`: G0+G1+G2+G3+G4+novelty pass; explicitly no human clinical response claim.
- `GO_CLINICAL`: additional clinically measured delta IOP with authorized GWAS summary or patient genotype+phenotype.
- `STOP`: incompatible phenotype, lack of inferential unit or no novelty after high-impact paper comparison.

## Priority and evidence gates
1. **P0** GSE313736: linked DEX/vehicle donor metadata, scRNA counts, demultiplexing maps; independent donor as unit, no cell pseudoreplication. SOFT-derived overlap is hypothesis only until demultiplex QC and count matrix matching.
2. **P0** GSE301517: control/DEX TM RNA; public CPM table is normalized, **not raw counts**. Biological donor count TM62/TM70/TM72 must not be conflated with 21 libraries.
3. **P0** GSE301516: control/DEX ATAC peak sets; pooled peak coordinates are insufficient for replicate differential accessibility.
4. **P0** phs004037.v1.p1: clinical steroid-induced delta IOP GWAS/WES, **controlled**. Do not download genotype files or mislabel ordinary glaucoma/IOP GWAS as steroid response.
5. **P0** DME nszmnfcwhs v1: published between-responder aqueous proteomics tables; patient x protein matrix, longitudinal OCT/BCVA and drug metadata are necessary for clinical prediction. If only DOCX group contrasts, DME therapy-response model = HOLD.

## CI run
The new GitHub Action performs bounded transfers of small publicly advertised GEO files, parses sample-level SOFT metadata, extracts inferred donor overlaps, checks compressed table headers and SHA-256, lists Mendeley files, tests download of small DOCX when possible, and writes `ocular-p0-audit.json`. It does not transfer massive single-cell MTX, FASTQ, genotype VCF, or individual clinical records.

Go/no-go decisions must use **actual workflow artifact**, not the existence of this markdown or CI YAML. If no workflow run/artifact is available, all transfer-dependent gates remain unverified.

## Next decision
Steroid molecular study = **HOLD pending byte transfer/donor-linkage QC** (candidate for conditional GO). Steroid clinical genomics = **HOLD controlled phenotype/genotype**. DME anti-VEGF prediction = **HOLD patient-level matrix**. Generic DR progression = distinct fallback topic; cannot substitute for DME treatment response.

Results: GitHub Actions workflow `ocular-p0-access-audit.yml`, artifact `ocular-p0-access-evidence`. Do not commit downloaded data, medical genotypes, credentials or private clinical records.
