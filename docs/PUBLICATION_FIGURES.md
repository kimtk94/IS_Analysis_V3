# Free scientific figure system (v1)

Optional Python/Matplotlib pipeline for aggregate-only, draft research graphics.
No paid SaaS, generative-image service, new scientific computation, server mutation or publication claim.
Use Inkscape (https://inkscape.org/) for editable SVG finishing, Cytoscape (https://cytoscape.org/) for complex biological networks, and original researcher-reviewed diagrams or separately licensed Servier Medical Art (https://smart.servier.com/) for medical illustrations. Verify each asset's current license and attribution. No third-party asset is bundled in v1.

## Quick start
    python -m pip install -r config/publication-figures-requirements.txt
    python scripts/publication_figures.py --manifest tests/fixtures/publication_figures/demo.json --input-root tests/fixtures/publication_figures --out /tmp/figure-demo

Demo numbers are synthetic, visibly watermarked. Never use them as research results.
An existing output folder is rejected. Outputs: editable SVG, vector PDF, 600-dpi PNG/TIFF, provenance JSON with input/code/output SHA256, references, package versions, row counts and raster dimensions.

## Inputs
Copy the example manifest outside Git; set data_status to aggregate ONLY for verified aggregate outputs; replace the fixtures and source references. Source paths must be relative to input root. Never upload real results, private/controlled records, donor identifiers or source archives to public GitHub or CI. The aggregate flag is a declaration, not automatic de-identification. Render real outputs locally or in a separately approved private execution environment.

| Panel | TSV columns / manifest fields | Scientific scope |
|---|---|---|
| Forest | label, effect, lower, upper; scale ratio/difference, effect_label | Existing estimates and intervals only; verify CI level and effect scale |
| Volcano | label, log2fc, p; p_label, p_threshold, fc_threshold | Existing effects and P/FDR, display thresholds only |
| Heatmap | label followed by numeric columns; units | No implicit normalization or clustering; zero-centered scale |
| Locuscompare | variant, p_gwas, p_qtl | Unique harmonized variants, descriptive only, no colocalization posterior |
| Schematic | nodes with id, label, x, y; edges with from/to/kind | workflow / documented / hypothesis; documented edges require references, hypotheses dashed |

Panel references mandatory. A documented edge reference is traceability, not proof that the citation supports the claim. Custom anatomical illustration and visual QA remain manual. Keep captions with cohort, ancestry, units, method, denominator, uncertainty, multiplicity, accession and exact analysis revision. Zero P values are rejected, never fabricated.

## CI and release gates
Path-scoped CI renders only synthetic fixtures, runs focused tests and existing repository smoke tests, uploads a seven-day demo artifact. No real-data adapters yet. Before submission: scientific and source review, visual quality, journal-specific font/width/DPI/file requirements, asset licenses and attribution. Keep editable originals. Researchers remove the draft watermark only after these gates.

## Next steps
1. Verify CI and visually inspect demo exports.
2. Inspect actual aggregate contracts on active IS/CKD and ocular branches; add read-only adapters without changing scientific methods or evidence gates.
3. Curate anatomical SVG assets with individual license and attribution manifests.
4. Add MR sensitivity, Manhattan, regional association and cell-state plots, then journal profiles.
