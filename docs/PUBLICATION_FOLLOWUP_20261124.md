# Publication follow-up — 2026-11-24

This is a follow-up to the first focused public-analysis run. No Radiomics, harness, private data uploads or SSH.

- Ocular: the original public P5d ZIP is pinned by SHA256, and the official TM workbook is retrieved again via the same parser. Export the 241 strong concordant and 64 sign-discordant genes only if the previous counts match. Exploratory Reactome overrepresentation uses the unique exactly matched gene symbols with external adjusted P available and Reactome annotation. Restrict term sizes to 10–500 background genes, perform one upper-tail hypergeometric family, include zero-overlap terms in BH. Preserve the current Reactome archive hash. Pathway significance is descriptive/conditional on selection, not a new independent experiment. Source: https://reactome.org/download/current/ReactomePathways.gmt.zip.
- Skin: use the two provider files actually listed by the first run, their exact named MD5 manifests and full row-level scan when obtainable. A bounded attempt can end in HOLD. Do not interpret row counts as independent loci or claim MR from a successful scan.
- Oral: use the five original GEO metadata labels documented by the Bioconductor GSE103322 vignette (https://bioconductor.org/packages/release/data/experiment/vignettes/GSE103322/inst/doc/GSE103322.html). Site/cancer/enzyme cell cross-tabs are descriptive only. No patient identifiers or expression are exported, and patient-level hypothesis testing is deferred until authoritative patient linking is verified.

All computed results remain Actions artifacts, not Git source files. The previous result gates and the other six research tracks remain unchanged.
