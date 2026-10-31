# GSE313736 — original investigators' **2024** public provenance and 59-vs-57 interpretation

**Audit date:** 2026-11-01 (KST).

## Executed source check

[GitHub Actions run #41220669040](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41220669040) **SUCCESS**; [artifact #13049344283](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41220669040/artifacts/13049344283) was downloaded, opened, and JSON verified; all PDF evidence includes public source SHA-256 and exact **printed/PDF page numbers**.

| Verified original public material | Provenance |
|---|---|
| Australian Vision Research (AVR) Annual Report **2024**, p. **20** | [Official 2024 report PDF](https://australianvisionresearch.org/wp-content/uploads/2025/08/2024_AVR-Annual-Report.pdf); **17,319,749 bytes**, SHA256 `4eac68f80d792f045ae0d3ebb5aa41b66796c9be670d92134c9c6324ea0de006`. Chief investigator **Alex W. Hewitt**; co-investigators **Joseph Powell, Seyhan Yazar**. Original single-cell eQTL/glucocorticoid work on 74 TM donors (45 men, 29 women), 316,799 singlets, **63 increased and 17 decreased** DEX-associated genes. Also reports 164,378 cis-eSNPs and 6,178 eGenes; reported as presented at the Australian and New Zealand Glaucoma Society meeting and **under review** at the time of this 2024 grant report. |
| University of Tasmania College of Health and Medicine Research Showcase **2024**, PDF p. **34** | [Public 2024 presentation abstract booklet](https://mooc-media.utas.edu.au/files/7/5/171877818175/CoHM%20Research%20Showcase%202024%20-%20Presentation%20Abstracts%20Booklet.pdf). Presenter **Thomas Guinan**, co-authors Helena Liang, Sandy Hung, Rajendra KC, Alexander Barnett, Primrose Mandalawatta, Kirsten Fairfax, David A. Mackey, Stuart MacGregor, Jamie E. Craig and Alex W. Hewitt. PDF SHA256 `923dc9f87c5038c4dfb1ab2ceaa8a94eec0714f54c70d9c4a00de6d8290b8a54`. Abstract reports **59 samples passed single-cell and genotype QC**, in nine pooled single-cell datasets, before genotype-based demultiplexing. |
| GEO accession **GSE313736**, public 2026-09-03 | [Official GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736), reports the same 74 donors/316,799 singlets and 63-up/17-down high-effect claim, but not the 80 individual source-author DEG gene identities in the official GEO Series supplemental file list. |

**IMPORTANT:** source SHA256 strings are frozen above and in the CI artifact. Only a small portion of publicly available copyright-protected PDF content is excerpted in the manifest; original PDFs themselves were not redistributed by GitHub Actions.

## Why the original 59 is not evidence our 57 are wrong

The 2024 UTAS abstract states **59 samples passed both single-cell and genotype QC before genotype-based demultiplexing**, without providing a final named donor list or the post-DEX BeamA classifier subset. Our distinct public-source reanalysis:
- uses paired Vehicle/DEX per donor;
- applies CellTypist **BeamA confidence ≥0.75**;
- demands **at least 30 matrix-joined cells per donor per condition**;
- retains **57 biological paired donors**; and
- performs a separate, more stringent matched **55-donor BeamA-vs-scPred** sensitivity.

**59 ≠ 57 is not automatically a contradiction.** Differences in filtering stages, cell types, demultiplexing or treatment pairing can explain the change, but without original authors' QC table the **identities and causes of any differences cannot be established**. Do not claim that exactly two particular donors were excluded from the author's original 59 or from specific source pools.

## What 2024 publication information establishes

The same source research, including the reported **63 increased and 17 decreased** response genes, was publicly described in the 2024 Australian Vision Research funding report **before the September 2026 GEO accession appeared**. The funding report identifies investigator and presentation provenance, not a peer-reviewed publication nor a gene-ID table. It reports a talk title roughly equivalent to **“Dissecting the Genetic Determinants of Intraocular Pressure”** and a 2024 presentation and source manuscript **under review at that time**. This does **not** verify the existence or status of a later peer-reviewed manuscript by 2026-11-01; any such claim requires a publication DOI and methods.

### Exact fold-change threshold ambiguity

The AVR report printed the high-effect genes as 63 increased with `log2 fold change >3` and **17 decreased also with `log2 fold change >3`**. The 2026 GEO abstract prints **downregulated `log2 fold change <3`**. Neither describes the same unambiguous **signed negative** threshold `log2FC <= -3`; likely editorial errors are possible but **not proven**. In particular, the exact original author's 80 gene set cannot be recreated by assuming the intended sign and applying our edgeR threshold to our separate test universe.

**Existing internal 83/88 hits are OUR gene sets**, not an 80-member match to the authors. 81/83 overlap between our own BeamA thresholds cannot resolve this separate source audit.

## Recommendation and manuscript-ready scope

- **GO:** label the author group's 2024 original presentation provenance and 2026 public GEO deposition separately; cite precise source method and count.
- **GO:** retain current 57/55 within-source computational results with explicit filters, not presume identical donor eligibility to the 2024 59.
- **HOLD:** exact author 63↑/17↓ Ensembl IDs and original significance/FDR model; original 59 single-cell/genotype QC sample identities; any published original paper DOI/title after 2024.
- **NO invented set overlap:** authors' identity-level overlap remains **UNDETERMINED**.
- **Specific evidence request:** source author group (as listed on official GEO and 2024 presentation) can supply the 80-row DEG table with Ensembl stable IDs, gene symbols, original effect scale/log2FC and BH-FDR, underlying contrast/annotation, QC-eligible donor IDs *in an appropriately anonymized form*, original code/version and publication/preprint URL if any. Also request the exact definition of “downregulated (log2FC <3)” and whether the actual cut was <−3.

**Next action:** if no published gene list can be retrieved, use a short official author request and keep thesis overlap as HOLD rather than forcing a false identity match.


## Additional exact ANZGS presentation provenance — 2024-02-16

**Independent primary-source update (2026-11-01):** The [official ANZGS 2024 congress program, PDF page 1](https://anzgsconference.com/wp-content/uploads/2024/02/ANZGS-2024-Program_FINAL.pdf) lists **Mr Thomas Guinan (University of Tasmania)** for **Friday, 16 February 2024, 17:00–17:15 AEDT**, talk **"Expression Quantitative Trait Mapping of Trabecular Meshwork Cells at Single Cell Resolution."** This independently verifies the talk was **scheduled** before the **3 July 2024** [UTAS Research Showcase abstract](https://mooc-media.utas.edu.au/files/7/5/171877818175/CoHM%20Research%20Showcase%202024%20-%20Presentation%20Abstracts%20Booklet.pdf) that described 59 QC-passing samples.

**Evidence scope:** the official program establishes speaker, date, talk title and slot — **not** the exact presentation's delivered content or that its 63 up/17 down results existed already by February. Those numeric DEG counts are supported separately by the **2024 reporting-period** AVR annual report. Its PDF is hosted at a URL dated **2025/08**; do not substitute the report's coverage year for a verified public publication date. The program, UTAS abstract and AVR report are three distinct archival resources, not three independent datasets or three replications.

**Automated re-audit:** [ocular-gse313736-2024-original-audit.yml](../../.github/workflows/ocular-gse313736-2024-original-audit.yml) now checks all three sources, their exact documentary markers, SHA-256 and compact excerpts, and writes a three-source manifest. [CI run #41224782224](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41224782224) was triggered on the ocular feature branch; its **final conclusion/artifact must be checked separately** before labeling three-source CI PASS. The earlier two-source [run #41220669040](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41220669040) remains valid historical evidence.

**Unchanged HOLD:** original 80 Ensembl gene IDs, exact effect/threshold model, author QC-eligible donor identities and specific reasons behind 59 vs 57. Do not infer original author gene-list overlap using our 83 large-effect genes.
