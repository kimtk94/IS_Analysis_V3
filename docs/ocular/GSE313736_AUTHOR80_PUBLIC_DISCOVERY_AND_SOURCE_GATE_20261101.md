# GSE313736 original-author DEG80 — public evidence gate (2026-11-01 KST)

## CI execution and boundary
- **Executed SUCCESS**: [GitHub Actions #41225415076](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41225415076), SHA `6face31a5296de2172d874c45a170969ab2cd0c1`; two jobs both passed.
- [Actual source-check artifact #13051518029](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41225415076/artifacts/13051518029), unexpired when audited.
- Source script: `scripts/ocular/ocular_gse313736_author80_public_gate.py` (Python standard library only).
- Workflow: `.github/workflows/ocular-gse313736-author80-public-gate.yml`. Manual and ocular-branch path-scoped push runs; no recurring 10x matrix redownload.
- **PASS means the public endpoint requests and the checker operated**, not that the original 80 genes were retrieved, or that a discovered paper is author-authenticated.

## Actual live CI log results
| Public check | Log | Evidence interpretation |
|---|---|---|
| Europe PMC `EXT_ID:GSE313736` | `CHECKED 0` | No indexed hits under this exact identifier query |
| Europe PMC `"GSE313736"` | `CHECKED 0` | No indexed hits under the literal string query |
| Europe PMC `"Thomas Guinan" AND "trabecular meshwork"` | `CHECKED 0` | No indexed hits under this restrictive author/subject query |
| Official GEO Series FTP supplement directory | `CHECKED []` | No candidate tabular filename from the directory listing |
| Structural validator | `AUTHOR80_GATE_SELF_TEST_PASS` | Synthetic IDs, duplicates, source-provenance requirement, 63/17 count, signed overlap worked |

These are **bounded negative searches**, not proof a work or a supplemental list does not exist anywhere. In particular the official GEO **GSM** files include raw count matrices and associated metadata, which are **not** an author DEG80 table. The GEO Series inventory checked here should never be conflated with individual GSM inventories. The 2024 AVR grant report remains the source of the *reported count*, not gene identities.

## Strict source-author comparison acceptance
A downstream author-vs-ours comparison requires BOTH actual tables:
1. Source-author DEG table with 80 unique Ensembl stable gene IDs, numerical log2FC, reported FDR, and **63 positive/17 negative** log2FC; file hash is frozen. It must be traceable to an author-owned/publication artifact rather than guessed source.
2. Our existing own gene-level table from [P5 corrected primary edgeR #41074586446](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446) with stable gene IDs and effects; our 83 high-effect genes **must not** be relabeled as the original author 80.
3. Original contrast/cell annotation, test universe, statistical method, significance/effect threshold and **actual signed downregulated cutoff** must be authenticated before a scientific reproduction claim.

Checker CLI, once an **authentic** author file has been obtained:

```bash
python scripts/ocular/ocular_gse313736_author80_public_gate.py \
  --out results/ocular/author80-gate \
  --author-table /path/to/authenticated_author80.tsv \
  --ours-table /path/to/our_p5_gene_level.tsv \
  --author-source-url 'https://official-source.example/real-file' \
  --author-method 'verified author contrast and criteria'
```

**The CLI fields are declarations, not proof of authenticity.** Output is limited to reviewable aggregate set overlap and signed direction, never GO based on a provided URL string. A study that provides gene symbols only, ambiguous thresholds, or no author-owned reproducible source remains **HOLD** until reviewed. No individual-level genotype, patient metadata, or per-donor raw count matrix is needed.

## Research decision
- **GO:** historical 2024 study source provenance (AVR, UTAS, ANZGS), original 57-donor edgeR computational reanalysis, matched 55-donor annotation sensitivity, new public-source discovery workflow operation.
- **HOLD:** original authors' 80 gene identities, independent authorization of an author table, 59 source QC donor identity differences, signed cutoff, true source-author gene overlap.
- **HOLD:** treatment versus separate source 10x library-batch causality; existing 85.46% FDR-significant rate remains a major threat.
- **STOP:** using clinical steroid-induced IOP phenotypes not observed in this cultured donor cohort as if they were measured, or treating same-cohort annotation sensitivity as independent replication.

## Required next action and responsibility
The original investigators' data request is prepared in [ORIGINAL_AUTHOR_EVIDENCE_REQUEST_PACKET_20261101.md](ORIGINAL_AUTHOR_EVIDENCE_REQUEST_PACKET_20261101.md). **Not sent.** Official GEO contact/affiliation must be checked before a user-approved outreach. Ask for the 80-row table with FDR and stable IDs, actual signed cutoff and model, anonymized QC eligibility counts, and any manuscript/preprint. There is no justification for repeating the 18 original raw 10x matrix downloads while author80 identity remains blocked.

## Follow-up 2026-11-02 — strict 80 versus 83 structural comparison

**CI PASS:** [Actions #41324876650](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41324876650), source commit `b88a6f1e49dd532ef0d0a5a473e0b7ecd95eb3f6`, two parallel jobs SUCCESS. Verified job log:
- `AUTHOR80_GATE_SELF_TEST_PASS: 63/17, derived 67/16 high-effect 83, tested-vs-DEG separation, zero effect, provenance, duplicate`.
- Three Europe PMC queries returned zero matching records; official GEO Series supplement directory yielded no candidate original author DEG table. This does not establish global absence.

**Important analytical correction:** `ours-table` can be the **full P5 tested-gene table**, not merely a prefiltered 83-row file. The validator now derives our high-effect set explicitly using `FDR < 0.05 AND abs(log2FC) >= 3`, then requires **83 = 67 up + 16 down** before permitting a reviewable original-author comparison. It reports **both** (a) number of 80 author genes appearing anywhere in our full tested universe, and (b) number appearing in our own 83 high-effect genes, plus direction concordance in the *latter intersection*. This avoids falsely treating detectable but not differentially expressed genes as matched DEG hits. Exact zero `log2FC` can occur in the full tested universe and is now allowed in that input; author 80 DEG set requires nonzero effects.

**Official author contact confirmed:** [NCBI GEO GSM9373290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9373290) lists Alex Hewitt, University of Tasmania / Menzies Institute, institutional e-mail `alex.hewitt@utas.edu.au`. The [author request packet](ORIGINAL_AUTHOR_EVIDENCE_REQUEST_PACKET_20261101.md) now names the official address. **No email has been sent.**

**Decision:** structural CI capability **GO**; actual author original DEG80 **HOLD**; neither same-cohort sensitivity nor a synthetic self-test can be described as external gene-level replication. Prior source-analyzed technical confounding remains **HOLD**.
