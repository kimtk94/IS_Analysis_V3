# Ocular seven-day pilot: technical readback and full-cohort resource decision — 2026-11-28

Status: TECHNICAL_PILOT_COMPLETE / HOLD_FULL_TRANSFER_APPROVAL / SCIENTIFIC_VALIDATION_FALSE.

Source: run 46130367058, commit ac15290f4a93185d1917cbe820034915059221c0. Actual job stdout and inventory readback; not an independent artifact ZIP extraction. Library metadata pairing source: docs/publication/evidence/ocular_seven_day_library_pairing_20261128.json (blob 20b51fb9cf482d98003160c7875b07fa0666b1db). Pair 18-39 is the smallest compressed paired FASTQ pair, not a biologically selected pair.

| Run | Arm | Processed fragments | Pseudoaligned | Percent | Gene rows | Elapsed seconds | Job ID |
|---|---|---:|---:|---:|---:|---:|---:|
| SRR14535091 | DEX | 29595850 | 25271172 | 85.4 | 62700 | 3940.23 | 140933987737 |
| SRR14535092 | ETH | 21004317 | 18270452 | 87.0 | 62700 | 3135.68 | 140933987736 |

Both QC reports are TECHNICAL_COMPLETE_QC_REVIEW_REQUIRED. Inventory job 140949941014 verifies gene table hashes; artifact 15259034100, ZIP SHA256 56cba154c851368ad9626d4188d53e41cb29bbb56a9ed8a4c7c5047e5a3f9f72. Gene table SHA256 DEX 13ab69af8ea76dac33e45d03664333502e2e3c03a6844c2e58fb0796c93ba8c5; ETH e9339bb061a55921c9e61ae4f7ef5f42609daea2ba5349eb3e9ce699be5d5585. Kallisto 0.48.0; GENCODE v44; singleton FASTQs excluded uniformly. Counts and pseudoalignment alone do not establish biological quality or treatment response.

## Resource inventory

Recalculated from the pinned ENA metadata blob 9aef45b75089f69304c3e850fe382df70a12009f and corroborated pilot total in job 140933940791. The two large paired mates per run are the compressed R1/R2 files (singleton excluded; original workflow checks exact file roles). Full 16-run sum 74883257207 bytes; pilot 8673362787 bytes; remaining 66209894420 bytes if pilot artifacts are reused under a future explicitly reviewed provenance-compatible plan. Current full workflow **reruns all 16** and requires approve_full_75gb=true; no bypass, no dispatch here. Pilot total processing time 7075.91 seconds across two jobs, not a measured full-cohort cost. Budget depends on actual runner billing and retry traffic; 75 GB is nominal paired FASTQ only, not total network transfer including reference/retries.

## Decision gate

- Pair metadata: 8/8 pairs identified, 4 responder and 4 nonresponder; pilot only 1/8 technical pairs, 2/16 libraries. Do not select or drop pairs using expression outcomes.
- 16h technical run SRR16961557 gene TPM aggregation confirmed separately in run 46130181827, job 140933382384; not an independent donor validation and not a seven-day run.
- Full-cohort sensitivity gate run 46130518397: contract job 140934407065 passed 9 synthetic checks, actual full-sensitivity skipped. No paired score or biological outcome read.
- Explicit full-transfer approval and branch protocol review remain required. Verify artifacts are still available, all 16 identities, hashes, reference, mapping, quantifier, library QC and complete pair coverage; preserve all exclusions. Current sensitivity workflow requires a successful full dispatch and a separate reviewed_qc_and_protocol approval. Do not quietly combine partial artifacts or run biological inference on one pair.
- No change to panel, thresholds, analysis code, schedules, or main. No new FASTQ transfer or CI claim. Technical feasibility improved; novelty/independence and scientific validation unchanged.
