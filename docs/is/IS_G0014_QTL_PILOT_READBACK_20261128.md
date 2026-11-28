# G0014 rs2767060 QTL pilot readback — 2026-11-28

Status: REGIONAL_INVENTORY_ONLY / SCIENTIFIC_VALIDATION_FALSE.

Source run 46138809556, job 140959965946, code 3db6ed324cbac5a32487d2a3a92dba72db58fb78. Actual aggregate job stdout readback; artifact 15261390343 metadata digest sha256:eaad7ffa3f8009c55e9422c0dc5444ac290b2d1eda50afd32c2dabcb9380e34d, ZIP not independently downloaded. First run 46138759270 failed before data access due to reserved Python package name; corrected with importlib. Preserve both runs.

Fixed GRCh38 chr10 35,333,363–36,333,363; nine positional genes from immutable 0bf646fd6f29322d17da878253f37def24c01c74, all retained. Original indexed source regional rows 33,328 Aorta, 36,291 Cortex, 37,796 ROSMAP DLPFC; sum 107,415, not unique variants.

| Gene | GTEx Aorta QTD000131 | GTEx Cortex QTD000171 | ROSMAP DLPFC QTD000434 |
|---|---:|---:|---:|
| CCNY | 3070 | 2826 | 3257 |
| FZD8 | 3407 | 3149 | 3630 |
| GJD4 | 0 | 3149 | 3630 |
| RP11-382K22.1 | 0 | 0 | 0 |
| RP11-309N24.1 | 0 | 0 | 0 |
| RP11-810B23.1 | 3407 | 0 | 0 |
| MIR4683 | 0 | 0 | 0 |
| RP11-324I22.4 | 3407 | 3149 | 3630 |
| RP11-425A6.5 | 0 | 3149 | 3630 |

6/9 genes and 14/27 gene-source cells have rows; all 14 have one anchor C/T position row and all counted rows pass the pilot numeric-field checks. The 13 empty cells are unobserved, not null molecular associations. No anchor beta/p inference, uniqueness, full cis, FASTA harmonization, sample N, sdY, ancestry-matched LD, or disease colocalization established. Formal coloc 0/27.

LD follow-up workflow .github/workflows/is-rs2767060-phase3-regional-ld-20261128.yml, code b48090f2717eeb515f64bc759c3aff34d2a9e9b7, run 46138964839 queued at readback. Initial chr10 source hashes recorded, not independently pinned. Reference LD is not GWAS-matched LD. No change to main, #362, or scheduled research. Next: verify LD actual readback, then independently validate source fingerprints, test regional coverage and harmonize all qualifying assays without selecting by QTL significance. External LAS replication still HOLD.
