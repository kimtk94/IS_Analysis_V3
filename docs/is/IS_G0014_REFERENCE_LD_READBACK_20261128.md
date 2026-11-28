# G0014 rs2767060 reference LD readback — 2026-11-28

Status: REFERENCE_LD_CONTEXT_ONLY / SCIENTIFIC_VALIDATION_FALSE.
Run 46138964839, job 140960450726, code b48090f2717eeb515f64bc759c3aff34d2a9e9b7; completed SUCCESS. Actual aggregate job stdout; artifact 15261171566 metadata SHA256 5180b4c2e629b02fa858b822f00fa9acaac50241856d11298a733da8888bb37c, ZIP not independently downloaded. Source mirror S3 1000 Genomes phase3, primary FTP 404. GRCh37 35,622,291–36,622,291. 29,844 VCF rows, one anchor, 1,103 excluded non-PASS/non-biallelic SNVs, 28,740 partners attempted per population.

| Population | T / haplotypes | finite r2 | monomorphic | r2 >= .8 |
|---|---:|---:|---:|---:|
| CHB | 16/206 | 4281 | 24459 | 22 |
| JPT | 16/208 | 3599 | 25141 | 31 |
| CEU | 32/198 | 4441 | 24299 | 37 |

55 unique high-LD proxies in union, 13 common to all three; rs2845915 C/G r2 approximately 1 in all three, rs2845914 G/A r2 0.93259 CHB and 1 JPT/CEU. Palindromic rs2845915 requires careful allele handling in summary-statistic joins. Distinct phase haplotype LD and integer contingency phi-squared match to maximum absolute error 3.33e-16; this is algorithm parity, not independent cohort replication. Small reference allele counts, reference population composition, not GWAS-matched LD. T frequency higher in CEU; do not infer cause of LAS effect contrast or an EAS-specific variant. Shared high-LD proxy list must be retrieved before novelty classification.

First-source fingerprints (not independent source validation):
- VCF header range 67992f7d02a936121f9928067d115f005cf7e7fff369203596f687074197d8e6
- index 0c3d42abd86cb059c778add85a2d1affa9e6aff54c0118d5beebdb0858f91ed8
- region stream fe8eebbdf045a6122d81a5ae5418021ae23f6f71206e164c9cd7335ca77129eb
- sample panel b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858
- Ensembl mapping 41e90272d974c92f073cfa57a1e70444cb75d32afa9ef33d94832c876f9c9389

Anchor eQTL follow-up run 46139159621, commit 513e864f7cd9648b5d4a3dbd4b262da3dd63b41e, queued at first readback. Fixed 27-cell family includes missing cells, conservative nominal p x 27, source hashes and denominators checked against run 46138809556; diagnostic only, no tissue-level eGene significance, no formal coloc. No main merge or schedule changes.
