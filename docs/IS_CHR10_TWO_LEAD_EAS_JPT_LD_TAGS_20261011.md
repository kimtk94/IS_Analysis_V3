# IS chr10: EAS/JPT two-lead LD-tag architecture

2026-10-11. BBJ IS locus lead: `10:105459834:T:C`, GIGASTROKE EAS AIS locus lead: `10:105599770:T:G` (source's opposite REF/ALT ID normalized). These lead positions have EAS504 r²≈0.00737; JPT104 r²≈0.02375. This is **public-reference genotype LD only**, not in-study conditional association or SuSiE credible sets.

The script `scripts/is/audit_is_chr10_two_lead_ld_tags.py` calculates ALT-oriented signed LD with each lead for every typed 1000G Phase 3 EAS504 SNP in the full-cis chr10 export. Monomorphic SNPs are excluded from LD tagging.

| Reference panel | Lead | Informative SNPs | Tagged SNPs with r²≥0.8 | r²≥0.2 | CALHM2 EBI QTD000161 nominal p<1e-4 and r²≥0.8 |
|---|---|---:|---:|---:|---:|
| EAS504 | BBJ | 2,778 | **6** | 122 | **3** |
| EAS504 | GIGASTROKE | 2,778 | **11** | 132 | **0** |
| JPT104 | BBJ | 2,666 | **8** | 101 | **3** |
| JPT104 | GIGASTROKE | 2,666 | **11** | 136 | **0** |

The QTL evidence alignment is suggestive of molecular variation near the BBJ lead, but it **cannot** validate colocalization without properly aligned molecular donor LD, validated effective expression N, prior sensitivity and multi-causal modelling. An arbitrary r²-tagging threshold is not a posterior credible set; tagged SNP counts are not independent associations. EAS includes JPT and is not an independent replication sample.

Results: `/srv/is-analysis/results/is/audits/is_chr10_two_lead_ld_tags_20261011_v1/`. Files: `IS_CHR10_EAS_JPT_TWO_LEAD_PER_SNP_LD.tsv`, `IS_CHR10_EAS_JPT_TWO_LEAD_LD_QC.tsv`, and `IS_CHR10_EAS_JPT_TWO_LEAD_LD_MANIFEST.json`. The prior 2,225-gene IS universe, 646 historical coloc outputs, and 120 ABF sensitivity scenarios remain intact.
