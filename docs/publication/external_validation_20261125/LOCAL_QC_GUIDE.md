# Local aggregate response QC — before scientific validation

The 24 x 4 template is immutable in marker/phenotype identity. This parser validates 96 rows, exact coordinates/alleles, orientation, explicit missingness, finite estimates, sample counts, and optional prelocked model-specific quality/count rules. It does not perform scientific hypothesis tests. Use approved local storage only for any external responses; the GitHub workflow uses synthetic tests and the blank template only. Do not copy external responses or normalized output into GitHub/Actions, and do not upload external data through this packet.

Run locally on a permitted machine (replace response and output paths with approved paths):

```sh
python scripts/is/external_validation_qc.py --packet docs/publication/external_validation_20261125 --response /approved/path/aggregate_response.tsv --out /approved/path/qc.json
```

Optional policy JSON is a model-ID mapping under "models"; each entry has phenotype, imputation_metric, min_quality in [0,1], min_cases, min_controls. No defaults are assumed, and no real model eligibility has been agreed. Directly genotyped markers require a quality metric and a documented policy; do not invent INFO. per_marker_n is total analyzed individuals for the declared case-control model, not effective sample size. If the collaborator uses a different definition, leave the record unresolved until metadata/protocol are updated. Use explicit statuses for inaccessible estimates. No authorizations, phenotype equivalence, cohort non-overlap, external direction, significance or heterogeneity are certified by automated QC. Only the separate prelocked protocol and approved investigator review can support those determinations.
