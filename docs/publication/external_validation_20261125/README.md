# External validation request packet — draft, not sent

24 previously frozen, nonpalindromic markers across 16 chromosomes. This is a request for authorized aggregate data only. Variant identity derived from verified EAS AIS anchor run 45534891070 job 138938960219 artifact 14986197928 (upload-reported SHA256 5c978c3fc19cb61c2a9f2fb7e8de99ba81ccf3a0d30b0f9d1371c344e9dc7aef); no new independent FASTA check, do not treat historical REF as newly certified reference. The marker file contains the exact requested analysis allele, not a risk allele assertion. Multiallelic rs7203740 must match C/T, not C/A. No proxy replacements.

## Files
- `marker_request.tsv`: all 24 marker identities, no discovery effects or P values, for the collaborator.
- `aggregate_response_template.tsv`: all 24 x 4 phenotype rows (96), blank returned statistics; NOT_RETURNED by default. A blank is not a zero or negative result. Return explicit reasons for unavailable markers, not just significant results.
- `cohort_metadata_dictionary.tsv`: questions for eligibility and interpretation.
- `internal_eas_ais_reference.tsv`: discovery reference for internal harmonization; not necessary to send to avoid directional selection.
- `OUTREACH_DRAFT.md`: unsent message.
- `PROTOCOL_DRAFT.md`: proposed endpoint and QC rules, to finalize before access/readout.
- `anchor_readback.json`: verified source aggregate (public GWAS estimates only).

All files in the packet are drafts for review. No external data obtained or sent, no person-level data requested, no automatic upload of external data. Keep authorized responses in approved storage; never commit restricted extracts, identifiers or credentials to GitHub. Independent replication remains HOLD.
