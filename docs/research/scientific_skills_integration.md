# Scientific Agent Skills CI pilot

This integration is **metadata-only**, not a scientific validation or automatic agent execution.

- Upstream: https://github.com/K-Dense-AI/scientific-agent-skills
- Pinned commit: `5f2bc0341354fd5590aef76a533ada1e1e5b4e62`
- Execution: GitHub Actions only; no SSH, local, or Colab jobs.
- Four initial mappings: IS → database-lookup; glaucoma → ontology-term-resolution; HCC → depmap; radiomics → imaging-data-commons.
- The audit reads each SKILL.md and records its SHA-256, provenance and existence; it **never executes third-party instructions or scripts**.
- No source datasets, patient identifiers, secrets or protected data are sent to upstream.
- PASS_METADATA_ONLY means the pinned skill document exists, not that GWAS, RNA-seq, HCC aliquot mapping, or radiomics are scientifically validated.
- Expansion to the other four studies and automated use of skill guidance require independent review, sandboxing, data provenance checks and existing study-specific scientific gates.
