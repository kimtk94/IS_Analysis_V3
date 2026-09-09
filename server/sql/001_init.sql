CREATE TABLE IF NOT EXISTS pipeline_runs (
    run_id TEXT PRIMARY KEY,
    started_at TIMESTAMPTZ NOT NULL DEFAULT now(),
    finished_at TIMESTAMPTZ,
    status TEXT NOT NULL,
    git_commit TEXT,
    ancestry TEXT,
    phenotype TEXT,
    message TEXT
);

CREATE TABLE IF NOT EXISTS source_artifacts (
    id BIGSERIAL PRIMARY KEY,
    dataset TEXT NOT NULL,
    gene TEXT,
    ancestry TEXT,
    source_file TEXT NOT NULL,
    source_uri TEXT,
    checksum_sha256 TEXT,
    size_bytes BIGINT,
    raw_rows BIGINT,
    standardized_rows BIGINT,
    status TEXT NOT NULL,
    run_id TEXT REFERENCES pipeline_runs(run_id),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT now(),
    UNIQUE (dataset, ancestry, source_file)
);

CREATE TABLE IF NOT EXISTS analysis_results (
    id BIGSERIAL PRIMARY KEY,
    run_id TEXT REFERENCES pipeline_runs(run_id),
    gene TEXT NOT NULL,
    protein_id TEXT,
    ancestry TEXT NOT NULL,
    phenotype TEXT NOT NULL,
    stage TEXT NOT NULL,
    metric TEXT NOT NULL,
    value_numeric DOUBLE PRECISION,
    value_text TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_results_gene_stage
    ON analysis_results(gene, stage);

CREATE TABLE IF NOT EXISTS pipeline_alerts (
    id BIGSERIAL PRIMARY KEY,
    run_id TEXT REFERENCES pipeline_runs(run_id),
    severity TEXT NOT NULL,
    alert_type TEXT NOT NULL,
    gene TEXT,
    message TEXT NOT NULL,
    resolved BOOLEAN NOT NULL DEFAULT false,
    created_at TIMESTAMPTZ NOT NULL DEFAULT now()
);
