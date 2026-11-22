#!/usr/bin/env python3
"""2026-11-02: eight independent bounded public-source probes; scientific gates fail closed.

Run ONLY inside GitHub Actions. Never publish patient/sample identifiers or raw matrices.
A source/metadata PASS does not imply a biological, causal, or clinical PASS.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

TRACKS = {
    "is": ("GWAS Catalog v2 GIGASTROKE outcome accession and ancestry provenance",
           "Full EUR/EAS GWAS summary-statistics integrity, harmonization and LD"),
    "glaucoma": ("GSE313736 official GEO source and raw-count availability",
                 "Nine leave-pool-out quasi-likelihood refits on verified integer counts"),
    "hcc": ("PDC000198/199 official case-sample-aliquot uniqueness and cross-study join",
            "Official quantitative column-to-aliquot and RNA-to-patient mapping"),
    "radiomics": ("LiverHccSeg Zenodo source file/version and imaging-access feasibility",
                 "Read actual 14 patient NIfTI headers and verify PV voxel spacing"),
    "ocular": ("PRJNA1243143 official ENA 46-run exposure-time metadata completeness",
               "Unique 1h/6h/48h run-to-strain mapping and donor-aware raw-count analysis"),
    "muscle": ("Locate exact 20 RT loci and query GTEx v10 variant resolution",
               "20-locus source-verified GTEx coverage and EAS LD proxy analysis"),
    "oral": ("GSE103322 official GEO metadata donor-linkage feasibility",
             "Download expression and patient metadata; verify 18 donor identities"),
    "skin": ("Bristol perceived-age GWAS readme and bounded summary-file header QC",
             "Full summary-statistics integrity, genome build, allele and effect-scale QC"),
}
MAX_SOURCE = 5_000_000
USER_AGENT = "MASTER-DEGREE-8track-CI-source-audit/2026-11-02 (read-only)"


class Hold(Exception):
    pass


class IntegrityFail(Exception):
    pass


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(url: str, limit: int = MAX_SOURCE, headers: dict | None = None,
          method: str = "GET", data: bytes | None = None) -> tuple[bytes, dict]:
    request = urllib.request.Request(
        url, data=data, method=method,
        headers={"User-Agent": USER_AGENT, "Accept": "*/*", **(headers or {})})
    with urllib.request.urlopen(request, timeout=28) as response:
        body = response.read(limit + 1) if method != "HEAD" else b""
        if len(body) > limit:
            raise Hold("RESPONSE_OVER_BOUNDED_LIMIT")
        return body, {"url": url, "http_status": response.status,
                      "bytes_read": len(body), "sha256": digest(body) if body else None,
                      "content_length_header": response.headers.get("Content-Length"),
                      "last_modified": response.headers.get("Last-Modified")}


def json_get(url: str, limit: int = MAX_SOURCE) -> tuple[dict, dict]:
    body, meta = fetch(url, limit=limit, headers={"Accept": "application/json"})
    try:
        obj = json.loads(body)
    except (ValueError, UnicodeDecodeError) as exc:
        raise Hold("SOURCE_NOT_JSON") from exc
    if not isinstance(obj, dict):
        raise Hold("UNEXPECTED_JSON_ROOT")
    return obj, meta


def official_gwas():
    url = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104540"
    obj, meta = json_get(url)
    if "GCST90104540" not in json.dumps(obj):
        raise IntegrityFail("GWAS_STUDY_ACCESSION_MISMATCH")
    ancestry_url = url + "/ancestries"
    try:
        ancestries, ancestry_meta = json_get(ancestry_url, limit=2_000_000)
        ancestry_status = "RESPONSE_RECEIVED_NOT_POPULATION_HARMONIZED"
        ancestry_digest = ancestry_meta["sha256"]
    except (Hold, urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
        ancestry_status = "HOLD_ANCESTRY_ENDPOINT_" + type(exc).__name__
        ancestry_digest = None
    return "PASS_STUDY_ID_ONLY", {
        "study": "GCST90104540", "metadata_sha256": meta["sha256"],
        "ancestry_status": ancestry_status, "ancestry_sha256": ancestry_digest,
        "api_version": "v2", "source": url,
    }, "No source SNP effects, allele harmonization or ancestry-matched LD validated"


def geo_source(accession: str):
    url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
           + urllib.parse.urlencode({"acc": accession, "targ": "self",
                                     "form": "text", "view": "quick"}))
    body, meta = fetch(url, limit=4_000_000)
    text = body.decode("utf-8", errors="replace")
    if not (re.search(r"\^SERIES\s*=\s*" + accession + r"\b", text, re.I)
            or re.search(r"!Series_geo_accession\s*=\s*" + accession + r"\b", text, re.I)):
        raise Hold("GEO_ACCESSION_NOT_CONFIRMED_IN_SOFT")
    return text, meta


def glaucoma():
    text, meta = geo_source("GSE313736")
    return "PASS_GEO_ACCESSION_ONLY", {
        "accession": "GSE313736", "soft_sha256": meta["sha256"],
        "series_sample_links": len(re.findall(r"^!Series_sample_id", text, re.M)),
        "supplementary_link_count": len(re.findall(r"^!Series_supplementary_file", text, re.M)),
    }, "No pool-resolved integer counts or nine leave-pool-out fits validated"


def validate_pdc(rows: list) -> tuple[set, dict]:
    if not isinstance(rows, list) or not rows:
        raise Hold("PDC_BIOSPECIMEN_EMPTY")
    pairs = {}
    sample_cases = {}
    for r in rows:
        if not isinstance(r, dict):
            raise Hold("PDC_RECORD_NOT_OBJECT")
        a, s, c = (r.get("aliquot_id"), r.get("sample_id"), r.get("case_id"))
        if not all(isinstance(x, str) and x.strip() for x in (a, s, c)):
            raise IntegrityFail("PDC_MISSING_CASE_SAMPLE_ALIQUOT_ID")
        if a in pairs and pairs[a] != (s, c):
            raise IntegrityFail("PDC_ALIQUOT_TO_MULTIPLE_PARENTS")
        if s in sample_cases and sample_cases[s] != c:
            raise IntegrityFail("PDC_SAMPLE_TO_MULTIPLE_CASES")
        pairs[a] = (s, c)
        sample_cases[s] = c
    if len(pairs) != len(rows):
        raise IntegrityFail("PDC_DUPLICATE_ALIQUOT_ROWS")
    canonical = sorted((a, s, c) for a, (s, c) in pairs.items())
    return set(canonical), {
        "aliquots": len(pairs), "samples": len(sample_cases),
        "cases": len(set(sample_cases.values())),
        "mapping_sha256": digest(json.dumps(canonical, separators=(",", ":")).encode()),
    }


def hcc():
    all_maps, stats = [], []
    for study in ("PDC000198", "PDC000199"):
        query = ('{ biospecimenPerStudy(pdc_study_id:"' + study
                 + '",acceptDUA:true){aliquot_id sample_id case_id sample_type} }')
        payload = json.dumps({"query": query}).encode()
        body, meta = fetch("https://pdc.cancer.gov/graphql", limit=5_000_000,
                           method="POST", data=payload,
                           headers={"Content-Type": "application/json",
                                    "Accept": "application/json"})
        obj = json.loads(body)
        if obj.get("errors"):
            raise Hold("PDC_GRAPHQL_ERROR_" + study)
        rows = obj.get("data", {}).get("biospecimenPerStudy")
        mapping, stat = validate_pdc(rows)
        stat["study"] = study
        stat["source_response_sha256"] = meta["sha256"]
        all_maps.append(mapping)
        stats.append(stat)
    if any(x["aliquots"] != 331 or x["cases"] != 171 for x in stats):
        raise IntegrityFail("PDC_CATALOG_DENOMINATOR_DRIFT_REQUIRES_REVIEW")
    if all_maps[0] != all_maps[1]:
        raise IntegrityFail("PDC_CROSS_ASSAY_ALIQUOT_MAPPING_MISMATCH")
    return "PASS_OFFICIAL_BIOSPECIMEN_JOIN_ONLY", {
        "studies": stats, "common_aliquot_mapping_count": len(all_maps[0]),
        "quant_column_mapping_verified": False,
        "rna_patient_join_verified": False,
    }, "Biospecimen identity does not map 330 quantitation columns or RNA samples"


def radiomics():
    obj, meta = json_get("https://zenodo.org/api/records/8179129")
    if str(obj.get("id")) != "8179129":
        raise IntegrityFail("ZENODO_RECORD_ID_MISMATCH")
    files = obj.get("files", [])
    if not isinstance(files, list) or not files:
        raise Hold("ZENODO_NO_FILES")
    return "PASS_ZENODO_RECORD_ONLY", {
        "record_id": 8179129, "record_sha256": meta["sha256"],
        "listed_files": len(files),
        "listed_bytes": sum(int(f.get("size") or 0) for f in files),
        "nifti_file_names_exposed": sum(".nii" in str(f.get("key", "")).lower() for f in files),
        "nifti_header_spacing_checked": 0,
    }, "Actual 14-patient PV NIfTI headers are not fetched by a metadata probe"


TIMEPOINT_RE = re.compile(r"(?<!\w)(?:1|6|48)\s*(?:h|hr|hrs|hour|hours)(?!\w)", re.I)


def ena_summary(tsv: str) -> dict:
    reader = csv.DictReader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    if not rows or not reader.fieldnames or "run_accession" not in reader.fieldnames:
        raise Hold("ENA_RUN_SCHEMA_UNAVAILABLE")
    runs = [str(x.get("run_accession", "")).strip() for x in rows]
    if len(runs) != len(set(runs)):
        raise IntegrityFail("ENA_DUPLICATE_RUN_ACCESSION")
    if len(runs) != 46:
        raise IntegrityFail("ENA_46_RUN_COUNT_DRIFT")
    fields = ("sample_alias", "sample_title", "sample_description",
              "experiment_title", "experiment_alias", "study_title")
    mapped = 0
    ambiguous = 0
    for r in rows:
        value = " ".join(str(r.get(f) or "") for f in fields)
        times = {m.group(0).replace(" ", "").lower() for m in TIMEPOINT_RE.finditer(value)}
        if len(times) == 1:
            mapped += 1
        if len(times) > 1:
            ambiguous += 1
    return {"unique_runs": len(runs),
            "unique_biosamples": len(set(r.get("sample_accession") for r in rows)),
            "timepoint_unique_text_matches": mapped,
            "timepoint_ambiguous_text_matches": ambiguous,
            "available_metadata_fields": [f for f in fields if f in reader.fieldnames],
            "fastq_md5_entries": sum(len((r.get("fastq_md5") or "").split(";"))
                                      for r in rows if r.get("fastq_md5"))}


def ocular():
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?"
           "accession=PRJNA1243143&result=read_run&format=tsv")
    body, meta = fetch(url, limit=2_000_000)
    stats = ena_summary(body.decode("utf-8-sig"))
    stats["source_sha256"] = meta["sha256"]
    stats["project"] = "PRJNA1243143"
    return "PASS_46_RUN_METADATA_ONLY", stats, (
        "Run aliases are not proof of unique donor or 1h/6h/48h exposure; "
        + str(46 - stats["timepoint_unique_text_matches"]) + " run timepoints unverified")


def locate_stage1():
    candidates = []
    for root in ("data", "results", "research", "docs"):
        path = pathlib.Path(root)
        if not path.exists():
            continue
        for p in path.rglob("*.tsv"):
            if "stage1" not in p.name.lower() or p.stat().st_size > 2_000_000:
                continue
            with p.open(encoding="utf-8-sig", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if not {"rsid", "study_id", "gene"}.issubset(set(reader.fieldnames or [])):
                    continue
                rows = list(reader)
            rsids = sorted({r["rsid"] for r in rows if re.fullmatch(r"rs\d+", r.get("rsid", ""))})
            if len(rsids) == 20:
                candidates.append((str(p), rsids))
    if len(candidates) != 1:
        raise Hold("MUSCLE_STAGE1_20_RSID_INPUT_MISSING_OR_AMBIGUOUS_" + str(len(candidates)))
    return candidates[0]


def muscle():
    path, rsids = locate_stage1()
    resolved, absent, errors = 0, 0, 0
    for rsid in rsids:
        url = ("https://gtexportal.org/api/v2/dataset/variant?"
               + urllib.parse.urlencode({"snpId": rsid, "datasetId": "gtex_v10",
                                         "itemsPerPage": 2000}))
        try:
            obj, _ = json_get(url, limit=1_000_000)
            records = obj.get("data", [])
            if not isinstance(records, list):
                errors += 1
            elif records:
                resolved += 1
            else:
                absent += 1
        except (Hold, urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            errors += 1
    return ("PASS_20_VARIANT_API_AUDIT" if errors == 0 else "HOLD_PARTIAL_GTEX_ACCESS"), {
        "input_path": path, "source_rsids": len(rsids),
        "resolved": resolved, "not_in_gtex": absent, "api_errors": errors,
        "dataset": "gtex_v10", "tissue_for_next_stage": "Muscle_Skeletal",
    }, "GTEx variant availability is not eQTL, EAS LD, or colocalization evidence"


def oral():
    text, meta = geo_source("GSE103322")
    sample_links = re.findall(r"^!Series_sample_id\s*=\s*(\S+)", text, re.M)
    return "PASS_GEO_ACCESSION_ONLY", {
        "accession": "GSE103322", "soft_sha256": meta["sha256"],
        "series_sample_links": len(sample_links),
        "unique_sample_links": len(set(sample_links)),
        "explicit_donor_tokens": len(re.findall(r"\b(?:patient|donor)[ _-]?\d+\b", text, re.I)),
    }, "Series-level SOFT cannot establish 18 patient-to-cell expression mappings"


def skin():
    readme_url = "https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt"
    body, meta = fetch(readme_url, limit=500_000)
    readme = body.decode("utf-8", errors="replace")
    expected = ("ALLELE1", "ALLELE0", "A1FREQ", "LOG_OR", "LOG_SE")
    missing = [x for x in expected if x not in readme.upper()]
    if missing:
        raise IntegrityFail("BRISTOL_README_REQUIRED_FIELDS_MISSING_" + "_".join(missing))
    url = ("https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/"
           "perceived_age_summary_results.txt")
    header = None
    try:
        part, m = fetch(url, limit=8192, headers={"Range": "bytes=0-8191",
                                                   "Accept-Encoding": "identity"})
        header = {"bytes": len(part), "sha256_prefix": m["sha256"],
                  "http_status": m["http_status"],
                  "first_line_columns": part.splitlines()[0].decode("utf-8", errors="replace")[:240]
                  if part else None}
    except (Hold, urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
        header = {"status": "HOLD_PARTIAL_SOURCE_ACCESS_" + type(exc).__name__}
    return "PASS_README_SCHEMA_ONLY", {
        "readme_sha256": meta["sha256"],
        "required_allele_effect_fields": list(expected),
        "build_explicit_in_readme": any(x in readme.upper() for x in ("GRCH37", "GRCH38", "HG19", "HG38")),
        "bounded_header_probe": header,
        "full_gwas_file_downloaded": False,
    }, "Genome build, full SNP rows, allele alignment and genome-wide QC not verified"


PROBES = {"is": official_gwas, "glaucoma": glaucoma, "hcc": hcc,
          "radiomics": radiomics, "ocular": ocular, "muscle": muscle,
          "oral": oral, "skin": skin}


def run_one(track: str, outdir: pathlib.Path) -> int:
    task, next_gate = TRACKS[track]
    result = {
        "schema": "master-eight-track-live-source-gate-20261102-v1",
        "date_kst": "2026-11-02",
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "track": track, "selected_task": task, "next_scientific_gate": next_gate,
        "scientific_gate": "HOLD", "no_biological_inference": True,
        "execution_mode": "GITHUB_ACTIONS_ONLY", "ckd": "HOLD_NO_EXECUTION",
    }
    try:
        status, evidence, blocker = PROBES[track]()
        result.update(probe_status=status, evidence=evidence, blocker=blocker)
    except IntegrityFail as exc:
        result.update(probe_status="FAIL_SOURCE_INTEGRITY", evidence={},
                      blocker=str(exc))
    except Hold as exc:
        result.update(probe_status="HOLD_SOURCE_OR_SCHEMA", evidence={},
                      blocker=str(exc))
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError,
            OSError, ValueError, json.JSONDecodeError) as exc:
        result.update(probe_status="HOLD_SOURCE_ACCESS", evidence={},
                      blocker=type(exc).__name__ + ": " + str(exc)[:160])
    except Exception as exc:
        result.update(probe_status="FAIL_PROBE_IMPLEMENTATION", evidence={},
                      blocker=type(exc).__name__ + ": " + str(exc)[:160])
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / (track + "_20261102.json")).write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("track", "probe_status", "scientific_gate", "blocker")},
                     ensure_ascii=False), flush=True)
    return 1 if result["probe_status"].startswith("FAIL") else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--track", choices=tuple(TRACKS), required=True)
    parser.add_argument("--outdir", default="results/research_operations/20261102")
    args = parser.parse_args()
    raise SystemExit(run_one(args.track, pathlib.Path(args.outdir)))


if __name__ == "__main__":
    main()
