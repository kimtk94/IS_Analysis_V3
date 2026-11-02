#!/usr/bin/env python3
"""Official GWAS Catalog v2 metadata check for five GIGASTROKE EAS phenotypes.

Public, bounded GET requests only. No summary-statistics download, SSH,
patient IDs, effect allele inference, adoption of case/control N, or coloc.
"""
import argparse
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"
EXPECTED = {
    "GCST90104544": {"phenotype": "any_stroke", "cases": 27413},
    "GCST90104545": {"phenotype": "any_ischemic_stroke", "cases": 19032},
    "GCST90104546": {"phenotype": "cardioembolic_stroke", "cases": 926},
    "GCST90104547": {"phenotype": "large_artery_stroke", "cases": 1735},
    "GCST90104548": {"phenotype": "small_vessel_stroke", "cases": 5532},
}
CONTROLS_PUBLISHED = 237242
MAX_BYTES = 1024 * 1024
IDS = ("accession_id", "accessionId", "accessionID", "accession", "study_accession")
TRAIT_KEYS = ("disease_trait", "diseaseTrait", "reported_trait", "reportedTrait")
ANCESTRY_KEYS = ("ancestral_group", "ancestralGroup", "ancestry", "ancestry_group",
                 "ancestryGroup", "ancestral_groups", "ancestralGroups")


def official_url(study, ancestry=False):
    if study not in EXPECTED:
        raise ValueError("Study not allowlisted")
    return BASE + "/" + study + ("/ancestries" if ancestry else "")


def fetch_json(url, fetcher=None):
    if not url.startswith(BASE + "/"):
        raise ValueError("Unapproved endpoint")
    if fetcher:
        return fetcher(url)
    req = urllib.request.Request(url, headers={
        "Accept": "application/json",
        "User-Agent": "IS-GIGASTROKE-v2-provenance-audit/1.0",
    })
    with urllib.request.urlopen(req, timeout=18) as response:
        final = urllib.parse.urlsplit(response.url)
        if final.scheme != "https" or final.hostname != "www.ebi.ac.uk":
            raise ValueError("Unexpected redirect origin")
        if response.headers.get_content_type() not in ("application/json", "application/hal+json"):
            raise ValueError("Unexpected response MIME")
        raw = response.read(MAX_BYTES + 1)
        if len(raw) > MAX_BYTES:
            raise ValueError("Oversize metadata response")
    return json.loads(raw.decode("utf-8"))


def candidates(document):
    if isinstance(document, dict):
        if any(k in document for k in IDS):
            yield document
        for key in ("data", "results", "items", "studies", "content", "_embedded"):
            if key in document:
                yield from candidates(document[key])
    elif isinstance(document, list):
        for value in document:
            yield from candidates(value)


def exact_study(document, expected):
    matches = []
    for item in candidates(document):
        values = [str(item[k]) for k in IDS if isinstance(item.get(k), (str,int))]
        if expected in values:
            matches.append(item)
    if len(matches) != 1:
        return None
    return matches[0]


def flat_metadata_values(document, keys):
    values = []
    def walk(obj):
        if isinstance(obj, dict):
            for key, val in obj.items():
                if key in keys and isinstance(val, (str,int,float)):
                    values.append(str(val))
                if key in keys and isinstance(val, dict):
                    for name in ("name", "label", "value", "ancestral_group", "ancestralGroup"):
                        if isinstance(val.get(name), str):
                            values.append(val[name])
                if isinstance(val, (dict,list)):
                    walk(val)
        elif isinstance(obj, list):
            for val in obj:
                walk(val)
    walk(document)
    return values


def matching_trait(expected, text):
    label = re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()
    if expected == "any_stroke":
        return label in ("stroke", "any stroke", "all stroke", "stroke any")
    if expected == "any_ischemic_stroke":
        return label in ("any ischemic stroke", "any ischaemic stroke",
                         "ischemic stroke", "ischaemic stroke")
    terms = {
        "cardioembolic_stroke": ("cardioembolic", "cardio embolic"),
        "large_artery_stroke": ("large artery", "large vessel", "atherosclerotic"),
        "small_vessel_stroke": ("small vessel", "small artery", "lacunar"),
    }[expected]
    return any(t in label for t in terms)


def inspect(study, study_doc, ancestry_doc):
    e = EXPECTED[study]
    record = exact_study(study_doc, study)
    groups = flat_metadata_values(ancestry_doc, ANCESTRY_KEYS)
    group_match = any(v.strip().casefold() in ("east asian", "east asian ancestry")
                      for v in groups)
    traits = flat_metadata_values(record, TRAIT_KEYS) if record else []
    trait_match = any(matching_trait(e["phenotype"], v) for v in traits)
    return {
        "accession": study,
        "published_reference_phenotype": e["phenotype"],
        "published_reference_cases_not_adopted": e["cases"],
        "published_reference_controls_not_adopted": CONTROLS_PUBLISHED,
        "published_reference_total_not_adopted": e["cases"] + CONTROLS_PUBLISHED,
        "exact_study_id": record is not None,
        "official_study_trait_match": trait_match if traits else None,
        "official_EAS_ancestry_match": group_match if groups else None,
        "study_metadata_sha256": hashlib.sha256(
            json.dumps(study_doc, sort_keys=True).encode()).hexdigest(),
        "ancestry_metadata_sha256": hashlib.sha256(
            json.dumps(ancestry_doc, sort_keys=True).encode()).hexdigest(),
        "metadata_gate": ("STUDY_AND_EAS_TRAIT_SCREEN_PASS" if
                          record is not None and trait_match and group_match
                          else "HOLD_METADATA_INCOMPLETE_OR_MISMATCH"),
        "effect_allele_definition_verified": False,
        "GWAS_original_file_hash_verified": False,
        "GWAS_SE_N_case_fraction_verified": False,
        "coloc_approved": False,
    }


def run(fetcher=None):
    reports = []
    for study in EXPECTED:
        try:
            studies = fetch_json(official_url(study), fetcher)
            ancestries = fetch_json(official_url(study, ancestry=True), fetcher)
            reports.append(inspect(study, studies, ancestries))
        except (ValueError, OSError, json.JSONDecodeError) as exc:
            reports.append({
                "accession": study,
                "metadata_gate": "HOLD_SOURCE_UNAVAILABLE",
                "failure_class": type(exc).__name__,
                "coloc_approved": False,
            })
    return {
        "schema": "is-gigastroke-five-eas-v2-provenance-v1",
        "official_api": BASE,
        "source_accessions_checked": list(EXPECTED),
        "reported_trait_reference": "Mishra 2022 GIGASTROKE; secondary published EAS subtype sample-size table",
        "accession_level_screen_complete": all(
            r["metadata_gate"] == "STUDY_AND_EAS_TRAIT_SCREEN_PASS" for r in reports),
        "scientific_status": "HOLD",
        "independent_EAS_replication_verified": False,
        "GWAS_original_files_downloaded": False,
        "variant_level_harmonization_performed": False,
        "MR_or_coloc_performed": False,
        "reports": reports,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    a = parser.parse_args()
    report = run()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print("IS_GIGASTROKE_EAS_V2=" + json.dumps(report, separators=(",",":"), sort_keys=True))


if __name__ == "__main__":
    main()
