#!/usr/bin/env python3
"""Bounded paper-feasibility analyses. Inputs stay outside Git; only derived results are emitted."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
import zlib
from pathlib import Path

POOL_SHA = "12c3c94c1f1346cf72eda7ea9b0d6b102e66bde326040389e7aeffaf58c472c3"
CKD_SHA = "7fcdfdb70f689a149beff7a583c98f2b8e52cfaf831707829427b1d1f6db4424"
ORAL_URLS = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103322&format=file&file=GSE103322_HNSCC_all_data.txt.gz",
]
SKIN_IDS = ("GCST90094903", "GCST90094904")
THRESHOLDS = (0.0, 0.25, 0.5, 1.0, 2.0, 3.0)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def probability(value):
    p = float(value)
    if not math.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("Invalid probability")
    return p


def number(value):
    x = float(value)
    if not math.isfinite(x):
        raise ValueError("Nonfinite effect")
    return x


def same_sign(a, b):
    return a * b > 0


def read_gene_table(z, member, effect, fdr, expected):
    with z.open(member) as raw, gzip.GzipFile(fileobj=raw) as gz:
        reader = csv.DictReader(io.TextIOWrapper(gz, encoding="utf-8"), delimiter="\t")
        if not {"gene_id", "gene_symbol", effect, fdr} <= set(reader.fieldnames or []):
            raise ValueError("Missing table fields")
        out = {}
        for r in reader:
            gid = r["gene_id"].split(".")[0]
            if not re.fullmatch(r"ENSG\d+", gid) or gid in out:
                raise ValueError("Invalid/duplicated gene ID")
            out[gid] = (number(r[effect]), probability(r[fdr]), r["gene_symbol"])
    if len(out) != expected:
        raise ValueError("Unexpected tested-gene count")
    return out


def ocular_summary(donor, pool, external):
    shared = set(donor) & set(pool) & set(external)
    rows = []
    for threshold in THRESHOLDS:
        eligible = [g for g in sorted(shared) if donor[g][1] < 0.05 and abs(donor[g][0]) >= threshold]
        stable = [g for g in eligible if pool[g][1] < 0.05 and abs(pool[g][0]) >= threshold
                  and same_sign(donor[g][0], pool[g][0])]
        def tally(gs):
            known = [g for g in gs if external[g]["fdr"] is not None]
            signed = [g for g in known if same_sign(donor[g][0], external[g]["log2fc"])]
            supported = [g for g in signed if external[g]["fdr"] < 0.05]
            opposite = [g for g in known if external[g]["fdr"] < 0.05
                        and donor[g][0] * external[g]["log2fc"] < 0]
            return {"n": len(gs), "external_fdr_available": len(known),
                    "external_fdr_missing": len(gs)-len(known),
                    "external_same_sign": len(signed),
                    "external_fdr05_same_sign": len(supported),
                    "external_fdr05_opposite_sign": len(opposite)}
        rows.append({"abs_log2fc_threshold": threshold,
                     "donor_fdr05": tally(eligible),
                     "donor_and_pool_fdr05_same_sign": tally(stable)})
    strongest = [g for g in sorted(shared) if donor[g][1] < .05 and pool[g][1] < .05
                 and abs(donor[g][0]) >= 1 and abs(pool[g][0]) >= 1
                 and same_sign(donor[g][0], pool[g][0])
                 and external[g]["fdr"] is not None and external[g]["fdr"] < .05
                 and same_sign(donor[g][0], external[g]["log2fc"])]
    return {"status": "DESCRIPTIVE_THREE_MODEL_COMPARISON",
            "donor_n": len(donor), "pool_n": len(pool), "external_n": len(external),
            "shared_gene_ids": len(shared),
            "thresholds": rows, "three_model_abs1_n": len(strongest),
            "interpretation": "Donor/pool are same study; external published TM estimates differ in dose and protocol. Thresholds describe already inspected outcomes, not a new independent confirmation or a statistical test across genes. No causal, clinical, novelty, or effect-equality claim."}


def ocular(args):
    if sha(args.pool_zip) != POOL_SHA:
        raise ValueError("Pinned P5d ZIP mismatch")
    with zipfile.ZipFile(args.pool_zip) as z:
        if z.testzip() is not None:
            raise ValueError("ZIP integrity error")
        donor = read_gene_table(z, "p5_all_gene_results.tsv.gz", "logFC", "FDR", 16050)
        pool = read_gene_table(z, "p5d_poolunit_all_gene_results.tsv.gz", "pool_logFC", "pool_FDR", 21081)
    sys.path.insert(0, str(args.ocular_source / "scripts/ocular"))
    from ocular_p5j_official_TM_overlap import download_zip, parse_official, TM_ALL
    with tempfile.TemporaryDirectory() as d:
        source = Path(d) / "supp.zip"
        download_zip(source)
        external = parse_official(source)[TM_ALL]
        external_sha = sha(source)
    if len(external) != 13313:
        raise ValueError("Published TM table changed; verify erratum/version")
    out = ocular_summary(donor, pool, external)
    if out["shared_gene_ids"] < 12000:
        raise ValueError("Unexpected cross-study gene intersection")
    out.update({"pool_artifact_id": 12995175346, "pool_zip_sha256": POOL_SHA,
                "external_supplement_sha256": external_sha,
                "external_parser_commit": "0d39a9b41aa500d39327fdd8709f9ccb179f8502",
                "external_source": "PMC12882097 iovs-67-2-13_s003.xlsx SuppTable4_TM2Days_All_DEGs"})
    return out


def ckd_summary(rows):
    ids = [r["protein_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicated protein assay")
    # This is a retrospective re-screen: EAS result/direction is NOT a selection rule.
    discovery = [r for r in rows if probability(r["eur_egfr_fdr"]) < .05
                 and int(r["support_same_risk_direction_n"]) >= 3
                 and int(r["support_p05_same_risk_direction_n"]) >= 2]
    n = len(discovery)
    result = []
    for r in discovery:
        item = {"protein_id": r["protein_id"], "gene": r["gene_symbol"],
                "anchor": r["eur_egfr_rsid"], "eur_beta": number(r["eur_egfr_beta"])}
        valid = (r["eas_egfr_available"] == "1" and r["eas_same_top_snp"] == "1"
                 and r["eas_egfr_rsid"] == r["eur_egfr_rsid"]
                 and r["eas_egfr_p"] != "" and r["eas_egfr_beta"] != "")
        if valid:
            p = probability(r["eas_egfr_p"])
            b = number(r["eas_egfr_beta"])
            item.update({"eas_beta": b, "eas_p": p, "eas_bonferroni": min(1.0, n*p),
                         "direction_concordant": same_sign(item["eur_beta"], b),
                         "passes": n*p < .05 and same_sign(item["eur_beta"], b)})
        else:
            item.update({"eas_beta": None, "eas_p": None, "eas_bonferroni": None,
                         "direction_concordant": None, "passes": False})
        result.append(item)
    result.sort(key=lambda x: (x["eas_p"] is None, x["eas_p"] if x["eas_p"] is not None else 1, x["protein_id"]))
    return {"status": "RETROSPECTIVE_EUR_ONLY_SELECTION",
            "all_assays": len(rows), "eur_selected_family": n,
            "eas_same_anchor_tested": sum(x["eas_p"] is not None for x in result),
            "eas_unavailable": sum(x["eas_p"] is None for x in result),
            "eas_bonferroni_direction_concordant": [x["gene"] for x in result if x["passes"]],
            "candidates": result, "interpretation": "Bonferroni uses every EUR-selected assay, including unavailable EAS tests. EAS outcome results were already inspected in prior work; shared EUR pQTL instrument and cohort overlap not ruled out. Not independent exposure replication, MR refit or colocalization."}


def ckd(args):
    if sha(args.ckd_table) != CKD_SHA:
        raise ValueError("Pinned September aggregate source mismatch")
    with args.ckd_table.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 1506:
        raise ValueError("Expected 1506 original assays")
    result = ckd_summary(rows)
    result["source_sha256"] = CKD_SHA
    return result


def get_prefix(url, limit=262144, range_request=True):
    headers = {"User-Agent": "publication-feasibility/1.0", "Accept-Encoding": "identity"}
    if range_request:
        headers["Range"] = f"bytes=0-{limit-1}"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=25) as r:
            data = r.read(limit)
            return {"url": url, "http_status": r.status, "bytes_read": len(data),
                    "content_range": r.headers.get("Content-Range"),
                    "content_type": r.headers.get("Content-Type")}, data
    except urllib.error.HTTPError as e:
        return {"url": url, "http_status": e.code, "bytes_read": 0}, b""
    except (OSError, ValueError) as e:
        return {"url": url, "error_type": type(e).__name__, "bytes_read": 0}, b""


def inspect_vcf_prefix(data):
    if not data.startswith(b"\x1f\x8b"):
        return {"gzip": False, "vcf_header": False}
    try:
        decoded = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data, 2 * 1024 * 1024)
        # BGZF may consist of multiple concatenated members: first member is enough for
        # a header if complete; otherwise this stays an incomplete prefix, not a negative result.
        lines = decoded.decode("utf-8", errors="replace").splitlines()
        fields = sorted(set(re.findall(r"^##FORMAT=<ID=([^,>]+)", "\n".join(lines), re.M)))
        return {"gzip": True, "vcf_header": any(x.startswith("##fileformat=VCF") for x in lines),
                "format_ids": fields, "has_beta_se": {"ES", "SE"} <= set(fields),
                "has_variant_column_header": any(x.startswith("#CHROM\t") for x in lines),
                "reference_headers": [x[:200] for x in lines if x.startswith("##reference=")][:3]}
    except zlib.error:
        return {"gzip": True, "vcf_header": False, "prefix_decode_error": True}


def skin(args):
    results = []
    for acc in SKIN_IDS:
        ebi = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/" + acc + "/"
        directory, body = get_prefix(ebi, 131072, False)
        if directory.get("http_status") == 200:
            links = re.findall(r'href="([^"]+)"', body.decode("utf-8", errors="replace"))
            directory["summary_file_links"] = [x for x in links if re.search(r"\.(?:tsv|txt|vcf)(?:\.gz)?$", x, re.I)][:20]
        opengwas = "https://gwas.mrcieu.ac.uk/files/ebi-a-" + acc + "/ebi-a-" + acc + ".vcf.gz"
        status, prefix = get_prefix(opengwas)
        status.update(inspect_vcf_prefix(prefix))
        results.append({"accession": acc, "ebi_ftp": directory, "opengwas_vcf": status})
    return {"status": "BOUNDED_SOURCE_ACCESS_ATTEMPT", "sources": results,
            "interpretation": "Only actual bounded HTTP reads; a missing/blocked prefix is not data absence. A VCF header does not establish full-file, allele, phenotype or instrument validity; no MR executed."}


def oral_matrix(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader)
        # Some public GEO text matrices have a blank trailing delimiter.
        if first and first[-1] == "":
            first = first[:-1]
            trailing = True
        else:
            trailing = False
        n = len(first) - 1
        if n != 5902 or len(set(first[1:])) != n:
            raise ValueError("Unexpected or duplicate cell columns")
        labels = {}
        numeric_rows = 0
        invalid_width = 0
        partial_numeric_check = 0
        for row in reader:
            if trailing and row and row[-1] == "":
                row = row[:-1]
            if len(row) != n + 1:
                invalid_width += 1
                continue
            key = row[0].strip().lower()
            # Retain only explicit metadata rows; do not infer patient identity from column names.
            if re.fullmatch(r"(?:patient(?:[ _-]?id)?|tumou?r(?:[ _-]?(?:id|site|type))?|sample(?:[ _-]?(?:id|site|type))?|site|tissue|cell[ _-]?type|malignant(?:[ _-]?status)?|malignancy|classification|cancer(?:[ _-]?cell)?|non[ _-]?cancer[ _-]?cell[ _-]?type|enzyme|rt[ _-]?enzyme)", key):
                labels[key] = row[1:]
                continue
            try:
                vals = [float(row[i]) for i in (1, n//2, n)]
                if any(not math.isfinite(x) or x < 0 for x in vals):
                    continue
                numeric_rows += 1
                partial_numeric_check += len(vals)
            except (ValueError, IndexError):
                # Unknown metadata row: keep its existence without leaking cell-level values.
                continue
        fields = {key: {"unique_nonempty_values": len({x for x in vals if x.strip()}),
                        "nonempty_cells": sum(bool(x.strip()) for x in vals)}
                  for key, vals in labels.items()}
        patients = next((v for k,v in labels.items() if re.fullmatch(r"patient(?:[ _-]?id)?", k)), None)
        sites = next((v for k,v in labels.items() if k in {"site", "sample site", "sample_site", "tumor site", "tumour site", "tumor_site", "tumour_site", "tissue"}), None)
        paired = None
        if patients is not None and sites is not None:
            by = collections.defaultdict(set)
            for p,s in zip(patients,sites):
                if p.strip():
                    s = s.strip().lower().replace("_"," ").replace("-"," ")
                    if s in {"primary", "primary tumor", "primary tumour"}:
                        by[p].add("primary")
                    elif s in {"ln", "lymph node", "lymph node metastasis", "metastatic lymph node"}:
                        by[p].add("ln")
            paired = sum(v == {"primary","ln"} for v in by.values())
        return {"cells": n, "numeric_expression_rows_sample_checked": numeric_rows,
                "numeric_values_checked": partial_numeric_check, "inconsistent_width_rows": invalid_width,
                "explicit_metadata_fields": fields, "explicit_primary_ln_paired_patients": paired,
                "interpretation": "Public processed TPM; checked three values per expression row, not all entries. Counts describe file structure, not differential expression. No donor identifiers or cell expression exported."}


def oral(args):
    attempts = []
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "processed.tsv.gz"
        for url in ORAL_URLS:
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent":"publication-feasibility/1.0","Accept-Encoding":"identity"}), timeout=35) as response, path.open("wb") as out:
                    if response.status != 200:
                        raise ValueError("Expected complete response")
                    size = 0
                    while True:
                        block = response.read(1024*1024)
                        if not block:
                            break
                        size += len(block)
                        if size > 120*1024*1024:
                            raise ValueError("Transfer limit")
                        out.write(block)
                result = oral_matrix(path)
                result.update({"status": "PROCESSED_TPM_STRUCTURE_READ", "url": url,
                               "bytes": size, "sha256": sha(path), "earlier_attempts": attempts})
                return result
            except urllib.error.HTTPError as e:
                attempts.append({"url": url, "http_status": e.code})
            except (OSError, ValueError, StopIteration, UnicodeError, csv.Error) as e:
                attempts.append({"url": url, "error_type": type(e).__name__})
            finally:
                path.unlink(missing_ok=True)
    return {"status": "SOURCE_TRANSFER_OR_STRUCTURE_UNRESOLVED", "attempts": attempts,
            "interpretation": "No patient-level inference attempted. Use official ExperimentHub EH5419 if direct GEO remains inaccessible."}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("analysis", choices=["ocular","ckd","skin","oral"])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--pool-zip", type=Path)
    p.add_argument("--ocular-source", type=Path)
    p.add_argument("--ckd-table", type=Path)
    args = p.parse_args()
    result = globals()[args.analysis](args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
