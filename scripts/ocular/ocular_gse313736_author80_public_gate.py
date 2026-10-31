#!/usr/bin/env python3
"""Source-bounded audit of GSE313736 author's DEG80 evidence.
Do not confuse a source-author result with the 83-gene independent reanalysis.
No personal donor/genotype information is accessed or published.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ACCESSION = "GSE313736"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313736/suppl/"
QUERIES = [
    'EXT_ID:GSE313736',
    '"GSE313736"',
    '"Thomas Guinan" AND "trabecular meshwork"',
]
ID_PATTERN = re.compile(r"^ENSG\d{11}(?:\.\d+)?$")
DATA_EXT = re.compile(r"\.(?:csv|tsv|txt|xlsx|xls|rds|rdata)(?:\.gz)?$", re.I)
CANDIDATE = re.compile(r"(?:deg|differential|fold.?change|dex|results|supp|table)", re.I)


def fetch(url, maximum=2_000_000):
    last = None
    for attempt in range(2):
        try:
            with urlopen(Request(url, headers={
                "User-Agent": "GSE313736-author80-provenance-audit/1.0 (research)",
                "Accept": "application/json, text/html, */*",
            }), timeout=25) as response:
                body = response.read(maximum + 1)
                if len(body) > maximum:
                    raise ValueError("SOURCE_SIZE_LIMIT")
                return body, response.url
        except Exception as exc:
            last = f"{type(exc).__name__}: {str(exc)[:180]}"
            if attempt == 0:
                time.sleep(1)
    raise RuntimeError(last)


def check_public():
    out = {"accession": ACCESSION, "checked_utc": datetime.now(timezone.utc).isoformat(),
           "geo_series": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736",
           "publication_search": [], "geo_series_supplement_inventory": {},
           "author_80_original_gene_table": "HOLD_NOT_AUTHENTICATED",
           "source_policy": "No author-gene claims from publication hits, our reanalysis, or file names alone"}
    successful = 0
    for query in QUERIES:
        url = EPMC + "?" + urlencode({"query": query, "format": "json", "pageSize": "15"})
        record = {"query": query, "url": url, "status": "HOLD_NOT_RETRIEVED"}
        try:
            body, final = fetch(url)
            data = json.loads(body)
            hits = data.get("resultList", {}).get("result", [])
            record.update(status="CHECKED", total_hit_count=data.get("hitCount", 0),
                          returned_count=len(hits), final_url=final,
                          hits=[{"id": h.get("id"), "doi": h.get("doi"),
                                 "title": h.get("title"), "authorString": h.get("authorString"),
                                 "firstPublicationDate": h.get("firstPublicationDate")}
                                for h in hits[:15]])
            successful += 1
        except Exception as exc:
            record["error"] = str(exc)[:250]
        out["publication_search"].append(record)
    inv = {"url": GEO_FTP, "status": "HOLD_NOT_RETRIEVED"}
    try:
        body, final = fetch(GEO_FTP)
        txt = body.decode("utf8", "replace")
        names = re.findall(r'<a\s+href="([^"]+)"', txt, re.I)
        files = sorted({x.rsplit("/", 1)[-1].split("?")[0] for x in names
                        if x not in ("../", "./") and not x.startswith("?")})
        # Directory listing is not a biological/author-certified result set.
        candidates = [x for x in files if CANDIDATE.search(x) and DATA_EXT.search(x)]
        inv.update(status="CHECKED", final_url=final,
                   listing_sha256=hashlib.sha256(body).hexdigest(),
                   file_count=len(files), candidate_table_filenames=candidates,
                   candidate_table_count=len(candidates))
        successful += 1
    except Exception as exc:
        inv["error"] = str(exc)[:250]
    out["geo_series_supplement_inventory"] = inv
    out["public_source_retrieval_status"] = (
        "PASS_AT_LEAST_ONE_AUTHORITY_CHECKED" if successful else "REVIEW_ALL_SOURCES_UNAVAILABLE")
    return out


def standard_gene(gene):
    gene = str(gene).strip()
    if not ID_PATTERN.fullmatch(gene):
        raise ValueError("INVALID_STABLE_ENSEMBL_ID")
    return gene.split(".")[0]


def read_gene_table(path, *, source_author):
    raw = Path(path).read_text(encoding="utf-8-sig")
    first = raw.splitlines()[0] if raw.splitlines() else ""
    delim = "\t" if first.count("\t") > first.count(",") else ","
    rows = list(csv.DictReader(io.StringIO(raw), delimiter=delim))
    if not rows:
        raise ValueError("EMPTY_TABLE")
    available = set(rows[0])
    def pick(*names):
        return next((n for n in names if n in available), None)
    gene_col = pick("ensembl_gene_id", "gene_id", "GeneID", "ensembl_id")
    fc_col = pick("log2FC", "logFC", "log2FoldChange", "log2fc")
    fdr_col = pick("FDR", "padj", "adj.P.Val", "fdr")
    if not gene_col or not fc_col or (source_author and not fdr_col):
        raise ValueError("MISSING_GENE_FC_OR_AUTHOR_FDR_COLUMNS")
    genes = {}
    for idx, row in enumerate(rows, 2):
        gene = standard_gene(row[gene_col])
        if gene in genes:
            raise ValueError(f"DUPLICATE_GENE_ID_AT_LINE_{idx}")
        v = float(row[fc_col])
        if not math.isfinite(v) or v == 0:
            raise ValueError(f"INVALID_NONZERO_LOG2FC_AT_LINE_{idx}")
        q = None
        if fdr_col and row[fdr_col].strip():
            q = float(row[fdr_col])
            if not math.isfinite(q) or not (0 <= q <= 1):
                raise ValueError(f"INVALID_FDR_AT_LINE_{idx}")
        elif source_author:
            raise ValueError(f"MISSING_AUTHOR_FDR_AT_LINE_{idx}")
        genes[gene] = {"log2FC": v, "FDR": q}
    return genes


def compare(source_path, ours_path, source_url, source_method):
    if not source_url or not source_url.startswith("https://") or not source_method:
        return {"status": "HOLD_SOURCE_PROVENANCE_OR_METHOD_UNSPECIFIED",
                "note": "A local filename is not proof of source-author ownership."}
    try:
        author = read_gene_table(source_path, source_author=True)
        ours = read_gene_table(ours_path, source_author=False)
        author_up = sum(z["log2FC"] > 0 for z in author.values())
        author_down = sum(z["log2FC"] < 0 for z in author.values())
        structural = len(author) == 80 and author_up == 63 and author_down == 17
        shared = sorted(set(author) & set(ours))
        concordant = sum(author[g]["log2FC"] * ours[g]["log2FC"] > 0 for g in shared)
        return {
            "status": "REVIEW_AUTHOR_FILE_PROVENANCE_AND_THRESHOLD" if structural else
                      "HOLD_AUTHOR_COUNTS_DO_NOT_MATCH_63_17",
            "author_table_source_url_claimed_not_independently_verified": source_url,
            "author_method_claimed_not_independently_verified": source_method,
            "author_table_sha256": hashlib.sha256(Path(source_path).read_bytes()).hexdigest(),
            "our_table_sha256": hashlib.sha256(Path(ours_path).read_bytes()).hexdigest(),
            "author_gene_count": len(author), "author_up": author_up, "author_down": author_down,
            "our_gene_count": len(ours), "shared_genes": len(shared),
            "same_direction": concordant, "opposite_direction": len(shared) - concordant,
            "coverage_of_author_80": len(shared) / 80,
            "direction_concordance_among_shared": concordant / len(shared) if shared else None,
            "caution": "Descriptive overlap only. Check published method, source identity and donor/library effects before GO.",
        }
    except (OSError, ValueError, KeyError) as exc:
        return {"status": "HOLD_INPUT_INVALID",
                "error": f"{type(exc).__name__}: {str(exc)[:200]}"}


def self_test():
    assert standard_gene("ENSG00000123456.12") == "ENSG00000123456"
    for bad in ["ABO", "ENSG1", "ENSG00000123456.b", ""]:
        try:
            standard_gene(bad)
            raise AssertionError(f"accepted {bad}")
        except ValueError:
            pass
    import tempfile
    with tempfile.TemporaryDirectory() as temp:
        base = Path(temp)
        src = base / "source.tsv"
        ours = base / "ours.tsv"
        names = [f"ENSG{i:011d}" for i in range(1, 81)]
        src.write_text("gene_id\tlog2FC\tFDR\n" + "".join(
            f"{gene}\t{4 if idx < 63 else -4}\t0.01\n"
            for idx, gene in enumerate(names)), encoding="utf8")
        ours.write_text("gene_id\tlogFC\tFDR\n" + "".join(
            f"{gene}\t{4 if idx < 62 else -4}\t0.02\n"
            for idx, gene in enumerate(names[:79])), encoding="utf8")
        a = compare(src, ours, "https://example.org/source", "specified contrast")
        assert a["status"] == "REVIEW_AUTHOR_FILE_PROVENANCE_AND_THRESHOLD", a
        assert (a["shared_genes"], a["same_direction"], a["opposite_direction"]) == (79, 78, 1)
        b = compare(src, ours, "", "")
        assert b["status"] == "HOLD_SOURCE_PROVENANCE_OR_METHOD_UNSPECIFIED"
        src.write_text(src.read_text() + f"{names[0]}\t3\t0.01\n", encoding="utf8")
        c = compare(src, ours, "https://example.org/source", "specified contrast")
        assert c["status"] == "HOLD_INPUT_INVALID"
    print("AUTHOR80_GATE_SELF_TEST_PASS: ID, duplicate, provenance, 63/17, overlap", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--author-table", type=Path)
    p.add_argument("--ours-table", type=Path)
    p.add_argument("--author-source-url", default="")
    p.add_argument("--author-method", default="")
    a = p.parse_args()
    if a.self_test:
        self_test()
        return
    a.out.mkdir(parents=True, exist_ok=True)
    result = check_public()
    if bool(a.author_table) != bool(a.ours_table):
        result["author_vs_our_gene_overlap"] = {"status": "HOLD_BOTH_GENE_TABLES_REQUIRED"}
    elif a.author_table and a.ours_table:
        result["author_vs_our_gene_overlap"] = compare(
            a.author_table, a.ours_table, a.author_source_url, a.author_method)
    else:
        result["author_vs_our_gene_overlap"] = {
            "status": "HOLD_ORIGINAL_AUTHOR_80_GENE_TABLE_NOT_SUPPLIED",
            "our_83_not_original_80": True,
            "note": "The existing 83 and 88 sets are independent source-cohort reanalysis outputs."}
    output = a.out / "gse313736_author80_public_discovery_gate.json"
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    print("GSE313736_AUTHOR80_DISCOVERY",
          result["public_source_retrieval_status"],
          "author_overlap", result["author_vs_our_gene_overlap"]["status"], flush=True)
    for r in result["publication_search"]:
        print("EPMC", r["query"], r["status"], r.get("total_hit_count"), flush=True)
    print("GEO_SUPPL", result["geo_series_supplement_inventory"]["status"],
          result["geo_series_supplement_inventory"].get("candidate_table_filenames", []), flush=True)
    if result["public_source_retrieval_status"] == "REVIEW_ALL_SOURCES_UNAVAILABLE":
        sys.exit("No public official endpoint could be fetched; review network access, do not infer absence.")


if __name__ == "__main__":
    main()
