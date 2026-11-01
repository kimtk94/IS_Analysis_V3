#!/usr/bin/env python3
"""Targeted public HCC quantitative-file discovery; metadata only, no downloads of result data."""
import argparse
import collections
import datetime as dt
import json
import pathlib
import urllib.parse

from probe_quant_matrix_manifest import PDC_IDS, PRIDE_IDS, PRIDE, read_bounded, request_graphql

def pdc_protein_assembly():
    result = {}
    for acc in PDC_IDS:
        rec = {"study": acc, "status": "UNVERIFIED"}
        try:
            q = ('{getPaginatedFiles(pdc_study_id:"' + acc +
                 '" data_category:"Protein Assembly" offset:0 limit:100)'
                 '{total files{file_name file_type file_format data_category file_id}}}')
            d = request_graphql(q)["getPaginatedFiles"]
            if isinstance(d, list):
                d = d[0] if d else {}
            rec.update({"status": "PUBLIC_PROTEIN_ASSEMBLY_FILE_INDEX",
                        "total": d.get("total"),
                        "files": [{k: v.get(k) for k in ("file_name", "file_type", "file_format",
                                                      "data_category", "file_id")}
                                  for v in (d.get("files") or [])[:100]]})
        except Exception as e:
            rec["error"] = str(e)[:240]
        result[acc] = rec
    return result

def pride_pages(max_pages=8):
    result = {}
    for acc in PRIDE_IDS:
        names = []
        category_counts = collections.Counter()
        pages = []
        errors = []
        for page in range(max_pages):
            url = f"{PRIDE}/projects/{acc}/files?pageSize=100&page={page}"
            try:
                answer = json.loads(read_bounded(url, limit=3 * 1024 * 1024))
                if isinstance(answer, dict):
                    files = (answer.get("_embedded") or {}).get("files") or answer.get("files") or answer.get("content") or []
                else:
                    files = answer
                if not isinstance(files, list):
                    raise ValueError("unexpected project files response")
                pages.append({"page": page, "files": len(files)})
                for f in files:
                    if not isinstance(f, dict):
                        continue
                    label = f.get("fileCategory")
                    label = label.get("value") if isinstance(label, dict) else str(label)
                    category_counts[label] += 1
                    filename = f.get("fileName") or f.get("name") or ""
                    if label in ("RESULT", "OTHER") or any(s in filename.lower() for s in
                          ("quant", "result", "protein", "phospho", "matrix", "report", "sdrf", "tsv", "csv", ".xlsx")):
                        names.append({"file_name": filename, "file_size_bytes": f.get("fileSizeBytes"),
                                      "file_category": label})
                if len(files) < 100:
                    break
            except Exception as e:
                errors.append({"page": page, "reason": str(e)[:220]})
                break
        result[acc] = {"pages_fetched": pages, "category_counts_scanned": dict(category_counts),
                       "candidate_count_in_scanned_pages": len(names),
                       "candidate_file_manifest": names[:150],
                       "possibly_truncated": bool(pages and pages[-1]["files"] == 100),
                       "errors": errors,
                       "status": "PARTIAL_METADATA_INDEX" if errors or (pages and pages[-1]["files"] == 100) else "SCANNED_TO_PAGE_END"}
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/hcc/feasibility/quant_candidates.json")
    args = parser.parse_args()
    r = {"checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
         "pdc_protein_assembly": pdc_protein_assembly(),
         "pride_project_files": pride_pages(),
         "scope": "Public filenames/file classifications; NO protein quant values, RNA expression, or sample join."}
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
       "pdc": {k: {"status": v["status"], "total": v.get("total"),
                   "file_names": [x["file_name"] for x in v.get("files", [])[:20]], "error": v.get("error")}
               for k,v in r["pdc_protein_assembly"].items()},
       "pride": {k: {"status": v["status"], "scanned": sum(p["files"] for p in v["pages_fetched"]),
                    "candidate_count": v["candidate_count_in_scanned_pages"],
                    "candidate_names": [x["file_name"] for x in v["candidate_file_manifest"][:40]],
                    "truncated": v["possibly_truncated"], "errors": v["errors"]}
                 for k,v in r["pride_project_files"].items()}
    }, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
