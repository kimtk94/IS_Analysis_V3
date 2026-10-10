#!/usr/bin/env python3
"""Conservative GRCh37-to-GRCh38 rsID mapping for TPMI 433.21.
Ensembl Variation GRCh38 mapping and GRCh38 sequence independently checked.
No TPMI GWAS downloaded; no automatic allelic complementation/liftover.
"""
import argparse
import csv
import json
import time
import urllib.request
from collections import Counter
from pathlib import Path

HOST = "https://rest.ensembl.org"

def get_json(path):
    req = urllib.request.Request(HOST + path, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=25) as response:
        return json.load(response)

def resolve(target, variation, sequence):
    chrom = target["chrom"]
    rsid = target["allele_verified_rsids"]
    if variation.get("name") != rsid:
        return "RSID_MISMATCH", ""
    maps = [m for m in variation.get("mappings", [])
            if m.get("assembly_name") == "GRCh38"
            and str(m.get("seq_region_name")) == chrom]
    if len(maps) != 1:
        return "MAPPING_NOT_UNIQUE", ""
    m = maps[0]
    if m.get("strand") != 1 or not isinstance(m.get("start"), int) or m.get("end") != m["start"]:
        return "UNSUPPORTED_STRAND_OR_INDEL", ""
    ref = str(sequence.get("seq", "")).upper()
    if "GRCh38" not in str(sequence.get("id", "")):
        return "SEQUENCE_BUILD_UNVERIFIED", ""
    possible = set(str(m.get("allele_string", "")).upper().split("/"))
    original = {target["ref"].upper(), target["alt"].upper()}
    if not original.issubset(possible) or ref not in original:
        return "SOURCE_ALLELES_NOT_SUPPORTED_GRCH38", ""
    if any(len(a) != 1 or a not in "ACGT" for a in original | {ref}):
        return "NON_SNV_ALLELE", ""
    alt = (original - {ref}).pop()
    return "VERIFIED_GRCH38_REF_ALT", f"{chrom}:{m['start']}:{ref}:{alt}"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    with args.input.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 11 or len({r["variant_id"] for r in rows}) != 11:
        raise ValueError("Require exactly 11 unique audited GRCh37 variants")
    if any(r["allele_mapping_status"] != "VERIFIED_ALLELE_SET" for r in rows):
        raise ValueError("GRCh37 VEP allele verification incomplete")
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError("Refuse to overwrite output")
    output = []
    for r in rows:
        rsid, chrom = r["allele_verified_rsids"], r["chrom"]
        item = {"locus": r["locus"], "rsid": rsid, "GRCh37": r["variant_id"],
                "GRCh38": "", "status": "NOT_QUERIED", "error": ""}
        try:
            if not rsid.startswith("rs") or ";" in rsid:
                raise ValueError("Nonunique rsID")
            v = get_json(f"/variation/human/{rsid}?content-type=application/json")
            maps = [m for m in v.get("mappings", [])
                    if m.get("assembly_name") == "GRCh38"
                    and str(m.get("seq_region_name")) == chrom]
            if len(maps) == 1 and isinstance(maps[0].get("start"), int):
                pos = maps[0]["start"]
                seq = get_json(f"/sequence/region/human/{chrom}:{pos}..{pos}:1?content-type=application/json")
            else:
                seq = {}
            item["status"], item["GRCh38"] = resolve(r, v, seq)
        except (OSError, ValueError, KeyError, TimeoutError) as exc:
            item["status"], item["error"] = "API_OR_INPUT_ERROR", str(exc)[:160]
        output.append(item)
        print("MAP", rsid, item["status"], item["GRCh38"], flush=True)
        time.sleep(0.16)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir/"TPMI_11_SNPS_GRCH38_MAPPING.tsv").open("w", newline="") as f:
        cols = ["locus", "rsid", "GRCh37", "GRCh38", "status", "error"]
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    summary = {"status": "REFERENCE_MAPPING_ONLY_NOT_TPMI_GWAS",
               "n_snps": len(output),
               "verified": sum(x["status"]=="VERIFIED_GRCH38_REF_ALT" for x in output),
               "unresolved": sum(x["status"]!="VERIFIED_GRCH38_REF_ALT" for x in output),
               "statuses": dict(Counter(x["status"] for x in output)),
               "TPMI_phenotype": "433.21",
               "TPMI_GWAS_analyzed": False,
               "fine_mapping_ALDH2": "BLOCKED", "fine_mapping_ADH1B": "EXPLORATORY"}
    (args.out_dir/"TPMI_GRCH38_MAPPING_MANIFEST.json").write_text(json.dumps(summary, indent=2)+"\n")
    print("TPMI_GRCH38_MAPPING_COMPLETE", json.dumps(summary), flush=True)

if __name__ == "__main__":
    main()
