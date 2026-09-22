#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import json
import lzma
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

GENO_EXTS = {
    ".bim": "plink1_bim",
    ".pvar": "plink2_pvar",
    ".vcf": "vcf",
}
PHENO_EXTS = {".csv", ".tsv", ".txt", ".gz", ".bz2", ".xz"}

ROLE_PATTERNS = {
    "participant_id": [r"(^|_)id$", r"subject", r"participant", r"person", r"sample"],
    "visit_time": [r"visit", r"wave", r"follow", r"exam", r"date", r"time", r"year"],
    "egfr": [r"egfr"],
    "creatinine": [r"creatin", r"(^|_)scr($|_)"],
    "ckd": [r"(^|_)ckd($|_)", r"kidney"],
    "age": [r"(^|_)age($|_)"],
    "sex": [r"(^|_)sex($|_)", r"gender"],
    "bmi": [r"(^|_)bmi($|_)"],
    "diabetes": [r"diabet", r"dm"],
    "hypertension": [r"hypert", r"htn", r"blood.?pressure", r"(^|_)sbp($|_)", r"(^|_)dbp($|_)"],
}


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows, fields=None):
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as fh:
        if not fields:
            return
        w = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t",
            lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def norm_chr(x):
    x = str(x or "").strip()
    x = re.sub(r"^chr", "", x, flags=re.I)
    return x.upper()


def norm_allele(x):
    return str(x or "").strip().upper()


def is_palindromic(a, b):
    return {norm_allele(a), norm_allele(b)} in ({"A", "T"}, {"C", "G"})


def open_text(path: Path):
    name = path.name.lower()
    if name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    if name.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    if name.endswith(".xz"):
        return lzma.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    if name.endswith(".zst"):
        if not shutil.which("zstdcat"):
            raise RuntimeError("zstdcat is required for .zst pvar files")
        proc = subprocess.Popen(
            ["zstdcat", str(path)],
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        return proc.stdout
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def load_anchor_panel(path: Path):
    rows = read_tsv(path)
    out = []
    for r in rows:
        x = dict(r)
        x["gene_symbol"] = x["gene_symbol"].upper()
        x["chrom_hg19"] = norm_chr(x.get("chrom_hg19"))
        x["pos_hg19"] = str(x.get("pos_hg19", "")).strip()
        x["ref_allele"] = norm_allele(x.get("ref_allele"))
        x["effect_allele_alt"] = norm_allele(x.get("effect_allele_alt"))
        x["palindromic"] = int(
            is_palindromic(x["ref_allele"], x["effect_allele_alt"])
        )
        out.append(x)
    return out


def discover_genotype_metadata(root: Path):
    found = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        low = p.name.lower()
        kind = None
        if low.endswith(".pvar.zst") or low.endswith(".pvar.gz"):
            kind = "plink2_pvar"
        elif p.suffix.lower() in GENO_EXTS:
            kind = GENO_EXTS[p.suffix.lower()]
        elif low.endswith(".bgen"):
            kind = "bgen"
        elif low.endswith(".pgen"):
            kind = "plink2_pgen"
        elif low.endswith(".bed"):
            kind = "plink1_bed"
        elif low.endswith(".vcf.gz"):
            kind = "vcf"
        if kind:
            found.append({
                "path": str(p),
                "kind": kind,
                "size_bytes": p.stat().st_size,
            })
    found.sort(key=lambda r: (r["kind"], r["path"]))
    return found


def parse_bim(path: Path):
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip():
                continue
            v = re.split(r"\s+", line.strip())
            if len(v) < 6:
                continue
            yield {
                "source": str(path), "source_kind": "plink1_bim",
                "line": line_no, "chrom": norm_chr(v[0]), "id": v[1],
                "pos": v[3], "a1": norm_allele(v[4]),
                "a2": norm_allele(v[5]), "ref": "", "alt": "",
            }


def parse_pvar(path: Path):
    with open_text(path) as fh:
        header = None
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            v = line.rstrip("\n").split("\t")
            if header and len(v) >= len(header):
                r = dict(zip(header, v))
                chrom = r.get("CHROM", "")
                pos = r.get("POS", "")
                vid = r.get("ID", "")
                ref = r.get("REF", "")
                alt = r.get("ALT", "")
            elif len(v) >= 5:
                chrom, pos, vid, ref, alt = v[:5]
            else:
                continue
            yield {
                "source": str(path), "source_kind": "plink2_pvar",
                "line": line_no, "chrom": norm_chr(chrom), "id": vid,
                "pos": str(pos), "ref": norm_allele(ref),
                "alt": norm_allele(alt), "a1": "", "a2": "",
            }


def parse_vcf_with_bcftools(path: Path, anchors):
    if not shutil.which("bcftools"):
        return [], "bcftools_missing"
    regions = sorted({
        (a["chrom_hg19"], a["pos_hg19"]) for a in anchors
        if a["chrom_hg19"] and a["pos_hg19"]
    })
    region_arg = ",".join(f"{c}:{p}-{p}" for c, p in regions)
    cmd = [
        "bcftools", "query", "-r", region_arg,
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\n", str(path)
    ]
    try:
        cp = subprocess.run(
            cmd, check=True, text=True, capture_output=True
        )
    except subprocess.CalledProcessError as exc:
        return [], f"bcftools_exit_{exc.returncode}"
    rows = []
    for i, line in enumerate(cp.stdout.splitlines(), 1):
        v = line.split("\t")
        if len(v) < 5:
            continue
        rows.append({
            "source": str(path), "source_kind": "vcf",
            "line": i, "chrom": norm_chr(v[0]), "pos": v[1],
            "id": v[2], "ref": norm_allele(v[3]),
            "alt": norm_allele(v[4]), "a1": "", "a2": "",
        })
    return rows, "queried"


def allele_status(anchor, rec):
    ref = anchor["ref_allele"]
    alt = anchor["effect_allele_alt"]
    if rec["source_kind"] == "plink1_bim":
        alleles = {rec["a1"], rec["a2"]}
        if {ref, alt}.issubset(alleles):
            return "compatible_unordered"
        return "allele_mismatch"
    rref = rec.get("ref", "")
    alts = set(rec.get("alt", "").split(","))
    if rref == ref and alt in alts:
        return "ref_alt_match"
    if rref == alt and ref in alts:
        return "ref_alt_swapped"
    return "allele_mismatch"


def match_anchors(records, anchors, build_confirmed):
    by_id = defaultdict(list)
    by_pos = defaultdict(list)
    for r in records:
        if r.get("id") and r["id"] != ".":
            by_id[r["id"]].append(r)
        by_pos[(r["chrom"], str(r["pos"]))].append(r)

    out = []
    for a in anchors:
        matches = list(by_id.get(a["rsid"], []))
        match_type = "rsid"
        if not matches:
            matches = list(by_pos.get(
                (a["chrom_hg19"], a["pos_hg19"]), []
            ))
            match_type = (
                "chrpos_hg19_confirmed"
                if build_confirmed else
                "chrpos_hg19_unconfirmed_build"
            )
        if not matches:
            out.append({
                "gene_symbol": a["gene_symbol"], "rsid": a["rsid"],
                "chrom_hg19": a["chrom_hg19"],
                "pos_hg19": a["pos_hg19"],
                "panel_ref": a["ref_allele"],
                "panel_effect_alt": a["effect_allele_alt"],
                "palindromic": a["palindromic"],
                "found": 0, "match_type": "not_found",
                "dataset_variant_id": "", "dataset_ref": "",
                "dataset_alt_or_a1a2": "", "allele_status": "",
                "source": "", "source_kind": "",
            })
            continue
        for r in matches:
            ds_alleles = (
                f'{r["a1"]}/{r["a2"]}'
                if r["source_kind"] == "plink1_bim"
                else f'{r["ref"]}/{r["alt"]}'
            )
            out.append({
                "gene_symbol": a["gene_symbol"], "rsid": a["rsid"],
                "chrom_hg19": a["chrom_hg19"],
                "pos_hg19": a["pos_hg19"],
                "panel_ref": a["ref_allele"],
                "panel_effect_alt": a["effect_allele_alt"],
                "palindromic": a["palindromic"],
                "found": 1, "match_type": match_type,
                "dataset_variant_id": r["id"],
                "dataset_ref": r.get("ref", ""),
                "dataset_alt_or_a1a2": ds_alleles,
                "allele_status": allele_status(a, r),
                "source": r["source"], "source_kind": r["source_kind"],
            })
    return out


def role_for_column(col):
    low = col.strip().lower()
    roles = []
    for role, pats in ROLE_PATTERNS.items():
        if any(re.search(p, low) for p in pats):
            roles.append(role)
    return ";".join(roles)


def inspect_header(path: Path):
    try:
        with open_text(path) as fh:
            first = fh.readline()
    except Exception:
        return []
    if not first:
        return []
    if "\t" in first:
        delim = "\t"
    elif "," in first:
        delim = ","
    else:
        return []
    cols = next(csv.reader([first.rstrip("\r\n")], delimiter=delim))
    return [
        {
            "source": str(path), "column_index": i + 1,
            "column_name": c, "candidate_role": role_for_column(c),
        }
        for i, c in enumerate(cols)
    ]


def discover_phenotype_schema(root: Path, limit=500):
    rows = []
    n = 0
    for p in root.rglob("*"):
        if n >= limit:
            break
        if not p.is_file():
            continue
        low = p.name.lower()
        if not any(low.endswith(ext) for ext in PHENO_EXTS):
            continue
        if any(x in low for x in [".pvar", ".vcf", ".bim", ".fam", ".psam"]):
            continue
        hdr = inspect_header(p)
        if hdr and any(r["candidate_role"] for r in hdr):
            rows.extend(hdr)
            n += 1
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--koges-root", type=Path, required=True)
    ap.add_argument("--anchor-panel", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument(
        "--genome-build", choices=["unknown", "hg19", "GRCh37", "hg38", "GRCh38"],
        default="unknown"
    )
    args = ap.parse_args()

    if not args.koges_root.is_dir():
        raise SystemExit(f"KoGES root not found: {args.koges_root}")

    anchors = load_anchor_panel(args.anchor_panel)
    inventory = discover_genotype_metadata(args.koges_root)
    args.output_root.mkdir(parents=True, exist_ok=True)

    records = []
    source_status = []
    for item in inventory:
        p = Path(item["path"])
        try:
            if item["kind"] == "plink1_bim":
                recs = list(parse_bim(p))
                records.extend(recs)
                status = f"parsed_variants={len(recs)}"
            elif item["kind"] == "plink2_pvar":
                recs = list(parse_pvar(p))
                records.extend(recs)
                status = f"parsed_variants={len(recs)}"
            elif item["kind"] == "vcf":
                recs, status = parse_vcf_with_bcftools(p, anchors)
                records.extend(recs)
            else:
                status = "metadata_only_sidecar_required"
        except Exception as exc:
            status = "error:" + type(exc).__name__
        x = dict(item)
        x["audit_status"] = status
        source_status.append(x)

    build_confirmed = args.genome_build in {"hg19", "GRCh37"}
    coverage = match_anchors(records, anchors, build_confirmed)
    pheno = discover_phenotype_schema(args.koges_root)

    write_tsv(
        args.output_root / "STAGE4_KOGES_INPUT_INVENTORY.tsv",
        source_status,
        ["path", "kind", "size_bytes", "audit_status"],
    )
    write_tsv(
        args.output_root / "STAGE4_KOGES_ANCHOR_COVERAGE.tsv",
        coverage,
        [
            "gene_symbol", "rsid", "chrom_hg19", "pos_hg19",
            "panel_ref", "panel_effect_alt", "palindromic", "found",
            "match_type", "dataset_variant_id", "dataset_ref",
            "dataset_alt_or_a1a2", "allele_status", "source", "source_kind"
        ],
    )
    write_tsv(
        args.output_root / "STAGE4_KOGES_PHENOTYPE_SCHEMA.tsv",
        pheno,
        ["source", "column_index", "column_name", "candidate_role"],
    )

    found_genes = sorted({
        r["gene_symbol"] for r in coverage
        if int(r["found"]) == 1 and r["allele_status"] != "allele_mismatch"
    })
    all_genes = sorted({a["gene_symbol"] for a in anchors})
    summary = {
        "stage": "CKD Stage 4 KoGES input audit",
        "koges_root": str(args.koges_root),
        "genome_build_user_supplied": args.genome_build,
        "anchor_genes": all_genes,
        "anchor_genes_covered_allele_compatible": found_genes,
        "missing_or_incompatible_anchor_genes":
            sorted(set(all_genes) - set(found_genes)),
        "genotype_metadata_files": len(source_status),
        "candidate_phenotype_columns": sum(
            bool(r["candidate_role"]) for r in pheno
        ),
        "privacy": (
            "Only local file metadata, variant metadata, and phenotype column "
            "names are written. Individual-level genotype/phenotype values "
            "are not exported."
        ),
        "build_warning": (
            "GRCh37/hg19 chromosome-position fallback is only trusted when "
            "--genome-build hg19/GRCh37 is explicitly supplied. rsID matches "
            "remain usable independent of the coordinate-build fallback."
        ),
    }
    (args.output_root / "STAGE4_KOGES_AUDIT.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print("CKD_STAGE4_KOGES_AUDIT_PASS")


if __name__ == "__main__":
    main()
