#!/usr/bin/env python3
"""Prepare EUR LD matrices and SuSiE inputs for CKD Stage 2C.

Inputs are the harmonized Stage 2B coloc files (GRCh37) plus Stage 2 candidate
anchors.  The script:
  * selects the EUR Berisa-Pickrell LD block containing each anchor (fallback
    +/-500 kb if the block file is unavailable or has no containing block),
  * downloads only the required chromosome-specific 1000 Genomes Phase 3
    GRCh37 PLINK2 files from the official PLINK resources page,
  * keeps EUR samples and removes known 1st-degree related samples,
  * matches reference variants by GRCh37 position + allele pair,
  * computes signed REF-based LD with PLINK2,
  * records the sign transform needed to orient LD to the pQTL effect allele.

The resulting LD matrix is intentionally kept outside the results directory so
large derived matrices do not need to be synchronized to Google Drive.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import html as html_lib
import json
import math
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
RELATED_URL = (
    "https://www.dropbox.com/s/0omyj2tyu7jmmw9/"
    "deg1_phase3.king.cutoff.out.id?dl=1"
)
VCF_BASE = "https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3"
LDETECT_URL = (
    "https://bitbucket.org/nygcresearch/ldetect-data/raw/master/"
    "EUR/fourier_ls-all.bed"
)

def phase3_vcf_url(chrom: str) -> str:
    if str(chrom).upper() == "X":
        name = "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
    else:
        name = (
            f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a."
            "20130502.genotypes.vcf.gz"
        )
    return f"{VCF_BASE}/{name}"


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace", newline="")


def read_tsv(path: Path):
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x):
    try:
        v = float(str(x).strip())
        return v if math.isfinite(v) else None
    except Exception:
        return None


def norm_chr(x):
    s = str(x).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    return s.upper()


def norm_a(x):
    return str(x).strip().upper()


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_tsv_gz(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def download(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 0:
        return
    part = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "curl", "-L", "--fail", "--retry", "5", "--retry-delay", "2",
        "-C", "-", "-o", str(part), url,
    ]
    subprocess.run(cmd, check=True)
    if not part.is_file() or part.stat().st_size == 0:
        raise RuntimeError(f"empty download: {url}")
    part.replace(dest)


def parse_related_ids(path: Path):
    out = set()
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals = re.split(r"\\s+", line)
            if vals:
                out.add(vals[-1])
    return out


def parse_ld_blocks(path: Path):
    blocks = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            vals = re.split(r"\s+", line.strip())
            if len(vals) < 3:
                continue
            chrom = norm_chr(vals[0])
            try:
                start0, end0 = int(float(vals[1])), int(float(vals[2]))
            except ValueError:
                continue
            if end0 > start0:
                blocks[chrom].append((start0 + 1, end0))
    for chrom in blocks:
        blocks[chrom].sort()
    return blocks


def choose_region(chrom, anchor_pos, input_positions, blocks, fallback_bp):
    block = next(
        ((s, e) for s, e in blocks.get(chrom, []) if s <= anchor_pos <= e),
        None,
    )
    in_lo, in_hi = min(input_positions), max(input_positions)
    if block:
        raw_lo, raw_hi = block
        lo, hi = max(in_lo, raw_lo), min(in_hi, raw_hi)
        source = "Berisa-Pickrell_EUR_hg19"
        truncated = int(lo != raw_lo or hi != raw_hi)
    else:
        raw_lo, raw_hi = max(1, anchor_pos - fallback_bp), anchor_pos + fallback_bp
        lo, hi = max(in_lo, raw_lo), min(in_hi, raw_hi)
        source = f"fallback_anchor_pm_{fallback_bp}"
        truncated = int(lo != raw_lo or hi != raw_hi)
    if lo >= hi:
        raise RuntimeError(f"{chrom}:{anchor_pos}: invalid selected region {lo}-{hi}")
    return lo, hi, source, raw_lo, raw_hi, truncated


def parse_pvar(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.rstrip("\n\r").lstrip("#").split("\t")
                continue
            if not line.strip():
                continue
            if header is None:
                raise RuntimeError(f"{path}: pvar header missing")
            vals = line.rstrip("\n\r").split("\t")
            r = dict(zip(header, vals))
            rows.append({
                "chrom": norm_chr(r.get("CHROM", "")),
                "pos": int(r["POS"]),
                "id": r["ID"],
                "ref": norm_a(r["REF"]),
                "alt": norm_a(r["ALT"]),
            })
    return rows


def match_reference(summary_rows, pvar_rows):
    by_pos = defaultdict(list)
    for r in pvar_rows:
        if "," in r["alt"]:
            continue
        by_pos[r["pos"]].append(r)
    matched = []
    missing = 0
    ambiguous = 0
    for s in summary_rows:
        pos = int(float(s["pos37"]))
        a0 = norm_a(s["allele0_pqtl"])
        ea = norm_a(s["allele1_pqtl"])
        hits = [
            r for r in by_pos.get(pos, [])
            if {r["ref"], r["alt"]} == {a0, ea}
        ]
        if len(hits) == 0:
            missing += 1
            continue
        if len(hits) > 1:
            # Prefer an exact rsID match if one exists.
            exact = [r for r in hits if r["id"] == s["snp"]]
            if len(exact) == 1:
                hits = exact
            else:
                ambiguous += 1
                continue
        r = hits[0]
        if ea == r["ref"]:
            sign = 1
        elif ea == r["alt"]:
            sign = -1
        else:
            raise RuntimeError("allele match invariant broken")
        matched.append({
            **s,
            "_ref_id": r["id"],
            "_ref_ref": r["ref"],
            "_ref_alt": r["alt"],
            "_ld_sign": sign,
        })
    return matched, missing, ambiguous


def discover_matrix_files(prefix: Path):
    candidates = list(prefix.parent.glob(prefix.name + "*.vcor1.bin"))
    if len(candidates) != 1:
        raise RuntimeError(f"{prefix}: expected one vcor1.bin, found {candidates}")
    matrix = candidates[0]
    vars_candidates = [
        matrix.with_suffix("").with_suffix(".vars"),
        Path(str(matrix).replace(".bin", ".vars")),
        prefix.parent / (prefix.name + ".unphased.vcor1.vars"),
    ]
    vars_path = next((p for p in vars_candidates if p.is_file()), None)
    if vars_path is None:
        extra = list(prefix.parent.glob(prefix.name + "*.vcor1.vars"))
        if len(extra) == 1:
            vars_path = extra[0]
    if vars_path is None:
        raise RuntimeError(f"{prefix}: matrix .vars companion not found")
    return matrix, vars_path


def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], check=True)


def count_psam(path: Path):
    n = 0
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--stage2b-input", type=Path, required=True)
    ap.add_argument("--work-root", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--plink2", default="plink2")
    ap.add_argument("--fallback-bp", type=int, default=500_000)
    ap.add_argument("--min-reference-maf", type=float, default=0.01)
    ap.add_argument("--min-ld-snps", type=int, default=100)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--memory-mb", type=int, default=2500)
    ap.add_argument("--cleanup-chromosome-cache", action="store_true")
    args = ap.parse_args()

    if shutil.which(args.plink2) is None:
        raise SystemExit("plink2 not found; on Ubuntu install with: sudo apt update && sudo apt install -y plink2")
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools not found; on Ubuntu install with: sudo apt update && sudo apt install -y bcftools")
    if shutil.which("curl") is None:
        raise SystemExit("curl is required")

    candidates = read_tsv(args.candidates)
    candidate_by_gene = {r["gene_symbol"].upper(): r for r in candidates}
    genes = sorted(candidate_by_gene)
    if not genes:
        raise SystemExit("no Stage 2 candidates")

    args.work_root.mkdir(parents=True, exist_ok=True)
    args.output_root.mkdir(parents=True, exist_ok=True)
    refroot = args.work_root / "reference_cache"
    ldroot = args.work_root / "ld"
    regionroot = args.work_root / "reference_regions"
    inputroot = args.output_root / "susie_input"
    for p in (refroot, ldroot, regionroot, inputroot):
        p.mkdir(parents=True, exist_ok=True)

    # Public reference metadata.
    panel = refroot / "integrated_call_samples_v3.20130502.ALL.panel"
    download(PANEL_URL, panel)
    eur_ids = []
    with panel.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if (r.get("super_pop") or "").upper() == "EUR":
                eur_ids.append(r["sample"])
    if len(eur_ids) < 400:
        raise RuntimeError(f"unexpectedly few 1KG EUR samples: {len(eur_ids)}")

    related = refroot / "deg1_phase3.king.cutoff.out.id"
    download(RELATED_URL, related)
    related_ids = parse_related_ids(related)
    eur_unrelated = [x for x in eur_ids if x not in related_ids]
    eur_keep = refroot / "EUR.unrelated.samples"
    eur_keep.write_text("\\n".join(eur_unrelated) + "\\n", encoding="utf-8")
    if len(eur_unrelated) < 450:
        raise RuntimeError(
            f"unexpectedly few unrelated 1KG EUR samples: {len(eur_unrelated)}"
        )

    ldetect = refroot / "EUR_fourier_ls-all.bed"
    blocks = {}
    try:
        download(LDETECT_URL, ldetect)
        blocks = parse_ld_blocks(ldetect)
    except Exception as exc:
        print(f"WARNING: EUR LD block download/parse failed ({exc}); using fallback anchor windows", file=sys.stderr)
        blocks = {}

    # Read and subset all Stage 2B harmonized inputs before chromosome downloads.
    per_gene = {}
    region_rows = []
    for gene in genes:
        path = args.stage2b_input / f"{gene}.tsv.gz"
        rows = read_tsv(path)
        if not rows:
            raise RuntimeError(f"{gene}: empty Stage2B coloc input")
        chroms = {norm_chr(r["chr37"]) for r in rows}
        if len(chroms) != 1:
            raise RuntimeError(f"{gene}: multiple chromosomes in input: {chroms}")
        chrom = next(iter(chroms))
        positions = [int(float(r["pos37"])) for r in rows]
        c = candidate_by_gene[gene]
        anchor_pos = int(float(c["anchor_pos_hg19"]))
        if norm_chr(c["anchor_chr_hg19"]) != chrom:
            raise RuntimeError(f"{gene}: candidate/input chromosome mismatch")
        lo, hi, source, raw_lo, raw_hi, truncated = choose_region(
            chrom, anchor_pos, positions, blocks, args.fallback_bp
        )
        subset = [r for r in rows if lo <= int(float(r["pos37"])) <= hi]
        if len(subset) < args.min_ld_snps:
            raise RuntimeError(f"{gene}: selected region has only {len(subset)} summary SNPs")
        per_gene[gene] = {"chrom": chrom, "lo": lo, "hi": hi, "rows": subset}
        fields = list(subset[0].keys())
        write_tsv_gz(inputroot / f"{gene}.tsv.gz", subset, fields)
        region_rows.append({
            "gene_symbol": gene,
            "anchor_rsid": c["anchor_rsid"],
            "chrom37": chrom,
            "anchor_pos37": anchor_pos,
            "region_start37": lo,
            "region_end37": hi,
            "region_source": source,
            "source_block_start37": raw_lo,
            "source_block_end37": raw_hi,
            "block_truncated_by_stage2b": truncated,
            "stage2b_input_snps": len(rows),
            "region_summary_snps": len(subset),
        })

    # Download/process only chromosomes actually required.
    genes_by_chr = defaultdict(list)
    for gene, info in per_gene.items():
        genes_by_chr[info["chrom"]].append(gene)

    qc = []
    provenance = {
        "panel_url": PANEL_URL,
        "related_ids_url": RELATED_URL,
        "vcf_base": VCF_BASE,
        "vcf_access": "bcftools indexed remote region extraction",
        "ld_block_url": LDETECT_URL,
        "ld_population": "1000 Genomes Phase 3 EUR",
        "ld_build": "GRCh37/hg19",
        "ld_statistic": "signed unphased dosage correlation, REF-based then sign-flipped to pQTL effect allele",
        "eur_samples_in_panel": len(eur_ids),
        "eur_samples_after_deg1_removal": len(eur_unrelated),
        "remove_related_resource": "deg1_phase3.king.cutoff.out.id",
        "min_reference_maf": args.min_reference_maf,
        "threads": args.threads,
        "memory_mb": args.memory_mb,
    }

    for chrom in sorted(genes_by_chr, key=lambda x: int(x) if x.isdigit() else 100):
        vcf_url = phase3_vcf_url(chrom)

        for gene in sorted(genes_by_chr[chrom]):
            info = per_gene[gene]
            regprefix = regionroot / gene
            region_vcf = regionroot / f"{gene}.1kg_eur.hg19.vcf.gz"

            # htslib performs indexed HTTP range requests against the public
            # chromosome VCF, so only this locus is transferred.
            if not region_vcf.is_file() or region_vcf.stat().st_size == 0:
                run([
                    "bcftools", "view",
                    "--regions", f"{chrom}:{info['lo']}-{info['hi']}",
                    "--samples-file", eur_keep,
                    "--min-alleles", "2",
                    "--max-alleles", "2",
                    "--types", "snps",
                    "--output-type", "z",
                    "--output-file", region_vcf,
                    vcf_url,
                ])
            run(["bcftools", "index", "--force", "--tbi", region_vcf])

            run([
                args.plink2,
                "--vcf", region_vcf,
                "--snps-only", "just-acgt",
                "--maf", args.min_reference_maf,
                "--set-all-var-ids", "@:#:$r:$a",
                "--make-pgen",
                "--threads", args.threads,
                "--memory", args.memory_mb, "require",
                "--out", regprefix,
            ])
            region_pvar = Path(str(regprefix) + ".pvar")
            region_psam = Path(str(regprefix) + ".psam")
            pvar_rows = parse_pvar(region_pvar)
            matched, missing, ambiguous = match_reference(info["rows"], pvar_rows)
            if len(matched) < args.min_ld_snps:
                raise RuntimeError(
                    f"{gene}: only {len(matched)} summary SNPs matched 1KG EUR reference"
                )

            # PLINK extraction is by reference variant ID.
            extract_ids = args.work_root / "extract" / f"{gene}.ids"
            extract_ids.parent.mkdir(parents=True, exist_ok=True)
            extract_ids.write_text(
                "\n".join(r["_ref_id"] for r in matched) + "\n",
                encoding="utf-8",
            )
            outprefix = ldroot / gene
            run([
                args.plink2,
                "--pfile", regprefix,
                "--extract", extract_ids,
                "--r-unphased", "square", "bin4", "ref-based",
                "--threads", args.threads,
                "--memory", args.memory_mb, "require",
                "--out", outprefix,
            ])
            matrix, vars_path = discover_matrix_files(outprefix)
            final_matrix = ldroot / f"{gene}.ld.bin"
            final_vars = ldroot / f"{gene}.ld.vars"
            if matrix != final_matrix:
                shutil.move(matrix, final_matrix)
            if vars_path != final_vars:
                shutil.move(vars_path, final_vars)

            var_order = [x.strip() for x in final_vars.read_text(encoding="utf-8").splitlines() if x.strip()]
            by_ref = {r["_ref_id"]: r for r in matched}
            if len(var_order) != len(set(var_order)):
                raise RuntimeError(f"{gene}: duplicate reference variant IDs in LD matrix")
            ordered = []
            for ref_id in var_order:
                if ref_id not in by_ref:
                    raise RuntimeError(f"{gene}: LD variant absent from mapping: {ref_id}")
                r = by_ref[ref_id]
                ordered.append({
                    "ref_id": ref_id,
                    "snp": r["snp"],
                    "pos37": r["pos37"],
                    "effect_allele": r["allele1_pqtl"],
                    "other_allele": r["allele0_pqtl"],
                    "reference_ref": r["_ref_ref"],
                    "reference_alt": r["_ref_alt"],
                    "effect_vs_ref_sign": r["_ld_sign"],
                })
            if len(ordered) < args.min_ld_snps:
                raise RuntimeError(f"{gene}: LD matrix only has {len(ordered)} variants")
            write_tsv(ldroot / f"{gene}.ld.meta.tsv", ordered, list(ordered[0]))

            expected_bytes = len(ordered) * len(ordered) * 4
            actual_bytes = final_matrix.stat().st_size
            if actual_bytes != expected_bytes:
                raise RuntimeError(
                    f"{gene}: LD binary size mismatch {actual_bytes} != {expected_bytes}"
                )
            qc.append({
                "gene_symbol": gene,
                "chrom37": chrom,
                "region_start37": info["lo"],
                "region_end37": info["hi"],
                "summary_snps_in_region": len(info["rows"]),
                "reference_region_variants": len(pvar_rows),
                "summary_reference_matched": len(matched),
                "summary_reference_missing": missing,
                "summary_reference_ambiguous": ambiguous,
                "ld_matrix_snps": len(ordered),
                "ld_reference_samples": count_psam(region_psam),
                "ld_matrix_bytes": actual_bytes,
                "ld_matrix": str(final_matrix),
                "ld_meta": str(ldroot / f"{gene}.ld.meta.tsv"),
            })

        if args.cleanup_chromosome_cache:
            # Regional VCFs are reproducible from the indexed UCSC/1000G mirror.
            for gene in genes_by_chr[chrom]:
                p = regionroot / f"{gene}.1kg_eur.hg19.vcf.gz"
                p.unlink(missing_ok=True)
                Path(str(p) + ".tbi").unlink(missing_ok=True)

    write_tsv(args.output_root / "STAGE2C_REGIONS.tsv", region_rows, list(region_rows[0]))
    write_tsv(args.output_root / "STAGE2C_LD_QC.tsv", qc, list(qc[0]))
    (args.output_root / "STAGE2C_LD_PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "genes": genes,
        "chromosomes": sorted(genes_by_chr),
        "eur_panel_n": len(eur_ids),
        "eur_unrelated_n": len(eur_unrelated),
        "ld_qc": qc,
    }, indent=2))
    print(f"CKD_STAGE2C_LD_PASS output={args.output_root} ld={ldroot}")


if __name__ == "__main__":
    main()
