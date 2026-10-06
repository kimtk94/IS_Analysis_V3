#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import shutil
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
BASE_SCRIPT = HERE / "prepare_ckd_stage2c_ld.py"

spec = importlib.util.spec_from_file_location(
    "ckd_stage2c_base",
    BASE_SCRIPT,
)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def read_tsv(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv_gz(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


def fnum(x):
    try:
        v = float(str(x).strip())
        return v if math.isfinite(v) else None
    except Exception:
        return None


def subset_window(rows, chrom, anchor_pos, window_bp):
    lo = max(1, anchor_pos - window_bp)
    hi = anchor_pos + window_bp

    out = [
        r for r in rows
        if (
            base.norm_chr(r["chr37"]) == base.norm_chr(chrom)
            and lo <= int(float(r["pos37"])) <= hi
        )
    ]

    return out, lo, hi


def unique_by_snp(rows):
    out = {}

    for r in rows:
        snp = r["snp"]

        if snp not in out:
            out[snp] = r
            continue

        old_p = fnum(out[snp].get("p_pqtl"))
        new_p = fnum(r.get("p_pqtl"))

        old_p = old_p if old_p is not None else 1.0
        new_p = new_p if new_p is not None else 1.0

        if new_p < old_p:
            out[snp] = r

    return out


def population_samples(panel_path: Path, related_path: Path):
    pops = {
        "EUR": [],
        "EAS": [],
    }

    with panel_path.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for r in reader:
            pop = (r.get("super_pop") or "").upper()

            if pop in pops:
                pops[pop].append(r["sample"])

    related = base.parse_related_ids(related_path)

    unrelated = {
        pop: [
            sample
            for sample in samples
            if sample not in related
        ]
        for pop, samples in pops.items()
    }

    return pops, unrelated



def phase3_vcf_urls(chrom):
    """
    Return equivalent indexed UCSC 1000G Phase3 mirrors.

    Both hosts expose the same GRCh37/hg19 Phase3 resource.
    """
    primary = base.phase3_vcf_url(chrom)
    filename = primary.rsplit("/", 1)[-1]

    mirrors = [
        primary,
        (
            "https://hgdownload.cse.ucsc.edu/"
            "gbdb/hg19/1000Genomes/phase3/"
            + filename
        ),
    ]

    # preserve order, remove duplicates
    return list(dict.fromkeys(mirrors))


def remote_vcf_samples_resilient(
    chrom,
    tries_per_url=3,
):
    errors = []

    for url in phase3_vcf_urls(chrom):

        for attempt in range(
            1,
            tries_per_url + 1,
        ):
            print(
                f"REMOTE_HEADER_ATTEMPT "
                f"chr{chrom} "
                f"{attempt}/{tries_per_url} "
                f"{url}",
                flush=True,
            )

            try:
                samples = base.remote_vcf_samples(
                    url
                )

                print(
                    f"REMOTE_HEADER_PASS "
                    f"chr{chrom} "
                    f"samples={len(samples)} "
                    f"url={url}",
                    flush=True,
                )

                return url, samples

            except Exception as exc:
                errors.append(
                    (
                        url,
                        attempt,
                        repr(exc),
                    )
                )

                print(
                    f"REMOTE_HEADER_RETRY "
                    f"chr{chrom} "
                    f"attempt={attempt} "
                    f"error={exc}",
                    flush=True,
                )

                time.sleep(
                    min(
                        20,
                        3 * attempt,
                    )
                )

    raise RuntimeError(
        f"chr{chrom}: all remote VCF "
        f"header attempts failed: {errors}"
    )


def extract_region_resilient(
    *,
    chrom,
    lo,
    hi,
    keep_file,
    region_vcf,
    preferred_url,
    tries_per_url=4,
):
    urls = [
        preferred_url,
        *phase3_vcf_urls(chrom),
    ]

    urls = list(
        dict.fromkeys(urls)
    )

    errors = []

    for url in urls:

        for attempt in range(
            1,
            tries_per_url + 1,
        ):
            region_vcf.unlink(
                missing_ok=True
            )

            Path(
                str(region_vcf) + ".tbi"
            ).unlink(
                missing_ok=True
            )

            print(
                f"REGION_DOWNLOAD_ATTEMPT "
                f"chr{chrom}:{lo}-{hi} "
                f"{attempt}/{tries_per_url} "
                f"url={url}",
                flush=True,
            )

            try:
                base.run([
                    "bcftools",
                    "view",

                    "--regions",
                    f"{chrom}:{lo}-{hi}",

                    "--samples-file",
                    keep_file,

                    "--min-alleles",
                    "2",

                    "--max-alleles",
                    "2",

                    "--types",
                    "snps",

                    "--output-type",
                    "z",

                    "--output-file",
                    region_vcf,

                    url,
                ])

                if (
                    region_vcf.is_file()
                    and region_vcf.stat().st_size > 0
                ):
                    print(
                        f"REGION_DOWNLOAD_PASS "
                        f"chr{chrom}:{lo}-{hi} "
                        f"bytes={region_vcf.stat().st_size} "
                        f"url={url}",
                        flush=True,
                    )

                    return url

                raise RuntimeError(
                    "bcftools returned success "
                    "but output is empty"
                )

            except Exception as exc:
                errors.append(
                    (
                        url,
                        attempt,
                        repr(exc),
                    )
                )

                print(
                    f"REGION_DOWNLOAD_RETRY "
                    f"chr{chrom}:{lo}-{hi} "
                    f"attempt={attempt} "
                    f"error={exc}",
                    flush=True,
                )

                region_vcf.unlink(
                    missing_ok=True
                )

                Path(
                    str(region_vcf) + ".tbi"
                ).unlink(
                    missing_ok=True
                )

                time.sleep(
                    min(
                        30,
                        5 * attempt,
                    )
                )

    raise RuntimeError(
        f"chr{chrom}:{lo}-{hi}: "
        f"all regional extraction attempts "
        f"failed: {errors}"
    )


def prepare_population_reference(
    *,
    pop,
    gene,
    chrom,
    lo,
    hi,
    summary_rows,
    samples,
    work_root,
    plink2,
    min_maf,
    threads,
    memory_mb,
):
    pop_root = work_root / pop

    region_root = pop_root / "reference_regions"
    pgen_root = pop_root / "pgen"
    freq_root = pop_root / "freq"
    sample_root = pop_root / "samples"

    for p in (
        region_root,
        pgen_root,
        freq_root,
        sample_root,
    ):
        p.mkdir(parents=True, exist_ok=True)

    vcf_url, vcf_sample_list = (
        remote_vcf_samples_resilient(
            chrom
        )
    )

    vcf_samples = set(
        vcf_sample_list
    )

    keep_samples = [
        x for x in samples
        if x in vcf_samples
    ]

    if len(keep_samples) < 400:
        raise RuntimeError(
            f"{gene}/{pop}: only {len(keep_samples)} reference samples"
        )

    keep_file = sample_root / f"chr{chrom}.{pop}.samples"

    base.write_sample_list(
        keep_file,
        keep_samples,
    )

    region_vcf = (
        region_root /
        f"{gene}.1kg_{pop.lower()}.hg19.vcf.gz"
    )

    if not base.cached_region_vcf_ok(
        region_vcf,
        keep_samples,
    ):
        region_vcf.unlink(missing_ok=True)
        Path(str(region_vcf) + ".tbi").unlink(missing_ok=True)

        vcf_url = extract_region_resilient(
            chrom=chrom,
            lo=lo,
            hi=hi,
            keep_file=keep_file,
            region_vcf=region_vcf,
            preferred_url=vcf_url,
        )

    if not base.cached_region_vcf_ok(
        region_vcf,
        keep_samples,
    ):
        raise RuntimeError(
            f"{gene}/{pop}: regional VCF validation failed"
        )

    base.run([
        "bcftools",
        "index",
        "--force",
        "--tbi",
        region_vcf,
    ])

    pgen_prefix = pgen_root / gene

    base.run([
        plink2,
        "--vcf",
        region_vcf,
        "--snps-only",
        "just-acgt",
        "--maf",
        min_maf,
        "--set-all-var-ids",
        "@:#:$r:$a",
        "--make-pgen",
        "--threads",
        threads,
        "--memory",
        memory_mb,
        "require",
        "--out",
        pgen_prefix,
    ])

    pvar = Path(str(pgen_prefix) + ".pvar")
    psam = Path(str(pgen_prefix) + ".psam")

    pvar_rows = base.parse_pvar(pvar)

    matched, missing, ambiguous = base.match_reference(
        summary_rows,
        pvar_rows,
    )

    freq_prefix = freq_root / gene

    base.run([
        plink2,
        "--pfile",
        pgen_prefix,
        "--freq",
        "--threads",
        threads,
        "--memory",
        memory_mb,
        "require",
        "--out",
        freq_prefix,
    ])

    afreq = base.parse_afreq(
        Path(str(freq_prefix) + ".afreq")
    )

    by_snp = {}
    ties = 0
    no_freq = 0

    for r in matched:
        fr = afreq.get(r["_ref_id"])

        if fr is None:
            no_freq += 1
            continue

        af = fr["alt_freq"]

        if abs(af - 0.5) < 1e-12:
            ties += 1
            continue

        major = fr["alt"] if af > 0.5 else fr["ref"]

        effect = base.norm_a(
            r["allele1_pqtl"]
        )

        if effect not in {
            fr["ref"],
            fr["alt"],
        }:
            continue

        x = dict(r)

        x["_reference_alt_freq"] = af
        x["_reference_maf"] = min(af, 1.0 - af)
        x["_reference_major"] = major
        x["_effect_vs_major_sign"] = (
            1 if effect == major else -1
        )

        if x["snp"] in by_snp:
            raise RuntimeError(
                f"{gene}/{pop}: duplicate SNP {x['snp']}"
            )

        by_snp[x["snp"]] = x

    return {
        "by_snp": by_snp,
        "pgen_prefix": pgen_prefix,
        "region_vcf": region_vcf,
        "reference_region_variants": len(pvar_rows),
        "matched_before_frequency": len(matched),
        "reference_missing": missing,
        "reference_ambiguous": ambiguous,
        "frequency_missing": no_freq,
        "major_ties_dropped": ties,
        "reference_samples": base.count_psam(psam),
    }


def create_ld(
    *,
    pop,
    gene,
    mapping,
    canonical_snps,
    pgen_prefix,
    work_root,
    output_root,
    plink2,
    threads,
    memory_mb,
):
    selected = [
        mapping[snp]
        for snp in canonical_snps
    ]

    extract_file = (
        work_root /
        pop /
        "extract" /
        f"{gene}.ids"
    )

    extract_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    extract_file.write_text(
        "\n".join(
            r["_ref_id"]
            for r in selected
        ) + "\n",
        encoding="utf-8",
    )

    raw_prefix = (
        work_root /
        pop /
        "ld_raw" /
        gene
    )

    raw_prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    base.run([
        plink2,
        "--pfile",
        pgen_prefix,
        "--extract",
        extract_file,
        "--r-unphased",
        "square",
        "bin4",
        "--threads",
        threads,
        "--memory",
        memory_mb,
        "require",
        "--out",
        raw_prefix,
    ])

    matrix, vars_path = base.discover_matrix_files(
        raw_prefix
    )

    ld_dir = output_root / "ld" / pop
    ld_dir.mkdir(parents=True, exist_ok=True)

    final_matrix = ld_dir / f"{gene}.ld.bin"
    final_vars = ld_dir / f"{gene}.ld.vars"
    final_meta = ld_dir / f"{gene}.ld.meta.tsv"

    final_matrix.unlink(missing_ok=True)
    final_vars.unlink(missing_ok=True)

    shutil.move(matrix, final_matrix)
    shutil.move(vars_path, final_vars)

    var_order = [
        x.strip()
        for x in final_vars.read_text(
            encoding="utf-8"
        ).splitlines()
        if x.strip()
    ]

    by_ref = {
        r["_ref_id"]: r
        for r in selected
    }

    ordered = []

    for ref_id in var_order:
        if ref_id not in by_ref:
            raise RuntimeError(
                f"{gene}/{pop}: LD ID missing from map: {ref_id}"
            )

        r = by_ref[ref_id]

        ordered.append({
            "ref_id": ref_id,
            "snp": r["snp"],
            "pos37": r["pos37"],
            "effect_allele": r["allele1_pqtl"],
            "other_allele": r["allele0_pqtl"],
            "reference_ref": r["_ref_ref"],
            "reference_alt": r["_ref_alt"],
            "reference_alt_freq": r["_reference_alt_freq"],
            "reference_maf": r["_reference_maf"],
            "reference_major_allele": r["_reference_major"],
            "effect_vs_ld_major_sign": r["_effect_vs_major_sign"],
        })

    if len(ordered) != len(canonical_snps):
        raise RuntimeError(
            f"{gene}/{pop}: LD SNP count mismatch "
            f"{len(ordered)} != {len(canonical_snps)}"
        )

    base.write_tsv(
        final_meta,
        ordered,
        list(ordered[0]),
    )

    expected_bytes = (
        len(ordered)
        * len(ordered)
        * 4
    )

    actual_bytes = final_matrix.stat().st_size

    if actual_bytes != expected_bytes:
        raise RuntimeError(
            f"{gene}/{pop}: LD binary size mismatch "
            f"{actual_bytes} != {expected_bytes}"
        )

    return {
        "ld_snps": len(ordered),
        "ld_matrix_bytes": actual_bytes,
    }


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--anchor", type=Path, required=True)
    ap.add_argument("--egfr-input", type=Path, required=True)
    ap.add_argument("--bun-input", type=Path, required=True)
    ap.add_argument("--work-root", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)

    ap.add_argument("--plink2", default="plink2")
    ap.add_argument("--window-bp", type=int, default=500_000)
    ap.add_argument("--min-reference-maf", type=float, default=0.01)
    ap.add_argument("--min-ld-snps", type=int, default=100)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--memory-mb", type=int, default=3500)

    args = ap.parse_args()

    for exe in (
        args.plink2,
        "bcftools",
        "curl",
    ):
        if shutil.which(exe) is None:
            raise SystemExit(
                f"{exe} not found"
            )

    args.work_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    refmeta = (
        args.work_root /
        "reference_metadata"
    )

    refmeta.mkdir(
        parents=True,
        exist_ok=True,
    )

    panel = (
        refmeta /
        "integrated_call_samples_v3.20130502.ALL.panel"
    )

    related = (
        refmeta /
        "deg1_phase3.king.cutoff.out.id"
    )

    base.download(
        base.PANEL_URL,
        panel,
    )

    base.download(
        base.RELATED_URL,
        related,
    )

    pop_all, unrelated = population_samples(
        panel,
        related,
    )

    print("===== 1000G POPULATION COUNTS =====")

    for pop in ("EUR", "EAS"):
        print(
            pop,
            "panel=",
            len(pop_all[pop]),
            "unrelated=",
            len(unrelated[pop]),
        )

        if len(unrelated[pop]) < 400:
            raise RuntimeError(
                f"{pop}: too few unrelated samples"
            )

    anchors = read_tsv(
        args.anchor
    )

    if len(anchors) != 9:
        raise SystemExit(
            f"expected 9 anchors; found {len(anchors)}"
        )

    qc = []

    for a in anchors:
        gene = (
            a["gene_symbol"]
            .strip()
            .upper()
        )

        chrom = base.norm_chr(
            a["chrom_hg19"]
        )

        anchor_pos = int(
            float(
                a["pos_hg19"]
            )
        )

        egfr_rows = read_tsv(
            args.egfr_input /
            f"{gene}.tsv.gz"
        )

        bun_rows = read_tsv(
            args.bun_input /
            f"{gene}.tsv.gz"
        )

        egfr_sub, lo, hi = subset_window(
            egfr_rows,
            chrom,
            anchor_pos,
            args.window_bp,
        )

        bun_sub, _, _ = subset_window(
            bun_rows,
            chrom,
            anchor_pos,
            args.window_bp,
        )

        egfr_by = unique_by_snp(
            egfr_sub
        )

        bun_by = unique_by_snp(
            bun_sub
        )

        if set(egfr_by) != set(bun_by):
            raise RuntimeError(
                f"{gene}: eGFR/BUN SNP sets differ"
            )

        summary_rows = list(
            egfr_by.values()
        )

        pop_obj = {}

        for pop in ("EUR", "EAS"):
            print()
            print(
                "============================================================"
            )
            print(
                gene,
                pop,
                f"{chrom}:{lo}-{hi}",
            )
            print(
                "============================================================"
            )

            pop_obj[pop] = prepare_population_reference(
                pop=pop,
                gene=gene,
                chrom=chrom,
                lo=lo,
                hi=hi,
                summary_rows=summary_rows,
                samples=unrelated[pop],
                work_root=args.work_root,
                plink2=args.plink2,
                min_maf=args.min_reference_maf,
                threads=args.threads,
                memory_mb=args.memory_mb,
            )

        common = (
            set(egfr_by)
            & set(bun_by)
            & set(pop_obj["EUR"]["by_snp"])
            & set(pop_obj["EAS"]["by_snp"])
        )

        canonical_snps = sorted(
            common,
            key=lambda snp: (
                int(
                    egfr_by[snp]["pos37"]
                ),
                snp,
            ),
        )

        if len(canonical_snps) < args.min_ld_snps:
            raise RuntimeError(
                f"{gene}: only {len(canonical_snps)} "
                "dual-reference SNPs"
            )

        input_root = (
            args.output_root /
            "susie_input"
        )

        for phenotype, source in (
            ("eGFR", egfr_by),
            ("BUN", bun_by),
        ):
            rows = [
                source[snp]
                for snp in canonical_snps
            ]

            write_tsv_gz(
                input_root /
                phenotype /
                f"{gene}.tsv.gz",
                rows,
                list(rows[0]),
            )

        ld = {}

        for pop in ("EUR", "EAS"):
            ld[pop] = create_ld(
                pop=pop,
                gene=gene,
                mapping=pop_obj[pop]["by_snp"],
                canonical_snps=canonical_snps,
                pgen_prefix=pop_obj[pop]["pgen_prefix"],
                work_root=args.work_root,
                output_root=args.output_root,
                plink2=args.plink2,
                threads=args.threads,
                memory_mb=args.memory_mb,
            )

        qc.append({
            "gene_symbol": gene,
            "chrom37": chrom,
            "anchor_pos37": anchor_pos,
            "window_start37": lo,
            "window_end37": hi,

            "summary_pm500k_snps": len(egfr_by),

            "eur_reference_matched":
                len(pop_obj["EUR"]["by_snp"]),

            "eas_reference_matched":
                len(pop_obj["EAS"]["by_snp"]),

            "dual_reference_common_snps":
                len(canonical_snps),

            "eur_reference_samples":
                pop_obj["EUR"]["reference_samples"],

            "eas_reference_samples":
                pop_obj["EAS"]["reference_samples"],

            "eur_ld_bytes":
                ld["EUR"]["ld_matrix_bytes"],

            "eas_ld_bytes":
                ld["EAS"]["ld_matrix_bytes"],
        })

        print(
            "DUAL_LD_GENE_PASS",
            gene,
            "common_snps=",
            len(canonical_snps),
        )

    base.write_tsv(
        args.output_root /
        "DUAL_LD_QC.tsv",
        qc,
        list(qc[0]),
    )

    provenance = {
        "stage":
            "CKD Stage4 cross-ancestry dual-LD",

        "build":
            "GRCh37/hg19",

        "window_bp":
            args.window_bp,

        "pQTL_LD":
            "1000 Genomes Phase 3 EUR",

        "outcome_LD":
            "1000 Genomes Phase 3 EAS",

        "canonical_snp_rule":
            (
                "Intersection of eGFR/BUN regional summary SNPs "
                "and allele-matched variants present in both "
                "1000G EUR and 1000G EAS."
            ),

        "qc":
            qc,
    }

    (
        args.output_root /
        "DUAL_LD_PROVENANCE.json"
    ).write_text(
        json.dumps(
            provenance,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print()
    print(
        "GENE\tSUMMARY\tEUR\tEAS\tCOMMON\tEUR_N\tEAS_N"
    )

    for r in qc:
        print(
            r["gene_symbol"],
            r["summary_pm500k_snps"],
            r["eur_reference_matched"],
            r["eas_reference_matched"],
            r["dual_reference_common_snps"],
            r["eur_reference_samples"],
            r["eas_reference_samples"],
            sep="\t",
        )

    print()
    print("CKD_STAGE4_DUAL_LD_PASS")


if __name__ == "__main__":
    main()
