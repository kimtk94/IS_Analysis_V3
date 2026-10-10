#!/usr/bin/env python3
"""Read-only Japanese-alcohol CS to European GIGASTROKE AIS transferability audit.

Input European "variant_pair_id" is an allele-pair identifier, NOT genome REF:ALT.
The GRCh37 FASTA independently determines REF, then outcome beta and EAF are
oriented to the alcohol GWAS' verified ALT. No MR, colocalization, or causal claim.
European ancestry results are cross-ancestry context, not EAS replication.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_EUR = "/srv/is-analysis/data/is/processed/gigastroke/broad_v1/GCST90104540_AIS_EUR_GRCh37.allele_pairs.tsv.gz"
DEFAULT_FASTA = "/srv/is-analysis/reference/grch37_1000g/human_g1k_v37.fasta"
DEFAULT_QC = DEFAULT_EUR + ".qc.json"

def finite_number(x):
    try:
        n = float(x)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None

def fasta_index(path):
    out = {}
    with open(str(path) + ".fai") as f:
        for line in f:
            v = line.rstrip("\n").split("\t")
            if len(v) < 5:
                raise ValueError("Malformed FASTA index")
            out[v[0]] = tuple(map(int, v[1:5]))
    return out

def fasta_base(fp, index, chrom, pos):
    name = chrom if chrom in index else "chr" + chrom
    if name not in index:
        raise ValueError("FASTA contig absent: " + chrom)
    length, offset, line_bases, line_width = index[name]
    if not (1 <= pos <= length):
        raise ValueError("FASTA coordinate out of bounds")
    k = pos - 1
    fp.seek(offset + (k // line_bases)*line_width + (k % line_bases))
    value = fp.read(1).decode("ascii").upper()
    if value not in "ACGT":
        raise ValueError("FASTA non-ACGT")
    return value

def harmonize_record(row, target):
    ea, oa = row["effect_allele"].upper(), row["other_allele"].upper()
    ref, alt = target["ref"].upper(), target["alt"].upper()
    if any(len(a) != 1 or a not in "ACGT" for a in (ea, oa, ref, alt)):
        raise ValueError("Non-SNV in target")
    if ea == oa or {ea, oa} != {ref, alt}:
        raise ValueError("GWAS effect/other allele set does not match REF/ALT")
    pair = row["variant_pair_id"].split(":")
    if len(pair) != 4 or pair[:2] != [target["chrom"], target["pos"]]:
        raise ValueError("GWAS pair variant coordinate conflict")
    if sorted(x.upper() for x in pair[2:]) != sorted((ea, oa)):
        raise ValueError("GWAS variant_pair_id allele set disagrees with EA/OA")
    pair_order = ("LEXICOGRAPHIC_ALLELE_PAIR" if pair[2:] == sorted([ea, oa]) else
                  "OTHER_EFFECT" if pair[2:] == [oa, ea] else
                  "EFFECT_OTHER" if pair[2:] == [ea, oa] else "UNEXPECTED")
    if pair_order == "UNEXPECTED":
        raise ValueError("GWAS variant_pair_id ambiguous allele order")
    if row["build"] != "GRCh37" or row["ancestry"] != "EUR" or row["phenotype"] != "AIS":
        raise ValueError("Unexpected GWAS context")
    beta = finite_number(row["beta"])
    se = finite_number(row["se"])
    pval = finite_number(row["p"])
    eaf = finite_number(row["eaf"])
    if beta is None or se is None or se <= 0 or pval is None or not 0 <= pval <= 1:
        raise ValueError("Invalid or missing GWAS association statistic")
    if eaf is None or not 0 <= eaf <= 1:
        raise ValueError("GWAS effect allele frequency invalid")
    sign = 1 if ea == alt else -1
    beta_alt = beta*sign
    eaf_alt = eaf if sign == 1 else 1-eaf
    maf = min(eaf_alt, 1-eaf_alt)
    status = ("EUR_ULTRA_RARE_MAF_LT_0_001" if maf < 0.001 else
              "EUR_RARE_MAF_LT_0_01" if maf < 0.01 else
              "EUR_REPORTED_COMMON_MAF_GE_0_01")
    return {
        "eur_beta_ALT": beta_alt, "eur_se": se, "eur_p": pval,
        "eur_eaf_ALT": eaf_alt, "eur_maf": maf, "eur_status": status,
        "eur_effect_allele": ea, "eur_other_allele": oa,
        "eur_reported_n": row.get("reported_n", ""),
        "eur_effect_allele_orientation": "EFFECT_ALT" if sign == 1 else "EFFECT_REF_FLIPPED",
        "eur_source_pair_id": row["variant_pair_id"],
        "eur_source_pair_id_order": pair_order,
    }

def scan_eur(path, targets):
    positions = {(str(v["chrom"]), int(v["pos"])) for v in targets}
    conditions = " || ".join(f'($5=="{chrom}" && $6=="{pos}")'
                            for chrom,pos in sorted(positions))
    expr = f"NR==1 || ({conditions})"
    gz = subprocess.Popen(["gzip", "-dc", str(path)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    filt = subprocess.Popen(["awk", "-F", "\t", expr],
                            stdin=gz.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    gz.stdout.close()
    data, err = filt.communicate()
    gz_error = gz.stderr.read()
    gz_return = gz.wait()
    gz.stderr.close()
    if gz_return != 0 or filt.returncode != 0:
        raise RuntimeError(f"gzip/awk failed: {gz_return},{filt.returncode} {gz_error[:150]!r} {err[:150]!r}")
    lines = data.decode("utf-8").splitlines()
    if not lines:
        raise ValueError("No GWAS header")
    columns = lines[0].rstrip("\r").split("\t")
    required = {"accession", "ancestry", "phenotype", "build", "chr", "pos",
                "variant_pair_id", "effect_allele", "other_allele", "beta",
                "se", "p", "eaf", "reported_n"}
    if not required.issubset(columns):
        raise ValueError(f"Unexpected EUR schema: missing {required-set(columns)}")
    hits = {}
    for line in lines[1:]:
        values = line.rstrip("\r").split("\t")
        if len(values) != len(columns):
            raise ValueError("Malformed GWAS record")
        row = dict(zip(columns, values))
        key = (row["chr"], int(row["pos"]))
        if key not in positions:
            continue
        hits.setdefault(key, []).append(row)
    return hits

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--audit-tsv", type=Path, required=True)
    p.add_argument("--alcohol-source", type=Path, required=True)
    p.add_argument("--eur", type=Path, default=Path(DEFAULT_EUR))
    p.add_argument("--eur-qc", type=Path, default=Path(DEFAULT_QC))
    p.add_argument("--fasta", type=Path, default=Path(DEFAULT_FASTA))
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    qc = json.loads(args.eur_qc.read_text())
    if (qc.get("accession") != "GCST90104540" or
        qc.get("ancestry") != "EUR" or
        qc.get("source_md5_verified") != "1bfe4ae8a5042fb24cb0a562b05b0d2f" or
        qc.get("phenotype") != "AIS"):
        raise ValueError("Unexpected EUR source QC metadata")
    with args.audit_tsv.open(newline="") as f:
        targets = list(csv.DictReader(f, delimiter="\t"))
    if len(targets) != 11 or len({x["variant_id"] for x in targets}) != 11:
        raise ValueError("Expected 11 unique upstream SNPs")
    if any(x["allele_mapping_status"] != "VERIFIED_ALLELE_SET" for x in targets):
        raise ValueError("Expected all 11 rsID / allele sets verified")
    variants = {}
    for gene in ("ADH1B", "ALDH2"):
        with (args.alcohol_source/gene/"variants.tsv").open(newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["ID"] in variants:
                    raise ValueError("Duplicate alcohol GWAS source SNP")
                variants[r["ID"]] = r
    index = fasta_index(args.fasta)
    with args.fasta.open("rb") as f:
        for t in targets:
            if fasta_base(f,index,str(t["chrom"]),int(t["pos"])) != t["ref"].upper():
                raise ValueError("GRCh37 FASTA REF mismatch: " + t["variant_id"])
    hits = scan_eur(args.eur, targets)
    results = []
    for t in targets:
        key = (t["chrom"], int(t["pos"]))
        e = variants.get(t["variant_id"])
        if e is None or e["reference_ALT"].upper() != t["alt"].upper():
            raise ValueError("Exposure ALT conflict or missing alcohol source")
        recs = hits.get(key, [])
        chosen = []
        for rec in recs:
            if {rec["effect_allele"].upper(),rec["other_allele"].upper()} == {t["ref"].upper(),t["alt"].upper()}:
                chosen.append(rec)
        if len(chosen) > 1:
            raise ValueError("Multiple matching EUR records; fail closed")
        base = {
            "locus":t["locus"], "rsid":t["allele_verified_rsids"],
            "variant_id_GRCh37_REF_ALT":t["variant_id"],
            "alcohol_beta_ALT":e["source_beta_ALT"],
            "alcohol_se":e["source_se"],
            "alcohol_p_recorded":e["source_original_p"],
            "Japanese_alcohol_ALT_EAF":e["Japanese_original_ALT_EAF"],
            "EUR_rows_same_position":len(recs),
            "EUR_rows_same_alleles":len(chosen),
        }
        if len(chosen) == 1:
            h=harmonize_record(chosen[0], t)
            pb=finite_number(e["source_beta_ALT"])
            h["exposure_outcome_sign"] = ("SAME" if pb*h["eur_beta_ALT"] > 0 else
                                          "OPPOSITE" if pb*h["eur_beta_ALT"] < 0 else "ZERO")
            base.update(h)
        elif recs:
            base["eur_status"] = "EUR_POSITION_PRESENT_ALLELE_PAIR_NOT_FOUND"
        else:
            base["eur_status"] = "EUR_VARIANT_NOT_REPORTED"
        results.append(base)
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError("Refusing to overwrite output")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cols=list(dict.fromkeys(col for row in results for col in row))
    dest=args.out_dir/"ALCOHOL_CS_GIGASTROKE_EUR_AIS_TRANSFER.tsv"
    with dest.open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(results)
    manifest = {
        "status":"CROSS_ANCESTRY_DESCRIPTIVE_ASSOCIATION_ONLY",
        "source":"GCST90104540_GIGASTROKE_EUR_AIS",
        "source_qc_md5_record":qc["source_md5_verified"],
        "source_n":qc["source_total_n"],
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "n_alcohol_snps":11,
        "n_eur_exact_allele_pair":sum(r["EUR_rows_same_alleles"]==1 for r in results),
        "n_eur_not_reported":sum(r["eur_status"]=="EUR_VARIANT_NOT_REPORTED" for r in results),
        "n_eur_low_maf_below_0_01":sum(r["eur_status"].startswith("EUR_RARE") or r["eur_status"].startswith("EUR_ULTRA_RARE") for r in results),
        "interpretation":"EUR ancestry transferability only; Japanese sample-overlap cannot be inferred from EUR ancestry alone",
        "population": "EUR_OUTCOME_VS_JAPANESE_ALCOHOL_EXPOSURE",
        "outcome_sample_overlap": "STUDY_COHORT_MEMBERSHIP_AND_PER_PARTICIPANT_OVERLAP_NOT_INDEPENDENTLY_AUDITED",
        "ALDH2_finemap":"BLOCKED",
        "ADH1B_finemap":"EXPLORATORY",
        "causal_inference":"NOT_PERFORMED"
    }
    (args.out_dir/"ALCOHOL_CS_GIGASTROKE_EUR_TRANSFER_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("EUR_TRANSFER_AUDIT_COMPLETE",json.dumps(manifest,sort_keys=True),flush=True)
    for r in results:
        print("RESULT",r["locus"],r["rsid"],r["eur_status"],
              "EUR_P",r.get("eur_p","NA"),"EUR_BETA_ALT",r.get("eur_beta_ALT","NA"),
              "EUR_MAF",r.get("eur_maf","NA"),flush=True)

if __name__ == "__main__":
    main()
