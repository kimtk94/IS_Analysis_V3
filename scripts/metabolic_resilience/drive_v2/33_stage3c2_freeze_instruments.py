#!/usr/bin/env python3
from pathlib import Path
import csv
import gzip
import json
import sys

ROOT = Path("/srv/is-analysis")
CLUMP = ROOT / "results/metabolic_resilience/stage3_full_pgwas/ld_clump/stage3c1"
VCF = ROOT / "data/metabolic_resilience/stage2_gwas/ld_reference_1kg_eur/stage3_regional_vcf"
OUT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/instruments"
AUDIT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit"

GENES = ["AOC1","OGN","TNFRSF6B","TFPI","SULT1A1","CTRL","IDUA","COMT"]
MODES = ["primary","stringent"]

OUT.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

COMP = str.maketrans("ACGT", "TGCA")
def complement(a):
    return a.translate(COMP)

def is_pal(a,b):
    a,b = a.upper(),b.upper()
    return len(a)==1 and len(b)==1 and {a,b} in ({"A","T"},{"C","G"})

def load_1000g_vcf(path):
    by_key = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            chrom = p[0].replace("chr","")
            try:
                pos = int(p[1])
            except Exception:
                continue
            rsid = p[2]
            ref = p[3].upper()
            for alt in p[4].upper().split(","):
                by_key[(chrom,pos,ref,alt)] = rsid
                by_key[(chrom,pos,alt,ref)] = rsid
    return by_key

all_rows = []
summary = []
fail = 0

for gene in GENES:
    vcf_path = VCF / f"{gene}.1000G_EUR.b37.vcf.gz"
    if not vcf_path.exists():
        print(f"[FAIL] missing VCF {vcf_path}")
        fail += 1
        continue

    rsmap = load_1000g_vcf(vcf_path)

    for mode in MODES:
        src = CLUMP / f"{gene}_{mode}_independent_instruments.tsv"
        if not src.exists():
            print(f"[FAIL] missing {src}")
            fail += 1
            continue

        rows = []
        with src.open("r", encoding="utf-8", newline="") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for r in rd:
                chrom = str(r["chrom_hg19"]).replace("chr","")
                pos = int(r["pos_hg19"])
                a0 = r["allele0"].upper()
                a1 = r["allele1"].upper()
                ld_ref = r.get("ld_ref","").upper()
                ld_alt = r.get("ld_alt","").upper()

                rsid = rsmap.get((chrom,pos,ld_ref,ld_alt), ".")
                rr = dict(r)
                rr["mode"] = mode
                rr["rsid"] = rsid
                rr["palindromic"] = int(is_pal(a0,a1))
                rr["exposure_effect_allele"] = r["effect_allele"].upper()
                rr["exposure_other_allele"] = r["other_allele"].upper()
                rr["exposure_eaf"] = r["eaf"]
                rr["exposure_beta"] = r["beta"]
                rr["exposure_se"] = r["se"]
                rr["exposure_F"] = r["f_stat"]
                rows.append(rr)

        k = len(rows)
        if k == 1:
            plan = "WALD_ONLY"
        elif 2 <= k < 3:
            plan = "IVW"
        elif 3 <= k < 10:
            plan = "IVW_WEIGHTED_MEDIAN"
        else:
            plan = "IVW_WEIGHTED_MEDIAN_EGGER"

        for r in rows:
            r["mr_plan"] = plan
            r["k_instruments"] = k

        dst = OUT / f"{gene}_{mode}_frozen_instruments.tsv"
        if rows:
            fields = list(rows[0].keys())
        else:
            fields = ["gene_symbol","mode","status"]

        with dst.open("w", encoding="utf-8", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
            wr.writeheader()
            wr.writerows(rows)

        all_rows.extend(rows)
        summary.append({
            "gene": gene,
            "mode": mode,
            "k": k,
            "rsid_available": sum(r["rsid"] not in ("",".") for r in rows),
            "palindromic": sum(int(r["palindromic"]) for r in rows),
            "min_F": min([float(r["exposure_F"]) for r in rows], default=None),
            "max_F": max([float(r["exposure_F"]) for r in rows], default=None),
            "mr_plan": plan,
            "status": "PASS" if k > 0 else "FAIL_NO_IV",
            "file": str(dst),
        })
        if k == 0:
            fail += 1

combined = OUT / "STAGE3C2_FROZEN_INSTRUMENTS_ALL.tsv"
if all_rows:
    fields = list(all_rows[0].keys())
    with combined.open("w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                            lineterminator="\n", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(all_rows)

sumfile = AUDIT / "STAGE3C2_INSTRUMENT_FREEZE.tsv"
with sumfile.open("w", encoding="utf-8", newline="") as f:
    fields = list(summary[0].keys()) if summary else ["status"]
    wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    wr.writeheader()
    wr.writerows(summary)

j = {
    "sets": len(summary),
    "pass": sum(r["status"]=="PASS" for r in summary),
    "fail": sum(r["status"]!="PASS" for r in summary),
    "primary_k_total": sum(r["k"] for r in summary if r["mode"]=="primary"),
    "stringent_k_total": sum(r["k"] for r in summary if r["mode"]=="stringent"),
    "egger_rule": "MR-Egger only when k>=10",
    "weighted_median_rule": "weighted median when k>=3",
    "wald_rule": "Wald ratio when k=1",
}
(AUDIT / "STAGE3C2_INSTRUMENT_FREEZE.json").write_text(
    json.dumps(j, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(j, indent=2))
sys.exit(0 if fail == 0 else 2)
