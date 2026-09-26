#!/usr/bin/env python3
from pathlib import Path
import csv
import gzip
import zipfile
import math
import json
import sys

ROOT = Path("/srv/is-analysis")
INST = ROOT / "results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C2_FROZEN_INSTRUMENTS_ALL.tsv"
REG = ROOT / "results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv"
OUT = ROOT / "results/metabolic_resilience/stage3_confirmatory_mr/standardized_outcomes"
AUDIT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit"
OUT.mkdir(parents=True, exist_ok=True)

def split_line(s, tab):
    return s.rstrip("\n").split("\t") if tab else s.split()

def open_text(path, member=""):
    path = Path(path)
    if path.suffix.lower() == ".zip":
        z = zipfile.ZipFile(path)
        raw = z.open(member)
        import io
        return io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), z
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace"), None
    return open(path, "r", encoding="utf-8", errors="replace"), None

# Union instrument keys.
with INST.open("r", encoding="utf-8", newline="") as f:
    inst_rows = list(csv.DictReader(f, delimiter="\t"))

rsids = {r["rsid"] for r in inst_rows if r.get("rsid") not in ("","." ,"NA")}
positions = {(str(r["chrom_hg19"]).replace("chr",""), int(r["pos_hg19"])) for r in inst_rows}

registry = []
with REG.open("r", encoding="utf-8", newline="") as f:
    registry = list(csv.DictReader(f, delimiter="\t"))

summary = []
fail = 0

for spec in registry:
    trait = spec["trait"]
    if spec["status"] != "PASS":
        print(f"[SKIP] {trait}: registry status={spec['status']}")
        summary.append({"trait":trait,"matched":0,"status":"SKIP_REGISTRY_REVIEW"})
        continue

    path = spec["path"]
    member = spec.get("member","")
    fh, z = open_text(path, member)
    matched = []
    header = None
    H = None
    tab = True

    try:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("##"):
                continue
            if header is None:
                tab = "\t" in line
                header = split_line(line.lstrip("#"), tab)
                H = {x:i for i,x in enumerate(header)}
                continue

            p = split_line(line, tab)
            if len(p) < len(header):
                continue

            def val(col):
                return p[H[col]] if col and col in H and H[col] < len(p) else ""

            rsid = val(spec["col_rsid"])
            chrom = val(spec["col_chr"]).replace("chr","") if spec["col_chr"] else ""
            pos = None
            if spec["col_pos"]:
                try:
                    pos = int(float(val(spec["col_pos"])))
                except Exception:
                    pos = None

            is_match = (rsid in rsids) or ((chrom,pos) in positions if chrom and pos else False)
            if not is_match:
                continue

            ea = val(spec["col_ea"]).upper()
            oa = val(spec["col_oa"]).upper()

            beta = None
            if spec["col_beta"]:
                try:
                    beta = float(val(spec["col_beta"]))
                except Exception:
                    beta = None
            elif spec["col_or"]:
                try:
                    orv = float(val(spec["col_or"]))
                    beta = math.log(orv)
                except Exception:
                    beta = None

            try:
                se = float(val(spec["col_se"]))
            except Exception:
                se = None

            if beta is None or se is None or se <= 0:
                continue

            def fnum(col):
                if not col:
                    return None
                try:
                    return float(val(col))
                except Exception:
                    return None

            matched.append({
                "trait": trait,
                "domain": spec["domain"],
                "outcome_type": spec["type"],
                "source_build": spec["build"],
                "rsid": rsid,
                "chrom_source": chrom,
                "pos_source": pos if pos is not None else "",
                "effect_allele": ea,
                "other_allele": oa,
                "beta": beta,
                "se": se,
                "p": fnum(spec["col_p"]),
                "eaf": fnum(spec["col_eaf"]),
                "n": fnum(spec["col_n"]),
                "source_path": path,
            })
    finally:
        fh.close()
        if z is not None:
            z.close()

    # Dedupe exact variant representation.
    uniq = {}
    for r in matched:
        key = (r["rsid"],r["chrom_source"],r["pos_source"],r["effect_allele"],r["other_allele"])
        uniq[key] = r
    matched = list(uniq.values())

    dst = OUT / f"{trait}.instrument_outcomes.tsv.gz"
    fields = [
        "trait","domain","outcome_type","source_build","rsid","chrom_source","pos_source",
        "effect_allele","other_allele","beta","se","p","eaf","n","source_path"
    ]
    with gzip.open(dst, "wt", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(matched)

    status = "PASS" if matched else "FAIL_NO_MATCH"
    if status != "PASS":
        fail += 1
    summary.append({"trait":trait,"matched":len(matched),"status":status,"file":str(dst)})
    print(f"[{status}] {trait}: matched={len(matched)}")

sumfile = AUDIT / "STAGE3D1_OUTCOME_INSTRUMENT_EXTRACTION.tsv"
with sumfile.open("w", encoding="utf-8", newline="") as f:
    fields = list(summary[0].keys()) if summary else ["status"]
    wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    wr.writeheader()
    wr.writerows(summary)

sys.exit(0 if fail == 0 else 2)
