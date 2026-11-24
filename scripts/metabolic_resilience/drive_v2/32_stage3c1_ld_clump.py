#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import subprocess
import sys

ROOT = Path("/srv/is-analysis")
PRE = ROOT / "results/metabolic_resilience/stage3_full_pgwas/ld_clump"
REF = ROOT / "data/metabolic_resilience/stage2_gwas/ld_reference_1kg_eur/stage3_regional_pgen"
OUT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/ld_clump/stage3c1"
AUDIT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit"

PLINK2 = ROOT / "tools/plink2-a6.39-20260919/plink2"
if not PLINK2.exists():
    PLINK2 = Path("plink2")

GENES = ["AOC1","OGN","TNFRSF6B","TFPI","SULT1A1","CTRL","IDUA","COMT"]

MODES = {
    "primary": {
        "file_suffix": "primary_preclump.tsv",
        "log10p_min": 7.301029995663981,  # -log10(5e-8)
    },
    "stringent": {
        "file_suffix": "stringent_preclump.tsv",
        "log10p_min": 10.769551078621726,  # -log10(1.7e-11)
    },
}

CLUMP_R2 = 0.01
CLUMP_KB = 10000

OUT.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)


def load_pvar(path):
    by_pos = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            if header is None:
                header = line.lstrip("#").split()
                H = {x:i for i,x in enumerate(header)}
                required = {"CHROM","POS","ID","REF","ALT"}
                missing = required - set(H)
                if missing:
                    raise RuntimeError(f"{path}: missing pvar columns {sorted(missing)}")
                continue
            p = line.split()
            if len(p) < len(header):
                continue
            try:
                chrom = str(p[H["CHROM"]]).replace("chr","")
                pos = int(p[H["POS"]])
            except Exception:
                continue
            ref = p[H["REF"]].upper()
            alt = p[H["ALT"]].upper()
            if "," in alt:
                continue
            by_pos.setdefault((chrom,pos), []).append({
                "ld_id": p[H["ID"]],
                "ref": ref,
                "alt": alt,
            })
    return by_pos


def map_candidates(candidate_path, pvar_map):
    rows = []
    unmapped = []
    with candidate_path.open("r", encoding="utf-8", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            chrom = str(r["chrom_hg19"]).replace("chr","")
            pos = int(r["pos_hg19"])
            a0 = r["allele0"].upper()
            a1 = r["allele1"].upper()
            refs = pvar_map.get((chrom,pos), [])
            hits = [x for x in refs if {x["ref"],x["alt"]} == {a0,a1}]
            if len(hits) == 1:
                x = hits[0]
                orientation = "DIRECT" if (x["ref"] == a0 and x["alt"] == a1) else "SWAPPED"
                rr = dict(r)
                rr["ld_id"] = x["ld_id"]
                rr["ld_ref"] = x["ref"]
                rr["ld_alt"] = x["alt"]
                rr["ld_orientation"] = orientation
                rows.append(rr)
            else:
                rr = dict(r)
                rr["n_ld_hits"] = len(hits)
                unmapped.append(rr)
    return rows, unmapped


def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else ["status"]
    with path.open("w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                            lineterminator="\n", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)


def parse_clumps(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        if rd.fieldnames and "ID" in rd.fieldnames:
            return [r for r in rd if r.get("ID")]
    # PLINK output can be whitespace-delimited.
    lines = [x.strip() for x in path.read_text(errors="replace").splitlines() if x.strip()]
    if not lines:
        return []
    header = lines[0].split()
    H = {x:i for i,x in enumerate(header)}
    if "ID" not in H:
        raise RuntimeError(f"{path}: ID absent from .clumps")
    out = []
    for line in lines[1:]:
        p = line.split()
        if len(p) <= H["ID"]:
            continue
        out.append({"ID": p[H["ID"]]})
    return out


summary = []
overall_fail = 0

for gene in GENES:
    pvar = REF / f"{gene}.1000G_EUR.b37.pvar"
    pgen = REF / f"{gene}.1000G_EUR.b37.pgen"
    psam = REF / f"{gene}.1000G_EUR.b37.psam"
    ref_prefix = REF / f"{gene}.1000G_EUR.b37"

    if not (pvar.exists() and pgen.exists() and psam.exists()):
        print(f"[FAIL] {gene}: missing regional PGEN")
        overall_fail += 1
        continue

    pmap = load_pvar(pvar)

    for mode, spec in MODES.items():
        src = PRE / f"{gene}_{spec['file_suffix']}"
        if not src.exists():
            print(f"[FAIL] {gene}/{mode}: missing {src}")
            overall_fail += 1
            continue

        mapped, unmapped = map_candidates(src, pmap)

        mapped_path = OUT / f"{gene}_{mode}_mapped.tsv"
        unmapped_path = OUT / f"{gene}_{mode}_unmapped.tsv"
        write_tsv(mapped_path, mapped)
        write_tsv(unmapped_path, unmapped)

        clump_input = OUT / f"{gene}_{mode}.clump_input.tsv"
        with clump_input.open("w", encoding="utf-8", newline="") as f:
            wr = csv.writer(f, delimiter="\t", lineterminator="\n")
            wr.writerow(["ID","LOG10P"])
            for r in mapped:
                wr.writerow([r["ld_id"], r["log10p"]])

        out_prefix = OUT / f"{gene}_{mode}"
        cmd = [
            str(PLINK2),
            "--pfile", str(ref_prefix),
            "--clump", str(clump_input),
            "--clump-log10", "input-only",
            "--clump-id-field", "ID",
            "--clump-p-field", "LOG10P",
            "--clump-log10-p1", str(spec["log10p_min"]),
            "--clump-r2", str(CLUMP_R2),
            "--clump-kb", str(CLUMP_KB),
            "--out", str(out_prefix),
        ]

        print("\n[RUN]", " ".join(cmd))
        cp = subprocess.run(cmd, text=True)
        rc = cp.returncode

        clumps_path = Path(str(out_prefix) + ".clumps")
        clumps = parse_clumps(clumps_path) if rc == 0 else []
        index_ids = {r["ID"] for r in clumps}

        by_ld = {r["ld_id"]: r for r in mapped}
        independent = []
        for ld_id in index_ids:
            if ld_id in by_ld:
                rr = dict(by_ld[ld_id])
                rr["clump_mode"] = mode
                rr["clump_r2"] = CLUMP_R2
                rr["clump_kb"] = CLUMP_KB
                independent.append(rr)

        independent.sort(key=lambda r: float(r["log10p"]), reverse=True)
        independent_path = OUT / f"{gene}_{mode}_independent_instruments.tsv"
        write_tsv(independent_path, independent)

        missing_id_path = Path(str(out_prefix) + ".clumps.missing_id")
        missing_id_n = 0
        if missing_id_path.exists():
            missing_id_n = max(0, len(
                [x for x in missing_id_path.read_text(errors="replace").splitlines() if x.strip()]
            ) - 1)

        status = "PASS" if rc == 0 and len(independent) > 0 else "FAIL"
        if status != "PASS":
            overall_fail += 1

        row = {
            "gene": gene,
            "mode": mode,
            "preclump_n": sum(1 for _ in src.open("r", encoding="utf-8")) - 1,
            "mapped_n": len(mapped),
            "unmapped_n": len(unmapped),
            "mapped_pct": round(100*len(mapped)/(len(mapped)+len(unmapped)),3)
                if (mapped or unmapped) else 0,
            "plink_rc": rc,
            "missing_id_n": missing_id_n,
            "independent_n": len(independent),
            "clump_r2": CLUMP_R2,
            "clump_kb": CLUMP_KB,
            "status": status,
            "independent_file": str(independent_path),
        }
        summary.append(row)
        print("[RESULT]", row)

summary_path = AUDIT / "STAGE3C1_LD_CLUMP_SUMMARY.tsv"
write_tsv(summary_path, summary)

summary_json = {
    "targets": len(summary),
    "pass": sum(r["status"]=="PASS" for r in summary),
    "fail": sum(r["status"]!="PASS" for r in summary),
    "primary_total_independent": sum(r["independent_n"] for r in summary if r["mode"]=="primary"),
    "stringent_total_independent": sum(r["independent_n"] for r in summary if r["mode"]=="stringent"),
    "clump_r2": CLUMP_R2,
    "clump_kb": CLUMP_KB,
    "ld_reference": "1000 Genomes Phase 3 EUR503 GRCh37 v5b",
    "plink_ld": "phased/haplotype r2 default",
}
(AUDIT / "STAGE3C1_LD_CLUMP_SUMMARY.json").write_text(
    json.dumps(summary_json, indent=2) + "\n", encoding="utf-8"
)

print("\n" + json.dumps(summary_json, indent=2))
sys.exit(0 if overall_fail == 0 else 2)
