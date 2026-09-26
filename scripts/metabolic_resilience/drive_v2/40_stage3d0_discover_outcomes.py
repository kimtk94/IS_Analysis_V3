#!/usr/bin/env python3
from pathlib import Path
import gzip
import zipfile
import csv
import json
import re

ROOT = Path("/srv/is-analysis")
DATA = ROOT / "data/metabolic_resilience/stage2_gwas"
OUT = ROOT / "results/metabolic_resilience/stage3_confirmatory_mr"
AUDIT = ROOT / "results/metabolic_resilience/stage3_full_pgwas/audit"
OUT.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

TRAITS = [
    {
        "trait":"HDL","domain":"lipid","type":"quant","build":"GRCh37",
        "patterns":["without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"],
        "favorable":"POSITIVE",
    },
    {
        "trait":"TG","domain":"lipid","type":"quant","build":"GRCh37",
        "patterns":["without_UKB_logTG_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"SBP","domain":"bp","type":"quant","build":"RSID_PREFERRED",
        "patterns":["*GCST90310294*"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"DBP","domain":"bp","type":"quant","build":"RSID_PREFERRED",
        "patterns":["*GCST90310295*"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"BMI","domain":"adiposity","type":"quant","build":"GRCh37",
        "patterns":["*BMI*2018*.gz","*BMI*.sumstats*.gz","*BMI*.gz"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"WHR","domain":"adiposity","type":"quant","build":"GRCh37",
        "patterns":["*WHRadjBMI*2018*.gz","*WHRadjBMI*.gz","*WHR*.gz"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"FG","domain":"glycemia","type":"quant","build":"GRCh37",
        "patterns":["*GCST90002232*"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"HBA1C","domain":"glycemia","type":"quant","build":"GRCh37",
        "patterns":["*GCST90002244*"],
        "favorable":"NEGATIVE",
    },
    {
        "trait":"T2D","domain":"disease_validation","type":"cc","build":"GRCh37",
        "patterns":["Mahajan2018b_T2D_noUKBB_EUR.zip","*T2D*noUKBB*.zip","*T2D-noUKBB*"],
        "favorable":"NEGATIVE",
    },
]

ALIASES = {
    "chr": ["CHR","CHROM","chrom","chr","Chromosome","chromosome","#CHROM"],
    "pos": ["POS","BP","pos","position","Position","base_pair_location","GENPOS"],
    "rsid":["SNP","rsid","RSID","rsID","MarkerName","variant_id","ID"],
    "ea":  ["EA","effect_allele","Effect_allele","A1","ALLELE1","ALT","effectAllele"],
    "oa":  ["NEA","other_allele","Other_allele","A2","ALLELE0","REF","otherAllele"],
    "beta":["BETA","beta","Beta","Effect","effect","estimate","Estimate"],
    "or":  ["OR","or","OddsRatio","odds_ratio"],
    "se":  ["SE","se","StdErr","stderr","standard_error"],
    "p":   ["P","p","Pvalue","P_VALUE","p_value","P-value"],
    "eaf": ["EAF","eaf","A1FREQ","AF","effect_allele_frequency","Freq1"],
    "n":   ["N","n","N_total","TotalN","samplesize"],
}

def candidates_for(patterns):
    out = []
    for pat in patterns:
        out.extend(DATA.rglob(pat))
    uniq = {}
    for p in out:
        if p.is_file():
            uniq[str(p.resolve())] = p
    return list(uniq.values())

def read_header(path):
    member = ""
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as z:
            members = [n for n in z.namelist() if not n.endswith("/")]
            preferred = [n for n in members if "T2D-noUKBB" in n or n.endswith(".txt")]
            if not preferred:
                return "", "", "ZIP_NO_MEMBER"
            member = preferred[0]
            with z.open(member) as raw:
                for b in raw:
                    s = b.decode("utf-8", errors="replace").strip()
                    if s and not s.startswith("##"):
                        return s, member, "PASS"
        return "", member, "HEADER_NOT_FOUND"

    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "rt"
    with opener(path, mode, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("##"):
                return s, member, "PASS"
    return "", member, "HEADER_NOT_FOUND"

def split_header(s):
    if "\t" in s:
        return s.lstrip("#").split("\t")
    return s.lstrip("#").split()

def detect(cols):
    out = {}
    lower = {c.lower():c for c in cols}
    for key, aliases in ALIASES.items():
        hit = ""
        for a in aliases:
            if a in cols:
                hit = a
                break
            if a.lower() in lower:
                hit = lower[a.lower()]
                break
        out[key] = hit
    return out

rows = []
for spec in TRAITS:
    cand = candidates_for(spec["patterns"])
    # Prefer exact basename matches before broad wildcard hits.
    exact_names = {x for x in spec["patterns"] if "*" not in x and "?" not in x}
    cand.sort(key=lambda p: (0 if p.name in exact_names else 1, len(str(p)), str(p)))
    selected = cand[0] if cand else None

    if selected is None:
        rows.append({**spec, "path":"", "member":"", "header":"",
                     "status":"FAIL_NOT_FOUND", **{f"col_{k}":"" for k in ALIASES}})
        continue

    header, member, hstatus = read_header(selected)
    cols = split_header(header) if header else []
    mapping = detect(cols)

    required = ["se"]
    if not mapping["beta"] and not mapping["or"]:
        required.append("beta_or_or")
    # At least one matching key is needed.
    if not mapping["rsid"] and not (mapping["chr"] and mapping["pos"]):
        required.append("variant_key")

    status = "PASS" if hstatus == "PASS" and len(required) == 1 else "REVIEW"
    # If multiple candidate files from broad wildcard, force review unless selected exact basename.
    if len(cand) > 1 and selected.name not in exact_names:
        status = "REVIEW_MULTIPLE_CANDIDATES"

    row = {
        **spec,
        "path": str(selected),
        "member": member,
        "n_candidate_files": len(cand),
        "header": header,
        "status": status,
    }
    for k,v in mapping.items():
        row[f"col_{k}"] = v
    rows.append(row)

registry = OUT / "STAGE3D0_OUTCOME_REGISTRY.tsv"
fields = [
    "trait","domain","type","build","favorable","path","member","n_candidate_files",
    "col_chr","col_pos","col_rsid","col_ea","col_oa","col_beta","col_or",
    "col_se","col_p","col_eaf","col_n","header","status"
]
with registry.open("w", encoding="utf-8", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                        lineterminator="\n", extrasaction="ignore")
    wr.writeheader()
    wr.writerows(rows)

summary = {
    "traits": len(rows),
    "pass": sum(r["status"]=="PASS" for r in rows),
    "review": sum(r["status"].startswith("REVIEW") for r in rows),
    "fail": sum(r["status"].startswith("FAIL") for r in rows),
    "registry": str(registry),
    "rule": "Do not run confirmatory MR until every intended trait is PASS or manually reviewed.",
}
(AUDIT / "STAGE3D0_OUTCOME_REGISTRY.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(summary, indent=2))
for r in rows:
    print(r["trait"], r["status"], r["path"])
