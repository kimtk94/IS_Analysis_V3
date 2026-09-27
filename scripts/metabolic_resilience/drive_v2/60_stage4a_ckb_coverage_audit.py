#!/usr/bin/env python3
from pathlib import Path
import gzip
import csv
import json

ROOT=Path("/srv/is-analysis")
CKB=ROOT/"data/metabolic_resilience/stage1_pqtl/ckb/CKB_SomaScan_MR_coloc.gz"
OUT=ROOT/"results/metabolic_resilience/stage4_eas"
OUT.mkdir(parents=True,exist_ok=True)

TARGETS=[
    ("AOC1","P19801"),
    ("OGN","P20774"),
    ("TNFRSF6B","O95407"),
    ("TFPI","P10646"),
    ("SULT1A1","P50225"),
    ("CTRL","P40313"),
    ("IDUA","P35475"),
    ("COMT","P21964"),
]

with gzip.open(CKB,"rt",encoding="utf-8",errors="replace") as f:
    header=f.readline().rstrip("\n")
    cols=header.split("\t") if "\t" in header else header.split()

print("HEADER:")
for i,c in enumerate(cols,1):
    print(i,c)

hits={g:[] for g,_ in TARGETS}
with gzip.open(CKB,"rt",encoding="utf-8",errors="replace") as f:
    h=f.readline()
    tab="\t" in h
    for line in f:
        blob=line
        for gene,uni in TARGETS:
            if gene in blob or uni in blob:
                hits[gene].append(line.rstrip("\n"))

out=OUT/"STAGE4A_CKB_TARGET_RAW_HITS.tsv"
with out.open("w",encoding="utf-8") as f:
    f.write("gene\traw_line\n")
    for gene,_ in TARGETS:
        for line in hits[gene]:
            f.write(gene+"\t"+line.replace("\t"," | ")+"\n")

summary={
    "targets":len(TARGETS),
    "covered":sum(bool(hits[g]) for g,_ in TARGETS),
    "hits_per_gene":{g:len(hits[g]) for g,_ in TARGETS},
    "note":"CKB file is a precomputed MR/coloc result resource; this step only audits target coverage and does not reinterpret its schema.",
}
(OUT/"STAGE4A_CKB_COVERAGE.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
