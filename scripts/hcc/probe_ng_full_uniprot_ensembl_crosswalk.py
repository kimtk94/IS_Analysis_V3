#!/usr/bin/env python3
"""Public Ng/Suter full-proteome UniProt -> Ensembl ID crosswalk audit.

Inputs: original source-pinned transient RNA/Protein feature ID labels.
Target: reviewed Homo sapiens UniProt cross-reference TSV (CURRENT annotation).
Only aggregate scientific QC is written; no sample IDs, accession lists,
gene identifiers, expression data or patient-level data are emitted.
"""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.parse
import urllib.request

COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
ENSG=re.compile(r"ENSG\d{11}")
UNIPROT=re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){2})")
LIMIT=35*1024*1024
URL="https://rest.uniprot.org/uniprotkb/stream?"+urllib.parse.urlencode({
  "query":"organism_id:9606 AND reviewed:true",
  "format":"tsv",
  "fields":"accession,xref_ensembl",
})

def read_ids(path):
    p=pathlib.Path(path)
    if not p.is_file() or p.stat().st_size>2_000_000:
        raise ValueError("transient source ID file missing or oversized")
    t,pids=set(),[]
    for line in p.read_text(encoding="utf-8").splitlines():
        kind,sep,value=line.partition("\t")
        if not sep or kind not in ("T","P") or not value:
            raise ValueError("invalid source feature ID row")
        if kind=="T":
            if not ENSG.fullmatch(value):raise ValueError("non-Ensembl transcript ID")
            t.add(value)
        else:pids.append(value)
    if len(t)!=20356 or len(pids)!=5585 or len(set(pids))!=5585:
        raise ValueError("source feature label cardinality changed")
    return t,pids

def canonical_protein(label):
    x=label.split("|")[1] if label.startswith(("sp|","tr|")) and "|" in label[3:] else label
    x=x.split("-")[0]
    return x if UNIPROT.fullmatch(x) else None

def parse_uniprot_tsv(data):
    raw=io.StringIO(data.decode("utf-8-sig"))
    reader=csv.DictReader(raw,delimiter="\t")
    names=reader.fieldnames or []
    if "Entry" not in names or not any(x.lower().startswith("ensembl") for x in names):
        raise ValueError("unexpected UniProt crossref TSV headers")
    ec=next(x for x in names if x.lower().startswith("ensembl"))
    out={}
    rows=0
    for item in reader:
        rows+=1
        key=item.get("Entry","").strip()
        if not UNIPROT.fullmatch(key):continue
        genes=set(ENSG.findall(item.get(ec,"") or ""))
        if key in out:out[key].update(genes)
        else:out[key]=genes
    return out,rows

def summarize(t,pids,refs,rows,checksum):
    aliases={p:canonical_protein(p) for p in pids}
    valid=[a for a in aliases.values() if a]
    matched={p:refs.get(a,set()) & t for p,a in aliases.items() if a}
    mapped={p:gs for p,gs in matched.items() if gs}
    unique={p:next(iter(gs)) for p,gs in mapped.items() if len(gs)==1}
    gene_protein_count={}
    for g in unique.values():gene_protein_count[g]=gene_protein_count.get(g,0)+1
    unambiguous_genes=sum(v==1 for v in gene_protein_count.values())
    comparable=len(set(g for gs in mapped.values() for g in gs))
    return {
        "schema":"hcc-ng-full-uniprot-ensembl-live-audit-v1",
        "source_git_revision":COMMIT,
        "annotation_source":"UniProtKB reviewed Homo sapiens current crossrefs",
        "annotation_snapshot_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "annotation_tsv_sha256":checksum,
        "annotation_records_read":rows,
        "source_RNA_ENSG_n":len(t),
        "source_proteome_features_n":len(pids),
        "source_protein_features_with_recognized_UniProt_accession_n":len(valid),
        "protein_features_with_Ensembl_gene_in_source_RNA_n":len(mapped),
        "protein_features_with_unique_Ensembl_gene_n":len(unique),
        "distinct_Ensembl_genes_corresponding_to_proteins_n":comparable,
        "one_gene_one_protein_unambiguous_genes_n":unambiguous_genes,
        "multi_gene_protein_feature_count":sum(len(gs)>1 for gs in mapped.values()),
        "unmapped_or_noncanonical_source_protein_n":len(pids)-len(mapped),
        "source_paper_comparable_gene_benchmark_n":5490,
        "source_paper_correlated_gene_benchmark_n":5481,
        "absolute_deviation_from_paper_comparable_n":abs(comparable-5490),
        "exact_original_paper_2022_mapping_reproduced":False,
        "clinical_effect_or_pathway_significance_tested":False,
        "patient_sample_gene_accession_or_expression_values_exported":False,
        "interpretation":"Current UniProt crossrefs vs GRCh37-era published Ng gene universe; replicate original mappings only with versioned source annotations and isoform policy."
    }

def request_mapping():
    req=urllib.request.Request(URL,headers={"Accept":"text/tab-separated-values",
                                            "User-Agent":"HCC-public-source-mapping-QC/1.0"})
    with urllib.request.urlopen(req,timeout=105) as r:
        chunks=[];size=0
        while True:
            chunk=r.read(128*1024)
            if not chunk:break
            size+=len(chunk)
            if size>LIMIT:raise ValueError("annotation download exceeds bound")
            chunks.append(chunk)
    data=b"".join(chunks)
    if len(data)<2000:raise ValueError("unexpected tiny annotation result")
    return data

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-features",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    t,proteins=read_ids(args.private_features)
    data=request_mapping()
    refs,rows=parse_uniprot_tsv(data)
    result=summarize(t,proteins,refs,rows,hashlib.sha256(data).hexdigest())
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    # A live-access/data schema error cannot silently become a biological PASS.
    if result["protein_features_with_Ensembl_gene_in_source_RNA_n"]==0:
        raise SystemExit("no source protein mapped to source RNA")
if __name__=="__main__":main()
