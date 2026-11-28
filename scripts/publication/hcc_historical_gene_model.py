#!/usr/bin/env python3
"""Historical GRCh37 gene-model sensitivity; keep the current-ID result as primary."""
import argparse,collections,gzip,hashlib,io,json,math,re
from pathlib import Path
import hcc_tumor_biopsy_concordance as hcc
base=hcc.base;require=hcc.require;ROOT=hcc.ROOT
ARCH=ROOT/"docs/publication/evidence/hcc_exact_tumor_readback_20261125.json"
GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
CURRENT="ENSG00000278540"
HISTORICAL="ENSG00000132142"
def archive():
    d=json.loads(ARCH.read_text());q=d["aggregate"]
    require(d["run_id"]==45613227417 and d["job_id"]==139201877743 and q["matched_tumor_biopsies"]==22 and q["patients"]==9,"archived result identity")
    require(q["mapping"]["ACACA"]["ensembl"]==CURRENT and q["results"]["ACACA"]["RNA_rows"]==0 and q["results"]["ACACA"]["protein_rows"]==1,"historical sensitivity contract")
    return q
def annotation():
    base.HOSTS.add("ftp.ebi.ac.uk")
    raw,src=base.get(GTF,65_000_000)
    found=collections.defaultdict(list);n=0;headers=[]
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
        for b in f:
            line=b.decode("utf-8").rstrip("\r\n")
            if line.startswith("#"):
                if len(headers)<8:headers.append(line[:200])
                continue
            v=line.split("\t");require(len(v)==9,"GTF row width")
            if v[2]!="gene":continue
            n+=1
            attr=dict(re.findall(r'([A-Za-z_]+) "([^"]+)"',v[8]))
            g=attr.get("gene_name")
            if g not in hcc.orth.PANEL:continue
            found[g].append({"gene_id":attr.get("gene_id",""),"gene_type":attr.get("gene_type"),"chrom":v[0],"strand":v[6],"start":int(v[3]),"end":int(v[4])})
    require(n>=50000,"incomplete historical GTF")
    return {g:found[g] for g in hcc.orth.PANEL},{**src,"gene_rows":n,"header":headers}
def analyze():
    a=archive();maps,src=annotation()
    out={"status":"HCC_HISTORICAL_GENE_MODEL_SENSITIVITY","historical_annotation":src,"historical_fixed_panel":maps,
         "current_ID_ACACA_status":a["results"]["ACACA"]["status"]}
    old=maps["ACACA"]
    if not (len(old)==1 and old[0]["gene_id"].split(".")[0]==HISTORICAL and old[0]["gene_type"]=="protein_coding" and old[0]["chrom"].removeprefix("chr")=="17" and old[0]["strand"]=="-"):
        out["hold_reason"]="Historical unique ACACA gene model not established";return out
    # Do not relabel other genes or pool different gene models. Use only the documented
    # historical symbol identity in a separately labelled sensitivity.
    original=hcc.hgnc
    mapping=json.loads(json.dumps(a["mapping"]))
    mapping["ACACA"]={**mapping["ACACA"],"ensembl":HISTORICAL,"current_ensembl":CURRENT,
                         "sensitivity_basis":"GENCODE v19 GRCh37 exact unique gene_name ACACA; original RNA aligned to GRCh37; exact original annotation release not verified"}
    hcc.hgnc=lambda:mapping
    try:q=hcc.analyze()
    finally:hcc.hgnc=original
    require(q.get("matched_tumor_biopsies")==22 and q.get("patients")==9,"biopsy cohort drift")
    for g,v in a["results"].items():
        if g!="ACACA":require(q["results"][g]==v,"non-ACACA prior result drift: "+g)
    out.update(unchanged_other_genes=7,matched_tumor_biopsies=22,patients=9,historical_ACACA_result=q["results"]["ACACA"])
    out["limits"]=["Historical GRCh37 symbol-level gene annotation sensitivity, not proven exact original release or one-to-one gene-model equivalence. Keep current-ID primary missing result. Same nine patients, no p-value, independent validation, tumor-normal protein contrast or new mechanistic claim."]
    return out
def selftest():
    a=archive();require(len(a["results"])==8,"fixed gene panel")
    require(HISTORICAL!=CURRENT,"separate gene models")
    print("Prior cohort and missing current-ID gate passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["selftest","analyze"]);p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:z=analyze()
    except (OSError,EOFError,ValueError,KeyError,TypeError,IndexError,UnicodeError) as e:z={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:200]}
    z["scientific_status"]="HOLD_NOVELTY_AND_INDEPENDENT_VALIDATION"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print(json.dumps(z,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
