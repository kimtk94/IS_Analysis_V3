#!/usr/bin/env python3
"""Exact C/T vs C/A nominal eQTL readback at rs7203740.
Post-readout source screen only. No regional coloc, no MR or causal gene claim.
"""
import argparse,collections,hashlib,importlib,json,math,time
from pathlib import Path
import pysam
meta=importlib.import_module("scripts.is.audit_official_p3a_qtl_catalogue_provenance")
POS=81775211
GENES={"ENSG00000197943":"PLCG2","ENSG00000153815":"CMIP","ENSG00000184860":"SDR42E1","ENSG00000086696":"HSD17B2"}
DATASETS={
"QTD000021":("QTS000002","BLUEPRINT","monocyte","191"),
"QTD000131":("QTS000015","GTEx","artery_aorta","387"),
"QTD000171":("QTS000015","GTEx","brain_cortex","205"),
"QTD000356":("QTS000015","GTEx","blood","670"),
"QTD000434":("QTS000025","ROSMAP","brain_naive","560")}
COLS="molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid".split()
def require(ok,msg):
    if not ok:raise ValueError(msg)
def parse(lines,header):
    require(set(COLS)<=set(header) and len(header)==len(set(header)),"Unexpected nominal schema")
    rows={};seen=0;alleles=collections.Counter();digest=hashlib.sha256()
    for line in lines:
        seen+=1;require(seen<=10000,"Single-position row bound")
        digest.update((line+"\n").encode())
        a=line.split("\t");require(len(a)==len(header),"QTL row width")
        d=dict(zip(header,a));require(d["chromosome"].removeprefix("chr")=="16" and int(d["position"])==POS,"Wrong locus")
        ref,alt=d["ref"],d["alt"];require(d["variant"]==f"chr16_{POS}_{ref}_{alt}","Variant key drift")
        alleles[ref+"/"+alt]+=1
        if d["gene_id"] not in GENES or (ref,alt) not in {("C","T"),("C","A")}:continue
        require(d["type"]=="SNP","Not SNP")
        require(d["rsid"]=="rs7203740","Exact rsID drift")
        beta,se,p,maf=map(float,[d["beta"],d["se"],d["pvalue"],d["maf"]])
        ac,an=map(int,[d["ac"],d["an"]])
        require(all(math.isfinite(x) for x in (beta,se,p,maf)) and se>0 and 0<=p<=1 and 0<maf<=.5 and 0<ac<an,"Invalid numeric row")
        key=(d["gene_id"],d["molecular_trait_id"],ref,alt)
        require(key not in rows,"Duplicate gene/allele/trait row")
        rows[key]={"gene":GENES[d["gene_id"]],"gene_id":d["gene_id"],"molecular_trait_id":d["molecular_trait_id"],"ref":ref,"alt":alt,"variant":d["variant"],"beta_ALT":beta,"se":se,"p":p,"maf":maf,"ac":ac,"an":an,"alt_frequency":ac/an,"genotyped_alleles_half_not_model_N":an/2}
    return {"position_rows":seen,"allele_pair_row_counts":dict(alleles),"position_stream_sha256":digest.hexdigest(),"candidate_rows":list(rows.values())}
def run():
    src=meta.get_official_sources();md=meta.table(src["meta"]);paths=meta.table(src["paths"])
    require("GRCh38 alternative allele (also the effect allele)" in src["columns"],"QTL orientation documentation missing")
    out={"status":"POST_READOUT_EXACT_ALLELE_QTL_SCREEN","GRCh38_position":"16:81775211","target":"C/T","other_allele":"C/A (separate descriptive control)","genes":GENES,"datasets":{},"limits":["Fixed five public nominal gene-count QTL sources chosen for blood/monocytes/vascular/brain context after outcome readout; not a genome-wide or independent confirmatory screen.","Nominal dataset N and AN/2 are not certified phenotype regression N or sdY. Population ancestry and overlap not established for matched EAS analysis.","Only a single position tested, no full cis universe, LD, colocalization, MR or causal claims. No row is not evidence of no QTL. Never substitute C/A for C/T."]}
    for i,(k,(study,label,group,n)) in enumerate(DATASETS.items()):
        if i:time.sleep(3)
        a,b=md[k],paths[k]
        for field,value in {"study_id":study,"study_label":label,"sample_group":group,"sample_size":n,"quant_method":"ge"}.items():require(a[field]==b[field]==value,"Pinned metadata mismatch "+k+" "+field)
        canonical=f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{study}/{k}/{k}.all.tsv.gz"
        require(b["ftp_path"]==canonical,"Unexpected source path")
        url=canonical.replace("ftp://","https://",1);d={"source":url,"study":label,"sample_group":group,"nominal_dataset_n":int(n),"pmid":a["pmid"]}
        try:
            with pysam.TabixFile(url) as t:
                ch="16" if "16" in t.contigs else "chr16"
                require(ch in t.contigs,"chr16 absent")
                headers=[x.lstrip("#").split("\t") for x in t.header]
                hh=[x for x in headers if set(COLS)<=set(x)]
                require(len(hh)<=1,"Ambiguous nominal header")
                d.update(parse(t.fetch(ch,POS-1,POS),hh[0] if hh else COLS),status="EXACT_POSITION_READ",header_origin="ORIGINAL" if hh else "PINNED_COLUMN_ORDER")
        except (OSError,ValueError,TypeError,KeyError,IndexError) as e:d.update(status="HOLD_SOURCE_OR_SCHEMA",reason=str(e)[:200])
        out["datasets"][k]=d
    out["exact_CT_candidate_rows"]=sum(v["alt"]=="T" for d in out["datasets"].values() for v in d.get("candidate_rows",[]))
    out["exact_CA_candidate_rows"]=sum(v["alt"]=="A" for d in out["datasets"].values() for v in d.get("candidate_rows",[]))
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run()
    except (OSError,ValueError,KeyError,TypeError) as e:x={"status":"HOLD_PINNED_SOURCE","reason":str(e)[:200]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
