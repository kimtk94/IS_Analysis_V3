#!/usr/bin/env python3
"""Complete remaining exact rs7203740 tissue screens and screen rs2767060.
Post-readout descriptive nominal QTL, bounded index/header/range retrieval.
"""
import argparse,collections,hashlib,importlib,json,math,time
from pathlib import Path
base=importlib.import_module("scripts.is.rs7203740_qtl_http_recovery")
http=base.http
prev=base.prev
meta=base.meta
COLS=base.COLS
POS10=35833363
GENES10={
"ENSG00000294862":"lncRNA ENSG00000294862","ENSG00000280719":"PCAT5","ENSG00000309863":"lncRNA ENSG00000309863",
"ENSG00000227313":"LINC02630","ENSG00000301879":"lncRNA ENSG00000301879","ENSG00000299958":"lncRNA ENSG00000299958",
"ENSG00000299999":"lncRNA ENSG00000299999","ENSG00000226646":"RPL7P37","ENSG00000287528":"lncRNA ENSG00000287528",
"ENSG00000288893":"lncRNA ENSG00000288893","ENSG00000177283":"FZD8","ENSG00000284265":"MIR4683",
"ENSG00000177291":"GJD4","ENSG00000273312":"GJD4-AS1","ENSG00000227475":"lncRNA ENSG00000227475"}
PLAN={"QTD000131":["rs7203740","rs2767060"],"QTD000171":["rs7203740"],"QTD000434":["rs7203740"],"QTD000356":["rs2767060"]}
def require(x,m):
    if not x:raise ValueError(m)
def parse10(lines,header):
    require(set(COLS)<=set(header) and len(header)==len(set(header)),"Nominal schema")
    out=[];seen=set();alleles=collections.Counter();gids=set();digest=hashlib.sha256()
    for line in lines:
        digest.update((line+"\n").encode())
        a=line.split("\t");require(len(a)==len(header),"Row width");d=dict(zip(header,a))
        require(d["chromosome"].removeprefix("chr")=="10" and int(d["position"])==POS10,"Position")
        ref,alt=d["ref"],d["alt"];require(d["variant"]==f"chr10_{POS10}_{ref}_{alt}","Variant key")
        alleles[ref+"/"+alt]+=1
        if (ref,alt)!=("C","T"):continue
        require(d["type"]=="SNP" and d["rsid"]=="rs2767060","Exact SNP and rsID")
        gid=d["gene_id"];gids.add(gid)
        if gid not in GENES10:continue
        key=(gid,d["molecular_trait_id"],ref,alt);require(key not in seen,"Duplicate trait");seen.add(key)
        beta,se,p,maf=map(float,[d["beta"],d["se"],d["pvalue"],d["maf"]])
        ac,an=map(int,[d["ac"],d["an"]])
        require(all(math.isfinite(v) for v in (beta,se,p,maf)) and se>0 and 0<=p<=1 and 0<maf<=.5 and 0<ac<an,"Invalid QTL effect")
        out.append({"gene":GENES10[gid],"gene_id":gid,"molecular_trait_id":d["molecular_trait_id"],"ref":ref,"alt":alt,"beta_ALT":beta,"se":se,"p":p,"maf":maf,"ac":ac,"an":an,"alt_frequency":ac/an,"genotyped_alleles_half_not_model_N":an/2})
    return {"position_rows":len(lines),"allele_pair_row_counts":dict(alleles),"position_stream_sha256":digest.hexdigest(),"all_CT_gene_ids":sorted(gids),"candidate_rows":out}
def load_index(url):
    idx,_=http.get(url+".tbi",8000000);t=http.parse_index(idx)
    require(t["seq"]==2 and t["begin"]==3 and t["end"] in (0,3),"Index coordinate schema")
    time.sleep(4)
    first,cr=http.get(url,65536,(0,65535));total=int(cr.rsplit("/",1)[1])
    firstrow=base.parse_bgzf_range(first,0,0)["data_lines"].readline().lstrip("#").rstrip("\r\n").split("\t")
    if set(COLS)<=set(firstrow):h=firstrow;origin="ORIGINAL_HEADER"
    else:
        require(len(firstrow)==len(COLS) and firstrow[1].removeprefix("chr").isdigit() and firstrow[2].isdigit() and firstrow[5]==f"chr{firstrow[1].removeprefix('chr')}_{firstrow[2]}_{firstrow[3]}_{firstrow[4]}","Pinned column order not corroborated")
        h=COLS;origin="PINNED_COLUMN_ORDER_WITH_FIRST_ROW_CHECK"
    require(len(set(h))==len(h) and h[1:3]==["chromosome","position"],"Nominal header")
    return t,h,total,{"header_origin":origin,"index_sha256":http.sha(idx),"header_range_sha256":http.sha(first),"source_total_bytes":total,"original_full_file_downloaded":False}
def read_exact(url,t,h,total,rs):
    ch,pos,parse=("16",prev.POS,prev.parse) if rs=="rs7203740" else ("10",POS10,parse10)
    span=http.indexed_span(t,ch,pos,pos,total)
    require(span is not None,"No indexed chunk at position")
    a,b,x,y=span;require(b-a+1<=4000000,"Single-position compressed budget")
    time.sleep(4)
    raw,cr=http.get(url,b-a+1,(a,b));stream=base.parse_bgzf_range(raw,x,y)["data_lines"]
    rows=[];n=0
    for line in stream:
        n+=1;require(n<=500000,"Decoded row cap")
        s=line.rstrip("\r\n");v=s.split("\t");require(len(v)==len(h),"Indexed row width")
        if v[1].removeprefix("chr")==ch and v[2]==str(pos):rows.append(s)
    require(len(rows)<=10000,"Exact row cap")
    result=parse(rows,h)
    return {**result,"status":"EXACT_POSITION_READ_HTTP_RANGE","position":f"{ch}:{pos}","rsid":rs,"source_range":cr,"source_range_sha256":http.sha(raw),"decoded_rows":n}
def run(root):
    # Prior anchored positions and gene identity must be fixed before further effects.
    for path,rs,pos,alleles in [("is_rs2767060_reference_locus_20261125.json","rs2767060",POS10,"C/T"),("is_rs7203740_reference_locus_20261125.json","rs7203740",prev.POS,"C/A/T")]:
        x=json.loads((root/path).read_text())["result"]
        require(x["identity_verified"] and x["rsid"]==rs and x["GRCh38_pos"]==pos and x["ensembl_allele_string"]==alleles and x["historical_alt"]=="T" and x["historical_ref"]=="C","Prior verified anchor identity")
        if rs=="rs2767060":require(set(GENES10)=={g["id"] for g in x["regional_genes"]},"Frozen positional gene set")
    prior=json.loads((root/"is_rs7203740_qtl_http_recovery_20261125.json").read_text())["result"]
    require(prior["sources_read"]==2 and set(prior["datasets"])=={"QTD000021","QTD000356"},"Prior monocyte/blood source recovery")
    src=meta.get_official_sources();md=meta.table(src["meta"]);paths=meta.table(src["paths"])
    require("GRCh38 alternative allele (also the effect allele)" in src["columns"],"Orientation contract")
    out={"status":"POST_READOUT_EXACT_TWO_LOCUS_QTL_EXTENSION","plan":PLAN,"chr10_genes":GENES10,"chr16_genes":prev.GENES,"prior_rs7203740_read_sources":2,"datasets":{},"limits":["Prior 2/5 rs7203740 sources are not independently repeated; this extension attempts remaining 3 and two rs2767060 tissues. Candidate groups fixed from verified reference annotation before these effect readouts.","Exact C/T single-position gene expression associations; C/A retained separately only for rs7203740. No proxy substitution. Nominal P values descriptive, no confirmatory significance claim. No EAS ancestry-matched molecular analysis, cis-region completeness, model N/sdY, MR or coloc. Missing source/gene is not negative evidence."]}
    for i,(k,rsids) in enumerate(PLAN.items()):
        if i:time.sleep(6)
        study,label,group,n=prev.DATASETS[k];a,b=md[k],paths[k]
        for field,value in {"study_id":study,"study_label":label,"sample_group":group,"sample_size":n,"quant_method":"ge"}.items():require(a[field]==b[field]==value,"Pinned metadata mismatch "+k+" "+field)
        canonical=f"ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{study}/{k}/{k}.all.tsv.gz"
        require(b["ftp_path"]==canonical,"Source path")
        url=canonical.replace("ftp://","https://",1);d={"source":url,"study":label,"sample_group":group,"nominal_dataset_n":int(n),"pmid":a["pmid"],"loci":{}}
        try:
            t,h,total,info=load_index(url);d.update(info)
            for rs in rsids:
                try:d["loci"][rs]=read_exact(url,t,h,total,rs)
                except (OSError,ValueError,KeyError,TypeError,IndexError) as e:d["loci"][rs]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:220]}
        except (OSError,ValueError,KeyError,TypeError,IndexError) as e:
            for rs in rsids:d["loci"][rs]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:220]}
        out["datasets"][k]=d
    out["read_source_locus_pairs"]=sum(z["status"]=="EXACT_POSITION_READ_HTTP_RANGE" for d in out["datasets"].values() for z in d["loci"].values())
    out["attempted_source_locus_pairs"]=5
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence-root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.evidence_root)
    except (OSError,ValueError,KeyError,TypeError) as e:x={"status":"HOLD_PINNED_SOURCE","reason":str(e)[:220]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
