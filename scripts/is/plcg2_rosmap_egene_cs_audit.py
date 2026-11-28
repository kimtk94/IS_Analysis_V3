#!/usr/bin/env python3
"""Official ROSMAP PLCG2 eGene/credible set evidence gate; no new coloc."""
import argparse,collections,csv,gzip,hashlib,importlib,json,math,re,time
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen
meta=importlib.import_module("scripts.is.audit_official_p3a_qtl_catalogue_provenance")
GENE="ENSG00000197943"
VAR="chr16_81775211_C_T"
DATASET="QTD000434"
MAX_BYTES=40000000
def require(ok,msg):
    if not ok:raise ValueError(msg)
def download(url,dest):
    require(urlparse(url).scheme=="https" and urlparse(url).hostname=="ftp.ebi.ac.uk","Unapproved source host")
    h=hashlib.sha256();n=0
    with urlopen(Request(url,headers={"Accept-Encoding":"identity","User-Agent":"IS-ROSMAP-official-eGene-CS-audit/1.0"}),timeout=65) as r,dest.open("wb") as f:
        require(r.status==200 and urlparse(r.geturl()).hostname=="ftp.ebi.ac.uk","Unexpected source status/host")
        size=r.headers.get("Content-Length")
        if size:require(0<int(size)<=MAX_BYTES,"Declared file exceeds cap")
        while True:
            b=r.read(262144)
            if not b:break
            n+=len(b);require(n<=MAX_BYTES,"File exceeds cap");h.update(b);f.write(b)
        if size:require(n==int(size),"Truncated source")
    require(n>100,"Empty source")
    return {"url":url,"compressed_bytes":n,"sha256":h.hexdigest(),"complete_file_retrieved":True}
def read_perm(path):
    data={};target=None
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        t=csv.DictReader(f,delimiter="\t")
        require({"molecular_trait_object_id","molecular_trait_id","n_variants","variant","chromosome","position","pvalue","p_beta","p_perm","beta"}<=set(t.fieldnames or []),"Official permutation schema")
        for row in t:
            require(len(data)<60000 and None not in row,"Permutation row cap/width")
            gid=row["molecular_trait_object_id"];require(re.fullmatch(r"ENSG\d+(?:\.\d+)?",gid) and row["molecular_trait_id"]==gid and gid not in data,"Expected unique gene-expression traits")
            p=float(row["p_beta"]);require(math.isfinite(p) and 0<=p<=1,"Invalid beta-approximated permutation P")
            data[gid]=p
            if gid==GENE:
                target={k:row[k] for k in ("molecular_trait_id","n_variants","variant","chromosome","position","pvalue","p_beta","p_perm","beta")}
                require(row["chromosome"].removeprefix("chr")=="16" and int(row["n_variants"])>0,"PLCG2 gene identity")
    require(len(data)>=5000,"Not a complete gene-expression permutation universe")
    vals=sorted(data.items(),key=lambda x:(x[1],x[0]))
    adj={};running=1.
    for i,(g,p) in reversed(list(enumerate(vals,1))):
        running=min(running,p*len(vals)/i,1.);adj[g]=running
    if target:target["BH_across_complete_source_genes"]=adj[GENE];target["passes_BH_0.05"]=adj[GENE]<.05
    return {"status":"COMPLETE_OFFICIAL_PERMUTATION_READ","tested_gene_count":len(data),"BH_lt_0.05_gene_count":sum(v<.05 for v in adj.values()),"PLCG2":target}
def read_cs(path):
    """Use the observed release header, not the older documented CS columns."""
    sets=collections.defaultdict(list);total=0;seen=set();skipped_other_genes=0
    fields="molecular_trait_id gene_id cs_id variant rsid cs_size pip pvalue beta se z cs_min_r2 region".split()
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        t=csv.DictReader(f,delimiter="\t")
        require(t.fieldnames==fields,"Unexpected observed-release CS schema: "+repr(t.fieldnames))
        for row in t:
            total+=1;require(total<=2500000 and None not in row and all(v is not None for v in row.values()),"CS row bound/width")
            if row["molecular_trait_id"]!=GENE and row["gene_id"]!=GENE:
                skipped_other_genes+=1;continue
            require(row["gene_id"]==row["molecular_trait_id"]==GENE,"Inconsistent target gene identity")
            m=re.fullmatch(r"chr16_([0-9]+)_([ACGTN]+)_([ACGTN]+)",row["variant"])
            require(m is not None,"Noncanonical target GRCh38 variant key: "+row["variant"][:120])
            pos,ref,alt=int(m[1]),m[2],m[3]
            require(pos>0 and ref!=alt,"Target coordinate/alleles")
            pip,p,beta,se,z,lo=map(float,[row["pip"],row["pvalue"],row["beta"],row["se"],row["z"],row["cs_min_r2"]])
            require(all(math.isfinite(v) for v in (pip,p,beta,se,z,lo)) and 0<=pip<=1 and 0<=p<=1 and se>0 and 0<=lo<=1,"CS effect/LD statistics")
            key=(row["cs_id"],row["variant"]);require(key not in seen and bool(row["cs_id"]),"Duplicate or empty CS identity");seen.add(key)
            require(int(row["cs_size"])>0 and bool(row["region"]),"CS size/region")
            if row["variant"]==VAR:
                require("rs7203740" in re.findall(r"rs[0-9]+",row["rsid"]),"Target C/T rsID mismatch")
                require(math.isclose(beta,.279006,rel_tol=1e-5) and math.isclose(se,.068636,rel_tol=1e-5) and math.isclose(p,.0000551034,rel_tol=1e-5),"Target CS nominal values disagree with pinned QTL")
            sets[row["cs_id"]].append({"variant":row["variant"],"rsid":row["rsid"],"position":pos,"ref":ref,"alt":alt,"pip":pip,"p":p,"beta_ALT":beta,"se":se,"z":z,"cs_size":int(row["cs_size"]),"cs_min_r2":lo,"region_as_published":row["region"]})
    require(total>0,"Empty credible sets source")
    summary=[]
    for k,rows in sorted(sets.items()):
        require(all(v["cs_size"]==len(rows) for v in rows),"Incomplete PLCG2 credible set: "+str([k,len(rows),sorted({v["cs_size"] for v in rows})]))
        require(len({v["cs_min_r2"] for v in rows})==1 and len({v["region_as_published"] for v in rows})==1,"Inconsistent CS metadata")
        lead=max(rows,key=lambda v:(v["pip"],v["variant"]))
        summary.append({"cs_id":k,"size":len(rows),"top_pip_variant":lead["variant"],"top_pip_rsid":lead["rsid"],"top_pip":lead["pip"],"target_CT":[v for v in rows if v["variant"]==VAR],"other_CA":[v for v in rows if v["variant"]=="chr16_81775211_C_A"],"cs_min_r2":lead["cs_min_r2"],"region_as_published":lead["region_as_published"],"members":rows})
    return {"status":"COMPLETE_OFFICIAL_CREDIBLE_SET_READ","actual_header":fields,"total_source_rows":total,"other_gene_rows":skipped_other_genes,"target_gene_rows":total-skipped_other_genes,"PLCG2_credible_sets":summary,"PLCG2_credible_set_count":len(summary),"target_CT_in_any_CS":any(s["target_CT"] for s in summary),"unavailable_in_release":["cs_avg_r2","cs_log10bf","cs_index","posterior_mean","posterior_sd"],"limits":["Coordinates and alleles parsed from official GRCh38 variant key; no inferred missing CS statistics. PIP and LD purity are source-reported. Absence from a published CS is not proof of no molecular effect."]}

def run(root,work):
    p=json.loads((root/"is_fixed24_two_locus_qtl_extension_20261125.json").read_text())
    x=p["result"];require(p["provenance"]["run_id"]==45549445784 and x["read_source_locus_pairs"]==5,"Prior extension identity")
    q=x["datasets"][DATASET]["loci"]["rs7203740"]
    h=[z for z in q["candidate_rows"] if z["gene_id"]==GENE and z["variant"]==VAR]
    require(len(h)==1 and math.isclose(h[0]["p"],.0000551034,rel_tol=1e-8),"Prior exact nominal row")
    d=meta.get_official_sources();md=meta.table(d["meta"])[DATASET];paths=meta.table(d["paths"])[DATASET]
    require(md["study_label"]=="ROSMAP" and md["sample_group"]=="brain_naive" and md["quant_method"]=="ge" and md["sample_size"]=="560","Official dataset identity")
    require("Column names of the permutation p-value files" in d["columns"] and "fine mapped credible set files" in d["columns"],"Pinned schema")
    nominal="ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz"
    cs="ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000025/QTD000434/QTD000434.credible_sets.tsv.gz"
    require(paths["ftp_path"]==nominal and paths["ftp_cs_path"]==cs,"Official source paths")
    # Permutation filename follows documented source naming; 404 is HOLD.
    perm=nominal.replace(".all.tsv.gz",".permuted.tsv.gz")
    out={"status":"POST_READOUT_OFFICIAL_PLCG2_EGENE_CS_AUDIT","dataset":DATASET,"gene_id":GENE,"target_variant":VAR,"nominal_qtl":h[0],"sources":{},"limits":["Official gene-level permutation/FDR and SuSiE credible-set evidence, not independent replication; same ROSMAP QTL dataset as nominal result.","Exact C/T not rsID alone. Presence in QTL credible set does not imply shared causal variant with stroke GWAS. Absence can reflect no published set or signal structure; not absence of all regulation.","No regional GWAS/QTL comparison, MR or colocalization. Nominal N/AN not automatically model N; ancestry and overlap not certified. Source failure remains HOLD."]}
    # Reuse the independently archived exact permutation readback; do not redownload.
    archived=json.loads((root/"is_plcg2_rosmap_egene_cs_20261125.json").read_text())
    q=archived["result"]["sources"]["permutation"]
    require(archived["provenance"]["run_id"]==45550125894 and q["status"]=="COMPLETE_OFFICIAL_PERMUTATION_READ" and q["tested_gene_count"]==23076 and q["source"]["sha256"]=="adf13cbe9f784cd56a832700b710ffd6abe76a962a38918667a1d80c51413091","Prior permutation provenance")
    out["sources"]["permutation"]={**q,"reused_verified_run":45550125894}
    for k,url,parser in (("credible_sets",cs,read_cs),):
        dst=work/(k+".gz");url=url.replace("ftp://","https://",1)
        try:
            source=download(url,dst)
            out["sources"][k]={"source":source}
            out["sources"][k].update(parser(dst))
        except (OSError,ValueError,KeyError,TypeError,EOFError) as e:out["sources"][k]={**out["sources"].get(k,{}),"status":"HOLD_SOURCE_OR_SCHEMA","url":url,"reason":str(e)[:1200]}
        finally:dst.unlink(missing_ok=True)
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence-root",type=Path,required=True);p.add_argument("--work",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    try:x=run(a.evidence_root,a.work)
    except (OSError,ValueError,KeyError,TypeError) as e:x={"status":"HOLD_PINNED_SOURCE","reason":str(e)[:220]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
