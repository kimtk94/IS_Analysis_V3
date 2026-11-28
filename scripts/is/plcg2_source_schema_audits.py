#!/usr/bin/env python3
"""Source audit: exact hash-pinned CS gene IDs and bounded Catalog v2 schema probes."""
import argparse,collections,csv,gzip,hashlib,importlib,json,re
from pathlib import Path
cs=importlib.import_module("scripts.is.plcg2_rosmap_egene_cs_audit")
cat=importlib.import_module("scripts.is.fixed24_two_locus_catalog_audit")
GENE=cs.GENE
CS_SHA="e21efac4551ae7e3558699111481053dcf05deb6557f236fbecb9f330264c360"
CS_URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000025/QTD000434/QTD000434.credible_sets.tsv.gz"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def root_id(s):
    m=re.fullmatch(r"(ENSG[0-9]+)(?:\.[0-9]+)?",s)
    return m[1] if m else None
def audit_ids(work):
    dest=work/"source.gz";out={"status":"HOLD_SOURCE_OR_SCHEMA","expected_source_sha256":CS_SHA}
    try:
        source=cs.download(CS_URL,dest);out["source"]=source
        require(source["sha256"]==CS_SHA and source["compressed_bytes"]==12466333,"Source differs from verified credible set release")
        total=0;types=collections.Counter();genes=set();targets=[];atpos=[];rsid=[]
        with gzip.open(dest,"rt",encoding="utf8",newline="") as f:
            reader=csv.DictReader(f,delimiter="\t")
            require(reader.fieldnames=="molecular_trait_id gene_id cs_id variant rsid cs_size pip pvalue beta se z cs_min_r2 region".split(),"Previously observed schema changed")
            for d in reader:
                total+=1;require(total<=500000 and None not in d and all(x is not None for x in d.values()),"Bounded row schema")
                for k in ("gene_id","molecular_trait_id"):
                    v=d[k];g=root_id(v)
                    types[k+("_exact_ensembl" if g==v else "_versioned_ensembl" if g else "_other")]+=1
                    if k=="gene_id" and g:genes.add(g)
                if root_id(d["gene_id"])==GENE or root_id(d["molecular_trait_id"])==GENE:
                    require(len(targets)<50000,"Target row cap")
                    targets.append({k:d[k] for k in ("gene_id","molecular_trait_id","cs_id","variant","rsid","cs_size","pip")})
                if d["variant"].startswith("chr16_81775211_"):
                    require(len(atpos)<1000,"Position row cap")
                    atpos.append({k:d[k] for k in ("gene_id","molecular_trait_id","cs_id","variant","rsid","pip")})
                if "rs7203740" in re.findall(r"rs[0-9]+",d["rsid"]):
                    require(len(rsid)<1000,"rsID row cap")
                    rsid.append({k:d[k] for k in ("gene_id","molecular_trait_id","cs_id","variant","rsid","pip")})
        require(total==377180,"Previous full source row count mismatch")
        out.update(status="COMPLETE_SAME_RELEASE_IDENTIFIER_AUDIT",source_rows=total,identifier_counts=dict(types),distinct_ensembl_gene_roots=len(genes),target_exact_or_version_suffix_rows=len(targets),target_rows=targets,exact_position_all_gene_rows=atpos,rs7203740_all_gene_rows=rsid,limits=["Same official hash-pinned release, not independent replication. Normalization strips only a numeric Ensembl version suffix; it does not infer aliases or map other identifier systems. No absent variant PIP assigned."])
    except (OSError,ValueError,TypeError,KeyError,EOFError) as e:out["reason"]=str(e)[:350]
    finally:dest.unlink(missing_ok=True)
    return out
def catalog():
    out={"status":"BOUNDED_LIVE_SCHEMA_PROBE_NOT_NOVELTY","responses":{},"limits":["Probe uses small first pages only; no assertion of complete coverage, association absence or novelty. Positive control rs2129981 included to validate exact filter and schema."]}
    md,src=cat.fetch("metadata")
    out["metadata"]={k:md.get(k) for k in ("version","data_release_date","gene_build","dbsnp_build")};out["metadata_source"]=src
    probes=[
        ("chr10_region","single-nucleotide-polymorphisms",{"chromosome":"10","bp_start":35333363,"bp_end":36333363}),
        ("chr16_region","single-nucleotide-polymorphisms",{"chromosome":"16","bp_start":81275211,"bp_end":82275211}),
        ("rs2767060","associations",{"rs_id":"rs2767060"}),
        ("rs7203740","associations",{"rs_id":"rs7203740"}),
        ("known_rs2129981","associations",{"rs_id":"rs2129981"})]
    for name,endpoint,params in probes:
        try:
            x,src=cat.fetch(endpoint,{**params,"page":0,"size":2})
            if isinstance(x,dict):
                description={k:type(v).__name__ for k,v in x.items()}
                # Retain the actual compact page and one or two complete rows;
                # never treat first-page counts as all records.
                snapshot=x if len(json.dumps(x))<=60000 else {"keys":list(x),"field_types":description}
            elif isinstance(x,list):
                description={"type":"list","length":len(x)};snapshot=x[:2]
            else:description={"type":type(x).__name__};snapshot=str(x)[:500]
            out["responses"][name]={"status":"HTTP_RESPONSE","source":src,"schema":description,"snapshot":snapshot}
        except (OSError,ValueError,KeyError,TypeError) as e:out["responses"][name]={"status":"HOLD_SOURCE","reason":str(e)[:350]}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["identifiers","catalog"]);p.add_argument("--work",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    try:x=audit_ids(a.work) if a.mode=="identifiers" else catalog()
    except (OSError,ValueError,KeyError,TypeError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:350]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
