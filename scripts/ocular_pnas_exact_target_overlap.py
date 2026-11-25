#!/usr/bin/env python3
"""Compare all explicitly labeled publisher S2 promoter targets to frozen IDs."""
import argparse,collections,json,math,re
from pathlib import Path
from ocular_frozen_241_external import read_panel,PIN
TABS={"POAG_Gharahkhani":39,"POAG_Han":88,"IOP_Gharahkhani":29,"IOP_Khawaja":37,"IOP_MacGregor":26}
HEADER=["Lead variant","Loci","Chr","Chromosome position (hg38)","ID","Promoter target gene","Promoter length","Distance of variant from promoter"]
LABELS={"Target gene for dexamethasone control":"control","Targte gene dexamethasone treatment":"dex","Target gene dexamethasone treatment":"dex","Target gene resource":"resource"}
def req(ok,msg):
    if not ok:raise ValueError(msg)
def parse(tables):
    records=[]
    section_counts={}
    for tab,n in TABS.items():
        rows=tables[tab];req(len(rows)==n and all(len(v)==8 for v in rows),"ROW_COUNT_OR_WIDTH "+tab)
        req(rows[0][0].startswith("Dataset S2.") and all(v is None for v in rows[0][1:]),"TITLE "+tab)
        section=None;expect_header=False;blocks=collections.Counter();counts=collections.Counter()
        for ix,row in enumerate(rows[1:],2):
            if all(v is None for v in row):continue
            label=row[0].strip() if isinstance(row[0],str) else ""
            if label in LABELS:
                req(not expect_header and all(v is None for v in row[1:]),"SECTION "+tab+":"+str(ix))
                section=LABELS[label];blocks[section]+=1;expect_header=True;continue
            if expect_header:
                req(row==HEADER,"HEADER "+tab+":"+str(ix));expect_header=False;continue
            req(section is not None and not expect_header,"DATA_BEFORE_HEADER "+tab)
            rs,locus,ch,pos,gene,symbol,prom,dist=row
            req(isinstance(rs,str) and re.fullmatch(r"rs[0-9]+",rs) and isinstance(ch,str) and re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)",ch),"VARIANT "+tab+":"+str(ix))
            req(type(pos) is int and pos>0 and type(prom) is int and prom>0 and type(dist) is int,"NUMERIC "+tab+":"+str(ix))
            req(isinstance(gene,str) and re.fullmatch(r"ENSG[0-9]+",gene) and isinstance(symbol,str) and bool(symbol.strip()),"GENE_IDENTITY "+tab+":"+str(ix))
            counts[section]+=1
            records.append({"tab":tab,"phenotype":tab.split("_")[0],"condition":section,"source_row":ix,"rsid":rs,"locus_as_published":locus,"chr":ch,"pos":pos,"gene_id":gene,"published_symbol":symbol,"promoter_length":prom,"distance":dist})
        req(not expect_header and blocks=={"control":1,"dex":1,"resource":1} and all(counts[k]>0 for k in blocks),"INCOMPLETE_SECTIONS "+tab)
        section_counts[tab]=dict(counts)
    return records,section_counts
def analyze(full,evidence,panel_zip):
    d=json.loads(full.read_text());e=json.loads(evidence.read_text())
    req(e["provenance"]["run_id"]==45577966619 and e["provenance"]["artifact_id"]==15004853464,"SOURCE_PROVENANCE")
    req({k:v for k,v in d.items() if k!="all_source_rows"}==e["result"],"INVENTORY_PARITY")
    req(d["status"]=="COMPLETE_PINNED_S2_TABLE_INVENTORY" and d["supplement_source"]["sha256"]=="9197476141d1ca99a40d83fa92137fca48657082f14b9ad5ddcb646babd95cf7","SOURCE_IDENTITY")
    req(set(d["all_source_rows"])==set(TABS),"TAB_IDENTITY")
    for tab,v in d["all_source_rows"].items():
        req(v[:6]==d["tables"][tab]["first_six_rows"] and v[-2:]==d["tables"][tab]["last_two_rows"],"ROWS_PREVIEW_PARITY")
    records,counts=parse(d["all_source_rows"]);panel=read_panel(panel_zip)
    gene=lambda rows:{z["gene_id"] for z in rows}
    rs=lambda rows:{z["rsid"] for z in rows}
    poag=[z for z in records if z["phenotype"]=="POAG"];iop=[z for z in records if z["phenotype"]=="IOP"]
    p,q=gene(poag),gene(iop);u=p|q
    observed={"POAG_genes":len(p),"IOP_genes":len(q),"shared":len(p&q),"union":len(u),"POAG_lead_rsids":len(rs(poag)),"IOP_lead_rsids":len(rs(iop))}
    expected={"POAG_genes":57,"IOP_genes":29,"shared":19,"union":67,"POAG_lead_rsids":46,"IOP_lead_rsids":23}
    result={"status":"PUBLISHED_COUNTS_RECONCILED_EXACT_ID_OVERLAP" if observed==expected else "HOLD_PUBLISHED_COUNTS_MISMATCH","observed":observed,"published":expected,"all_source_records":len(records),"section_counts":counts,"panel_sha256":PIN,"panel_size":len(panel),"ECM_size":sum(v["ecm_member"] for v in panel.values()),"source_s2":d["supplement_source"],"article_xml":d["article_xml"],"groups":{},"limits":["All explicitly labeled resource/control/DEX promoter targets across five S2 tabs; exact stable Ensembl ID only, no symbol alias substitutions. No subset selected to achieve reported counts.","Publisher S2 targets are GWAS-linked promoter contacts, not a DEX DEG panel or proven causal genes. GSE301517 is RNA-seq subseries of same final PNAS study; no independent validation claim.","If counts do not reconcile, report literal source union but HOLD attribution to final 67-gene list. No enrichment significance: a comparable tested gene universe is not established. No claim that absent overlap establishes novelty."]}
    symbol=collections.defaultdict(set)
    for v in records:symbol[v["gene_id"]].add(v["published_symbol"])
    for phenotype,group in [("all",records),("POAG",poag),("IOP",iop)]:
        for cond in ("all","control","dex","resource"):
            rows=[v for v in group if cond=="all" or v["condition"]==cond]
            ids=gene(rows);hits=sorted(ids&set(panel))
            members=[{"gene_id":k,"panel_symbol":panel[k]["symbol"],"publisher_symbols":sorted(symbol[k]),"frozen_direction":"UP" if panel[k]["acute_log2fc"]>0 else "DOWN","ECM16":panel[k]["ecm_member"]} for k in hits]
            result["groups"][phenotype+"_"+cond]={"genes":len(ids),"lead_rsids":len(rs(rows)),"rows":len(rows),"panel_overlap_count":len(hits),"ecm_overlap_count":sum(v["ECM16"] for v in members),"overlap":members}
    result["source_gene_ids"]=sorted(u)
    result["source_symbol_conflicts"]={k:sorted(v) for k,v in symbol.items() if len(v)>1}
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--input",type=Path,required=True);p.add_argument("--evidence",type=Path,required=True);p.add_argument("--panel",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=analyze(a.input,a.evidence,a.panel)
    except (OSError,ValueError,KeyError,TypeError) as ex:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
