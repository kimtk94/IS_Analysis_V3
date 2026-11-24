#!/usr/bin/env python3
"""Compare exact reported Japanese conditional cis records with published stroke fine-mapping tables.

Selected rows, not full QTL or causal colocalization. Preserve all publisher
table identities and exact matched rsIDs, without interpreting unlabeled fields.
"""
import argparse,collections,gzip,hashlib,io,json,math,re,zipfile
from pathlib import Path

ZIP_SHA="b878a51af3d58311734bc3c84a0af1941d6bfe3b67548416d071c031286206ee"
WORKBOOK_SHA="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
GENES={"ENSG00000138172":"CALHM2","ENSG00000107957":"SH3PXD2A","ENSG00000107954":"NEURL1"}
HEADER=("Gene_id","Gene_name","CHR","TSS_position","Number_of_variants_cis","Variant_ID","Variant_CHR",
"Variant_position_start","Variant_position_end","Rank_of_association","Forward_nominal_P","Forward_slope","Backward_P","Backward_slope")
KNOWN=("rs4918058","rs55983834","rs1565706","rs7091346","rs2295786")
def digest(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def require(b,msg):
    if not b:raise ValueError(msg)
def member(z,name):
    hits=[s for s in z.namelist() if s.rsplit("/",1)[-1]==name]
    require(len(hits)==1,"Missing archive member "+name)
    return hits[0]
def qtl(path):
    require(digest(path)==ZIP_SHA,"Previous public artifact digest changed")
    d=collections.defaultdict(lambda:collections.defaultdict(list))
    counts=collections.Counter()
    with zipfile.ZipFile(path) as z:
        require(z.testzip() is None,"Archive CRC")
        meta=json.loads(z.read(member(z,"japanese.json")))
        require(meta["member_count"]==28 and meta["total_source_lines"]==33193329,"Original source identity")
        members={x["member"]:x for x in meta["members"]}
        require(all(tuple(x["first_lines"][0].split("\t"))==HEADER for x in members.values()),"Published conditional header changed")
        with z.open(member(z,"conditional_candidate_rows.tsv.gz")) as f,gzip.GzipFile(fileobj=f) as gz,io.TextIOWrapper(gz,encoding="utf-8") as t:
            require(t.readline().strip()=="member\toriginal_summary_row","Derived table header")
            for l in t:
                a=l.rstrip("\r\n").split("\t")
                require(len(a)==15 and a[0] in members,"Invalid conditional row")
                r=dict(zip(HEADER,a[1:]));g=GENES.get(r["Gene_id"].split(".")[0])
                if g is None:continue
                require(r["Gene_name"]==g and r["CHR"]==r["Variant_CHR"]=="chr10","Gene identity mismatch")
                rs=r["Variant_ID"]
                p=float(r["Backward_P"]);pos=int(r["Variant_position_start"]);end=int(r["Variant_position_end"]);rank=int(r["Rank_of_association"])
                require(math.isfinite(p) and 0<=p<=1 and rank>=0 and pos>0 and end>=pos,"Invalid source association")
                counts[g]+=1
                if not re.fullmatch(r"rs[0-9]+",rs):continue
                d[g][rs].append({"cell":a[0].split("_conditional")[0],"position":pos,"end":end,"rank":rank,"conditional_p":p})
    require(dict(counts)=={"CALHM2":3227,"SH3PXD2A":351,"NEURL1":290},"Fixed gene row totals")
    for g,rs in d.items():
        require(all(len({(x["position"],x["end"]) for x in rr})==1 for rr in rs.values()),"Ambiguous source rsID coordinate")
    return d,counts
def exact_tokens(row):
    # Full-cell rsID tokens only (allow semicolon/comma/whitespace lists), never substring matching.
    out=set()
    for v in row:
        if not isinstance(v,str):continue
        toks=re.split(r"[\s;,]+",v.strip())
        out.update(t for t in toks if re.fullmatch(r"rs[0-9]+",t))
    return out
def main():
    from openpyxl import load_workbook
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    q,counts=qtl(a.root/"qtl.zip")
    wpath=a.root/"supp.xlsx";require(wpath.stat().st_size==756910 and digest(wpath)==WORKBOOK_SHA,"Previously verified publisher XLSX differs")
    w=load_workbook(wpath,read_only=True,data_only=True)
    inventory=[];found={};hits=[];gene_mentions=[]
    union=set().union(*(set(v) for v in q.values()))
    lookup=union|set(KNOWN)
    for sh in w.worksheets:
        m=re.fullmatch(r"Suppl[.]?\s*Table[.]?\s*(17|18|19|27|38)(?:[a-zA-Z])?",sh.title,re.I)
        if not m:continue
        number=int(m.group(1));found.setdefault(number,[]).append(sh.title)
        first=[];matched_rs=set();n=0;selected=[]
        for i,row in enumerate(sh.iter_rows(values_only=True),1):
            n+=1
            fields=[str(x) if x is not None else "" for x in row]
            if i<=8:first.append(fields[:45])
            if number in (17,18,19):
                rs=exact_tokens(row)&lookup
                if rs:
                    matched_rs.update(rs)
                    entry={"sheet":sh.title,"row":i,"rsids":sorted(rs),"source_fields":fields[:60]}
                    selected.append(entry);hits.append(entry)
            elif number in (27,38):
                genes={g for g in GENES.values() if g in fields} | ({"NEURL1"} if "NEURL" in fields else set())
                if genes:gene_mentions.append({"sheet":sh.title,"row":i,"genes":sorted(genes),"source_fields":fields[:60]})
        inventory.append({"sheet":sh.title,"table":number,"rows_scanned":n,"first_rows":first,"matched_unique_rsids":sorted(matched_rs),"matched_rows":len(selected)})
    require({17,18,19,27,38}<=set(found),"Published table name mismatch "+str(w.sheetnames))
    w.close()
    overlap=[]
    for sh in inventory:
        if sh["table"] not in (17,18,19):continue
        sr=set(sh["matched_unique_rsids"])
        for gene,rr in q.items():
            ids=sorted(sr&set(rr))
            overlap.append({"sheet":sh["sheet"],"gene":gene,"unique_rsid_overlap":len(ids),"rsids":ids,
                            "cell_types_with_overlap":sorted({r["cell"] for rs in ids for r in rr[rs]}),
                            "minimum_conditional_p":min((r["conditional_p"] for rs in ids for r in rr[rs]),default=None)})
    pairs=[]
    for i,g in enumerate(sorted(q)):
        for h in sorted(q)[i+1:]:
            shared=set(q[g])&set(q[h])
            pairs.append({"gene1":g,"gene2":h,"shared_rsid_count":len(shared),"shared_rsids":sorted(shared),
                          "shared_same_cell_count":sum(bool({x["cell"] for x in q[g][rs]}&{x["cell"] for x in q[h][rs]}) for rs in shared)})
    known={g:{rs:q[g][rs] for rs in KNOWN if rs in q[g]} for g in q}
    out={"status":"PUBLISHED_FINE_MAPPING_TABLE_RSID_OVERLAP_NOT_COLOC","source_workbook_sha256":WORKBOOK_SHA,"qtl_artifact_sha256":ZIP_SHA,
         "source_gene_rows":dict(counts),"unique_rsid_counts":{g:len(v) for g,v in q.items()},"gene_pairs":pairs,
         "tables":inventory,"overlap_by_sheet_gene":overlap,"matched_source_rows":hits,"prior_gene_mentions":gene_mentions,
         "known_lead_conditional_records":known,
         "limits":["Exact rsID overlap in selected publisher sheets; inspect actual column, locus, trait, ancestry and credible-set membership before calling any row a credible-set variant.",
                   "No QTL effect allele or SE; no harmonized effect comparison, new significance, LD/PIP weighting or causal colocalization.",
                   "Same donor sets across immune cells; different gene-specific tested variant universes; absence from selected conditional source is not absence of regulation.",
                   "Published gene rows in Tables 27/38 require interpretation of analysis type and phenotype, not automatic novelty claim."]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
