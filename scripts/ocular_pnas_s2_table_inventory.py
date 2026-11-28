#!/usr/bin/env python3
"""Inspect exact source S2 headers/rows and article definitions before overlap.
Publisher PNAS CC BY; no individual genotype data; no automatic 67-member inference.
"""
import argparse,collections,hashlib,io,json,re
from pathlib import Path
from xml.etree import ElementTree as ET
import openpyxl
import ocular_final_publication_inventory as source
EXPECTED_TABS={"POAG_Gharahkhani":39,"POAG_Han":88,"IOP_Gharahkhani":29,"IOP_Khawaja":37,"IOP_MacGregor":26}
def req(ok,msg):
    if not ok:raise ValueError(msg)
def cell(v):
    if v is None:return None
    if isinstance(v,(str,int,float,bool)):return v
    return str(v)
def run(evidence):
    x=json.loads(evidence.read_text())
    req(x["provenance"]["run_id"]==45577692199 and x["provenance"]["artifact_id"]==15004134568,"Inventory provenance")
    a=x["result"]["papers"]["PNAS"]["articles"]
    req(len(a)==1 and a[0]["prefix"]=="PMC13778663.1/" and a[0]["metadata"]["license_code"]=="CC BY","Article identity/license")
    a=a[0]
    obj=a["xml_source"];raw,meta=source.fetch(obj["url"],5000000)
    req(meta["sha256"]==obj["sha256"] and hashlib.md5(raw).hexdigest()==obj["published_md5"],"XML source drift")
    root=ET.fromstring(raw)
    # Keep only short author-defined clause excerpts about Dataset S2/target counting.
    snippets=[]
    for p in root.findall(".//body//p"):
        v=" ".join("".join(p.itertext()).split())
        if not re.search(r"(?:Dataset\s+S0?2|67\s+(?:genes|target genes)|29\s+(?:target\s+)?genes|57\s+(?:target\s+)?genes|19\s+(?:genes|shared))",v,re.I):continue
        for m in re.finditer(r"(?:Dataset\s+S0?2|67\s+(?:genes|target genes)|29\s+(?:target\s+)?genes|57\s+(?:target\s+)?genes|19\s+(?:genes|shared))",v,re.I):
            # At most 25 words per snippet for human source interpretation; no entire paragraph.
            start=max(0,m.start()-100);end=min(len(v),m.end()+100)
            frag=" ".join(v[start:end].split()[:25])
            if frag not in snippets:snippets.append(frag)
            if len(snippets)>=16:break
        if len(snippets)>=16:break
    matches=[w for w in a["workbooks"] if w.get("source",{}).get("key","").endswith(".sd02.xlsx")]
    req(len(matches)==1,"S2 identity")
    obj=matches[0]["source"];req(obj["sha256"]=="9197476141d1ca99a40d83fa92137fca48657082f14b9ad5ddcb646babd95cf7","Pinned S2 hash")
    raw,ms=source.fetch(obj["url"],1000000)
    req(ms["sha256"]==obj["sha256"] and hashlib.md5(raw).hexdigest()==obj["published_md5"],"S2 drift")
    book=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    allrows={};tables={}
    try:
        req(set(book.sheetnames)==set(EXPECTED_TABS),"Sheet identities")
        for tab,expected in EXPECTED_TABS.items():
            sheet=book[tab];req(sheet.max_row==expected and sheet.max_column==8,"Sheet dimensions "+tab)
            rows=[[cell(v) for v in row] for row in sheet.iter_rows(min_row=1,max_row=expected,max_col=8,values_only=True)]
            req(len(rows)==expected and "Dataset S2." in str(rows[0][0]),"Title/row count")
            allrows[tab]=rows
            columns=[]
            for j in range(8):
                vals=[r[j] for r in rows[4:] if r[j] is not None]
                unique=collections.Counter(str(v) for v in vals)
                columns.append({"column":j+1,"nonempty_after_row4":len(vals),"unique":len(unique),"all_values_if_leq_15":dict(unique) if len(unique)<=15 else None,"examples":list(unique)[:4]})
            tables[tab]={"rows":expected,"first_six_rows":rows[:6],"last_two_rows":rows[-2:],"columns":columns}
    finally:book.close()
    return {"status":"COMPLETE_PINNED_S2_TABLE_INVENTORY","article_xml":meta,"supplement_source":ms,"author_definition_excerpts":snippets,"tables":tables,"all_source_rows":allrows,"source_lineage":{"PNAS_final":"GSE301525","RNA_subseries":"GSE301517","GEO_super_series":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301525","IOVS_corrected_PRJNA":"PRJNA1243143","IOVS_corrected_source":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12882097/"},"limits":["All five S2 sheets preserved verbatim as public CC BY source data with attribution (Singh et al., PNAS 2026 DOI 10.1073/pnas.2609070123; NLM PMC Cloud). No source gene set inferred from worksheet names, row counts or the abstract alone. Need explicit interpretation of control versus DEX columns and author gene inclusion rules before exact overlap.","No donor-level data or genetic records. No novelty result, no inferential enrichment test, no independent validation."]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.evidence)
    except (OSError,ValueError,KeyError,TypeError,ET.ParseError) as ex:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in x.items() if k!="all_source_rows"},indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
