#!/usr/bin/env python3
"""Source-pinned Ng 2022 Fig2b original gene-wise rho: runner-private ONLY.

Never print/export gene IDs or per-gene effect values. Exact paper-side Ensembl
set digest and original sample-free rho median are audited in public aggregate.
"""
import argparse,hashlib,io,json,pathlib,re,statistics,urllib.request,xml.etree.ElementTree as ET,zipfile
from audit_ng_stage27_original_figure2_source import URL,SHA
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RX=re.compile(r"^ENSG[0-9]{11}(?:\.[0-9]+)?$")
SET_DIGEST="61ac8fabc8144fc29e46f1d932f4400d86e22f533173743b678ab4dfcce356a7"
def fetch():
    req=urllib.request.Request(URL,headers={"User-Agent":"HCC-Stage29-original-figure2b-gene-effects/1.0"})
    with urllib.request.urlopen(req,timeout=60) as f:b=f.read(7_000_000)
    if hashlib.sha256(b).hexdigest()!=SHA:raise ValueError("source paper XLSX SHA256 mismatch")
    return b
def parse(b):
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        names=set(z.namelist())
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        relationships=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        m={}
        for r in relationships:
            match=re.search(r"(?:^|/)sheet(\d+)\.xml$",r.get("Target",""))
            if match:m[r.get("Id")]=int(match.group(1))
        rid="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        sheet={s.get("name"):m.get(s.get(rid)) for s in wb.iter(NS+"sheet")}
        if sheet.get("Fig2b")!=8 or sheet.get("Fig2a")!=7:
            raise ValueError("original Excel workbook relationship drift")
        shared=[]
        if "xl/sharedStrings.xml" in names:
            with z.open("xl/sharedStrings.xml") as f:
                for _,e in ET.iterparse(f,events=("end",)):
                    if e.tag==NS+"si":
                        shared.append("".join(t.text or "" for t in e.iter(NS+"t")))
                        e.clear()
        def value(c):
            typ=c.get("t")
            if typ=="inlineStr":return "".join(t.text or "" for t in c.iter(NS+"t"))
            e=c.find(NS+"v");v=e.text if e is not None else ""
            if typ=="s" and v is not None and v.isdecimal():
                j=int(v);return shared[j] if j<len(shared) else ""
            return v or ""
        entries=[]
        header=[]
        with z.open("xl/worksheets/sheet8.xml") as f:
            for _,row in ET.iterparse(f,events=("end",)):
                if row.tag!=NS+"row":continue
                cells={}
                for c in row.findall(NS+"c"):
                    idx=re.match(r"[A-Z]+",c.get("r",""))
                    if idx and idx.group() in ("A","B","C","D"):
                        cells[idx.group()]=value(c)
                if not header:
                    header=[cells.get(x,"") for x in ("A","B","C","D")]
                else:
                    g=cells.get("A","").strip()
                    if not RX.fullmatch(g):
                        raise ValueError("Fig2b first-column source Ensembl ID missing")
                    g=g.split(".")[0]
                    try:r=float(cells.get("B",""))
                    except ValueError:raise ValueError("Fig2b published Spearman rho is not numeric") from None
                    if not -1<=r<=1:raise ValueError("published rho outside [-1,1]")
                    entries.append((g,r))
                row.clear()
    ids=[g for g,_ in entries]
    if len(ids)!=5481 or len(set(ids))!=5481:
        raise ValueError("not exact 5481 unique original Figure2b ENSG rows")
    digest=hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()
    if digest!=SET_DIGEST:
        raise ValueError("published Figure2b Ensembl digest drift")
    median=statistics.median(r for _,r in entries)
    if abs(median-0.287059)>1e-5:
        raise ValueError("published original Figure2b rho median changed")
    if not any(w in header[1].lower() for w in ("spearman","correl","rho")):
        raise ValueError("Fig2b source correlation header not explicit")
    return entries,{"schema":"hcc-stage29-original-fig2b-rho-extraction-v1",
        "source_xlsx_sha256":SHA,"source_sheet":"Fig2b",
        "physical_source_worksheet_index":8,
        "published_ENSG_n":5481,
        "source_gene_set_sha256":digest,
        "source_spearman_median":round(median,9),
        "published_rho_numeric_n":len(entries),
        "source_gene_and_per_gene_effect_values_exported":False,
        "original_Ng_2022_quantitation_reproduced":False}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-rhos-out",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    rows,report=parse(fetch())
    target=pathlib.Path(a.private_rhos_out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text("".join(f"{g}\t{r:.16g}\n" for g,r in sorted(rows)),encoding="utf-8")
    target.chmod(0o600)
    path=pathlib.Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
