#!/usr/bin/env python3
"""Source-preserving OASIS Japan Omics Bokeh QTL audit, not a coloc substitute."""
import argparse,base64,csv,hashlib,json,re,urllib.request
from collections import Counter
from pathlib import Path
import numpy as np
GENES={"ALDH2":"ENSG00000111275","BRAP":"ENSG00000089234","RPH3A":"ENSG00000089169"}
TARGETS={"rs11066015","rs671","rs11066132","rs77768175"}
def parse(html,gene):
 m=re.search(r'<script[^>]*type="application/json"[^>]*>(.*?)</script>',html,re.S)
 if not m:raise ValueError("Missing embedded Bokeh data")
 obj=json.loads(m.group(1)); sources=[]
 def walk(x):
  if isinstance(x,dict):
   if x.get("name")=="ColumnDataSource" and x.get("attributes",{}).get("data",{}).get("entries"):
    sources.append(x["attributes"]["data"]["entries"])
   for y in x.values():walk(y)
  elif isinstance(x,list):
   for y in x:walk(y)
 walk(obj)
 expected={"variant_id_hg38","rsid","gene_name","pval_nominal","effect_size","effect_size_SE","maf"}
 source=None
 for s in sources:
  if expected.issubset({x[0] for x in s}):
   if source is not None:raise ValueError("Multiple matching data sources")
   source=s
 if source is None:raise ValueError("No complete OASIS association source")
 columns={}
 for name,v in source:
  if v.get("type")!="ndarray":continue
  arr=v["array"]
  if isinstance(arr,list):cols=arr
  else:
   dt=v["dtype"];endian="<" if v.get("order")=="little" else ">"
   raw=base64.b64decode(arr["data"],validate=True)
   cols=np.frombuffer(raw,dtype=np.dtype(endian+{"float64":"f8","int32":"i4","int64":"i8","float32":"f4","uint32":"u4"}.get(dt,"f8"))).tolist()
  columns[name]=cols
 n=len(columns["variant_id_hg38"])
 if any(len(col)!=n for col in columns.values()):raise ValueError("Bokeh column misalignment")
 rows=[{key:columns[key][i] for key in expected} for i in range(n)]
 if any(row["gene_name"]!=gene for row in rows):raise ValueError("Unexpected gene")
 if len({row["variant_id_hg38"] for row in rows})!=n:raise ValueError("Duplicate variant in cell-page source")
 return rows
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--outdir",required=True,type=Path)
 ap.add_argument("--cells",default="Mono-L1,B_Activated-L2")
 a=ap.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
 manifests=[]; allrows=[]
 for gene,ensg in GENES.items():
  for cell in a.cells.split(","):
   url=f"https://japan-omics.jp/gene/OASIS-{cell}?input_value={ensg}"
   snap=a.outdir/f"OASIS_{gene}_{cell}_original.html"
   if snap.exists():
    content=snap.read_bytes()
   else:
    req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req,timeout=22) as res:
     content=res.read(2_000_001)
   if len(content)>2_000_000:raise ValueError("Page too big")
   html=content.decode("utf8","replace")
   try:
    rows=parse(html,gene)
    source_status="BOUNDED_BOKEH_ROWS_PARSED"
   except ValueError as e:
    if "Missing embedded Bokeh data" not in str(e):raise
    rows=[]
    source_status="NO_EMBEDDED_BOKEH_DATA_NOT_BIOLOGICAL_NEGATIVE"
   snap=a.outdir/f"OASIS_{gene}_{cell}_original.html"
   if not snap.exists():snap.write_bytes(content)
   found=[r for r in rows if r["rsid"] in TARGETS]
   non_sig=sum(float(r["pval_nominal"])>.05 for r in rows)
   manifests.append(dict(gene=gene,cell=cell,n=len(rows),non_sig=non_sig,source_status=source_status,
                         target_rsid_found=[r["rsid"] for r in found],
                         sha256=hashlib.sha256(content).hexdigest(),url=url))
   for r in found:allrows.append(dict(gene=gene,cell=cell,**r))
   print(gene,cell,"N",len(rows),"p_gt05",non_sig,"TARGETS",[(r["rsid"],r["pval_nominal"]) for r in found],flush=True)
 out={"audit":"G0022_OASIS_GENE_BOKEH_V1","sources":manifests,
      "associated_four_cs_rows":allrows,
      "science_gate":"BROWSER_SOURCE_AUDIT_NOT_VALIDATED_FULL_CIS_NOT_COLOC",
      "full_locus_coloc_ready":False,"causal_gene_confirmed":False,
      "warning":"Individual OASIS cell pages may filter/sort SNPs; non-significant rows on a page do not prove all cis variants included. Missing target on browser is not biological negative."}
 (a.outdir/"G0022_OASIS_BROWSER_AUDIT.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
if __name__=="__main__":main()
