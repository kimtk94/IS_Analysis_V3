#!/usr/bin/env python3
"""GSE65240 blinded sample-label audit; no expression response readout."""
import argparse,collections,csv,gzip,hashlib,io,json,re,time
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlparse,unquote
from html.parser import HTMLParser
EXPECTED={f"GSM159060{3+2*i+j}":(lot,arm) for i,lot in enumerate(("2584","3423","4973")) for j,arm in enumerate(("control","dex"))}
PLATS=("GPL14550","GPL17077")
ARM={"control":"control","ctrl":"control","cont":"control","untreated":"control","naive":"control","dex":"dex"}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def fetch(url,cap):
    host=urlparse(url).hostname
    require(host in {"ftp.ncbi.nlm.nih.gov","www.ebi.ac.uk"} and urlparse(url).scheme=="https","Unapproved source")
    with urlopen(Request(url,headers={"Accept-Encoding":"identity","User-Agent":"GSE65240-blinded-sample-label-audit/1.0"}),timeout=45) as r:
        require(r.status==200 and urlparse(r.geturl()).hostname==host,"Unexpected source response")
        raw=r.read(cap+1)
    require(0<len(raw)<=cap,"Source byte bound")
    return raw,{"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)}
def meta(raw):
    out=collections.defaultdict(list);n=0
    with gzip.open(io.BytesIO(raw),"rt",encoding="utf-8-sig",newline="") as f:
        for line in f:
            n+=1;require(n<=10000,"Metadata line cap")
            if line.strip()=="!series_matrix_table_begin":return out
            if line.startswith("!"):
                row=next(csv.reader([line.rstrip("\r\n")],delimiter="\t"))
                if row:out[row[0]].append(row[1:])
    raise ValueError("No matrix metadata boundary")
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag=="a":
            href=dict(attrs).get("href")
            if href:self.links.append(href)
def filename(name,gsm):
    name=unquote(name.rsplit("/",1)[-1])
    m=re.fullmatch(re.escape(gsm)+r"_(2584|3423|4973)_([A-Za-z]+)\.txt\.gz",name,re.I)
    if not m:return None
    arm=ARM.get(m[2].lower())
    return {"filename":name,"lot":m[1],"arm":arm}
def run(evidence):
    e=json.loads(evidence.read_text());old=e["result"]
    require(e["provenance"]["run_id"]==45538483133 and old["source_name_conflict_count"]==4 and old["same_platform_lot_pairs"]==3,"Prior blinded source identity")
    out={"status":"BLINDED_SAMPLE_LABEL_CROSSCHECK","series":"GSE65240","original_paper":"https://pubmed.ncbi.nlm.nih.gov/26066748/","platforms":{},"samples":{},"limits":["Original paper reports three commercial primary TM cell batches, 100 nM DEX, 14 days; no proof of donor independence.","Inspect only matrix metadata above expression table and public supplementary directory filenames, never expression outcomes. Filenames are depositor-supplied evidence, not an independent cohort or verified chain of custody. A matching filename cannot prove that a source-name field is erroneous.","Any unresolved source-name discrepancy stays explicit. No scoring, cross-platform pooling, donor independence or post-readout label changes."]}
    for plat in PLATS:
        try:
            raw,src=fetch(f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65240/matrix/GSE65240-{plat}_series_matrix.txt.gz",12000000)
            require(src["sha256"]==old["platforms"][plat]["matrix_sha256"],"Matrix hash changed")
            m=meta(raw);ids=m["!Sample_geo_accession"][0]
            require(set(ids)=={x["gsm"] for x in old["platforms"][plat]["samples"]} and len(ids)==len(set(ids)),"Sample identity")
            require(m["!Sample_platform_id"]==[[plat]*len(ids)],"Platform identity")
            fields={}
            for key in ("!Sample_title","!Sample_source_name_ch1","!Sample_characteristics_ch1","!Sample_supplementary_file","!Sample_description","!Sample_data_processing"):
                for row in m.get(key,[]):require(len(row)==len(ids),"Metadata width "+key)
                fields[key]=m.get(key,[])
            for j,gsm in enumerate(ids):
                d={key:[row[j] for row in value] for key,value in fields.items()}
                lot,arm=EXPECTED[gsm]
                require(d["!Sample_title"]==[f"{lot} TM cell "+("control" if arm=="control" else "DEX treated")],"Prior title changed")
                require(f"cell line lot #: Lot {lot}".lower() in [v.lower() for v in d["!Sample_characteristics_ch1"]],"Prior explicit lot changed")
                out["samples"][gsm]={"platform":plat,"expected_lot":lot,"expected_arm":arm,"matrix_metadata":d,"matrix_file_matches":[z for u in d["!Sample_supplementary_file"] if (z:=filename(u,gsm)) is not None]}
            out["platforms"][plat]={"status":"PINNED_MATRIX_METADATA_READ","source":src,"sample_count":len(ids)}
        except (OSError,ValueError,KeyError,TypeError,IndexError,EOFError) as ex:out["platforms"][plat]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
        time.sleep(3)
    for gsm,d in sorted(out["samples"].items()):
        try:
            raw,src=fetch(f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1590nnn/{gsm}/suppl/",2000000)
            h=Links();h.feed(raw.decode("utf-8"))
            names=sorted({unquote(v.rsplit("/",1)[-1]) for v in h.links if unquote(v.rsplit("/",1)[-1]).startswith(gsm+"_")})
            matches=[z for name in names if (z:=filename(name,gsm)) is not None]
            d["supplementary_directory"]={"status":"LISTING_READ","source":src,"filenames":names,"parsed":matches}
            d["raw_filename_corroborates_lot_and_arm"]=len(matches)==1 and matches[0]["lot"]==d["expected_lot"] and matches[0]["arm"]==d["expected_arm"]
            d["matrix_supplementary_and_directory_agree"]=len(matches)==1 and d["matrix_file_matches"]==matches
        except (OSError,ValueError,KeyError,TypeError) as ex:d["supplementary_directory"]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
        time.sleep(3)
    # Optional public ArrayExpress sample metadata file. No repeated request/fallback.
    try:
        raw,src=fetch("https://www.ebi.ac.uk/biostudies/files/E-GEOD-65240/E-GEOD-65240.sdrf.txt",2000000)
        t=csv.reader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
        h=next(t);require(4<=len(h)<=120 and "Source Name" in h,"Unexpected SDRF schema")
        rows=[]
        for row in t:
            require(len(row)==len(h) and len(rows)<1000,"SDRF row width/cap")
            if any(gsm in v for gsm in EXPECTED for v in row):rows.append(row)
        out["arrayexpress"]={"status":"PUBLIC_SDRF_READ","source":src,"header":h,"gsm_linked_rows":rows}
    except (OSError,ValueError,TypeError,IndexError,StopIteration) as ex:out["arrayexpress"]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(ex)[:250]}
    out["matrix_samples_read"]=len(out["samples"])
    out["filename_corroborated_samples"]=sum(d.get("raw_filename_corroborates_lot_and_arm") is True for d in out["samples"].values())
    out["matrix_and_directory_filename_agreement"]=sum(d.get("matrix_supplementary_and_directory_agree") is True for d in out["samples"].values())
    out["scoring_status"]="HOLD_REVIEW_LABEL_CONFLICT"
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.evidence)
    except (OSError,ValueError,TypeError,KeyError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    # Print only pertinent identity fields, not large protocols.
    brief={**x,"samples":{k:{**v,"matrix_metadata":{a:b for a,b in v["matrix_metadata"].items() if a not in ("!Sample_data_processing","!Sample_description")}} for k,v in x.get("samples",{}).items()}}
    print(json.dumps(brief,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
