#!/usr/bin/env python3
"""GSE37474 ex-vivo source/panel feasibility, blinded to response direction.
No score until original sample pairing and annotation mapping are reviewed.
"""
import argparse,collections,csv,gzip,hashlib,json,math,re
from pathlib import Path
from ocular_frozen_241_external import read_panel,PIN
from publication_ocular_ckd_sources import read_table
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def metadata(meta,key):
    return [next(csv.reader([l],delimiter="\t"))[1:] for l in meta if l.startswith(key+"\t")]
def anno(p,panel):
    allowed=collections.defaultdict(list)
    for g,v in panel.items():allowed[v["symbol"].strip().upper()].append(g)
    allowed={s:v[0] for s,v in allowed.items() if s and len(v)==1}
    rows=collections.defaultdict(set);labels=collections.defaultdict(set);n=0;h=None
    with gzip.open(p,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
        for l in f:
            if l.startswith(("!platform_table_end","!annotation_table_end")):break
            a=next(csv.reader([l.rstrip("\r\n")],delimiter="\t"))
            if h is None:
                if "ID" not in a or "Gene ID" not in a:continue
                h=a;s=[v for v in h if v.strip().casefold()=="gene symbol"]
                require(len(s)==1 and len(h)==len(set(h)),"Unexpected annotation headers")
                sym=s[0];continue
            require(len(a)==len(h),"Annotation width")
            d=dict(zip(h,a));n+=1;require(n<=100000,"Annotation row cap")
            probe=d["ID"].strip();ent=d["Gene ID"].strip();symbol=d[sym].strip().upper()
            # Conservative unique-Entrez U133 probe sets; exclude cross-hybridizing _x_at and shared _s_at.
            if not re.fullmatch(r"[0-9]+_at",probe) or not re.fullmatch(r"[0-9]+",ent):continue
            rows[probe].add(ent)
            if symbol:labels[probe].add(symbol)
    require(h is not None and n>=10000,"Missing GPL570 annotation")
    eligible={p:next(iter(v)) for p,v in rows.items() if len(v)==1}
    symids=collections.defaultdict(set)
    for p,e in eligible.items():
        if len(labels[p])==1:symids[next(iter(labels[p]))].add(e)
    mapping={}
    for p,e in eligible.items():
        if len(labels[p])!=1:continue
        s=next(iter(labels[p]))
        if len(symids[s])==1 and s in allowed:mapping[p]=allowed[s]
    return eligible,mapping,n
def scan(root):
    panel=read_panel(root/"panel.zip");eligible,mapping,n=anno(root/"annot.gz",panel)
    meta,table=read_table(root/"matrix.gz","!series_matrix_table_begin","!series_matrix_table_end")
    h=next(table);ids=h[1:]
    require(h[0]=="ID_REF" and len(ids)==10 and len(set(ids))==10 and all(re.fullmatch(r"GSM\d+",x) for x in ids),"Expected 10 distinct GSE37474 arrays")
    require(metadata(meta,"!Series_geo_accession")==[["GSE37474"]],"Series mismatch")
    require(metadata(meta,"!Sample_geo_accession")==[ids],"Sample order mismatch")
    require(metadata(meta,"!Sample_platform_id")==[["GPL570"]*10],"Platform mismatch")
    titles=metadata(meta,"!Sample_title");require(len(titles)==1 and len(titles[0])==10,"Sample title mismatch")
    # Study sample titles are public source metadata; suppress potential arbitrary donor tokens.
    # Report only explicit treatment/paired-design hints, not guessed pairing from array order.
    hints=[]
    for accession,title in zip(ids,titles[0]):
        t=title.lower()
        hints.append({"GSM":accession,"title_tokens":[v for v in re.findall(r"[a-z]+",t) if v in {"dex","dexamethasone","control","untreated","treated","vehicle","pair","donor","eye","left","right","tm","sample","experimental","ctrl","con"}],
          "numeric_token_count":len(re.findall(r"\d+",t))})
    keys=collections.Counter()
    for row in metadata(meta,"!Sample_characteristics_ch1"):
        for v in row:
            if ":" in v:keys[v.split(":",1)[0].strip().lower()]+=1
    seen=set();observed=set();complete=0;nonfinite=0
    for row in table:
        require(len(row)==len(h) and row[0] not in seen,"Duplicate or malformed matrix row")
        probe=row[0];seen.add(probe);require(len(seen)<=100000,"Matrix bound")
        if probe not in eligible:continue
        try:values=[float(v) for v in row[1:]]
        except ValueError:nonfinite+=1;continue
        if not all(math.isfinite(v) for v in values):nonfinite+=1;continue
        complete+=1
        if probe in mapping:observed.add(mapping[probe])
    out={"status":"SOURCE_AND_BLINDED_COVERAGE_ONLY","study":"GSE37474","platform":"GPL570","panel_sha256":PIN,"matrix_sha256":sha(root/"matrix.gz"),"annotation_sha256":sha(root/"annot.gz"),"arrays":10,"reported_donor_pairs":5,"pairing_verified":False,"title_design_hints":hints,"characteristic_keys":dict(keys),"matrix_rows":len(seen),"annotation_rows":n,"eligible_unique_entrez_nonshared_probe_sets":len(eligible),"complete_rank_universe":complete,"nonfinite_exclusions":nonfinite,"panels":{}}
    for name,genes in [("all241",list(panel)),("ECM16",[g for g,v in panel.items() if v["ecm_member"]])]:
        mapped=[g for g in genes if g in observed]
        signs={str(d):{"frozen":sum((panel[g]["acute_log2fc"]>0)==(d==1) for g in genes),"mapped":sum((panel[g]["acute_log2fc"]>0)==(d==1) for g in mapped)} for d in (-1,1)}
        out["panels"][name]={"frozen":len(genes),"mapped":len(mapped),"direction_coverage":signs,"coverage_gate":len(mapped)>=math.ceil(.75*len(genes)) and all(x["mapped"]>=math.ceil(.7*x["frozen"]) for x in signs.values())}
    out["limits"]=["No expression response computed. Exact donor-arm identities, GEO source status and overlap with discovery/selection must be reviewed before scoring.","Ten arrays are five reported pairs, not ten independent participants. Tissue includes corneoscleral tissue; not pure TM or a clinical glaucoma diagnosis.","Strict probe mapping excludes shared/cross-hybridizing probes. Annotation revision is hashed, not assumed identical to other platforms.","Five pairs cannot attain two-sided exact sign P below .05 even if all positive. No post-hoc pooling or clinical prediction."]
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=scan(a.root)
    except (OSError,ValueError,KeyError,TypeError,StopIteration,EOFError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:240]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
