#!/usr/bin/env python3
"""GSE65240 source and blinded coverage audit only.
Three commercial primary TM lots across two Agilent platforms; no donor independence
assumed, no outcome analysis, and no cross-platform array pooling.
"""
import argparse,collections,csv,gzip,hashlib,json,math,re
from pathlib import Path
from ocular_frozen_241_external import read_panel,PIN
from publication_ocular_ckd_sources import read_table
from ocular_exvivo_source_gate import metadata,require,sha

PLATFORMS=("GPL14550","GPL17077")
EXPECTED={f"GSM159060{3+2*i+j}":(lot,arm) for i,lot in enumerate(("2584","3423","4973")) for j,arm in enumerate(("control","dex"))}
def annotation(path,panel):
    symbols=collections.defaultdict(list)
    for gid,v in panel.items():symbols[v["symbol"].strip().upper()].append(gid)
    symbols={s:v[0] for s,v in symbols.items() if s and len(v)==1}
    labels=collections.defaultdict(set);entrez=collections.defaultdict(set);n=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        h=None
        for line in f:
            if line.startswith(("!annotation_table_end","!platform_table_end")):break
            a=next(csv.reader([line.rstrip("\r\n")],delimiter="\t"))
            if h is None:
                if "ID" not in a or "Gene ID" not in a:continue
                h=a;ss=[v for v in h if v.strip().casefold()=="gene symbol"]
                require(len(ss)==1 and len(h)==len(set(h)),"Unexpected compact annotation schema")
                sym=ss[0];continue
            require(len(a)==len(h),"Annotation width")
            n+=1;require(n<=100000,"Annotation bound")
            d=dict(zip(h,a));p=d["ID"].strip();e=d["Gene ID"].strip();g=d[sym].strip().upper()
            if not re.fullmatch(r"A_[0-9]+_P[0-9]+",p) or not re.fullmatch(r"[0-9]+",e) or not g:continue
            entrez[p].add(e);labels[p].add(g)
    require(h is not None and n>=10000,"Missing platform annotation")
    eligible={p:next(iter(v)) for p,v in entrez.items() if len(v)==1 and len(labels[p])==1}
    symids=collections.defaultdict(set)
    for p,e in eligible.items():symids[next(iter(labels[p]))].add(e)
    mp={p:symbols[s] for p in eligible if (s:=next(iter(labels[p]))) in symbols and len(symids[s])==1}
    return eligible,mp,n
def design(meta,ids,platform):
    require(metadata(meta,"!Series_geo_accession")==[["GSE65240"]],"Wrong series")
    require(metadata(meta,"!Sample_geo_accession")==[ids] and metadata(meta,"!Sample_platform_id")==[[platform]*len(ids)],"Sample/platform mismatch")
    titles=metadata(meta,"!Sample_title");src=metadata(meta,"!Sample_source_name_ch1")
    require(len(titles)==1 and len(titles[0])==len(ids),"Title width")
    require(len(src)==1 and len(src[0])==len(ids),"Source width")
    chars={}
    for row in metadata(meta,"!Sample_characteristics_ch1"):
        require(len(row)==len(ids),"Characteristics width")
        pairs=[v.split(":",1) for v in row]
        require(all(len(v)==2 for v in pairs),"Missing characteristic")
        keys={v[0].strip().lower() for v in pairs}
        require(len(keys)==1,"Mixed characteristic keys")
        k=next(iter(keys));require(k not in chars,"Duplicate characteristic")
        chars[k]=[v[1].strip().lower() for v in pairs]
    require({"cell line lot #","cultured with"}<=set(chars),"Missing explicit lot/treatment")
    out=[]
    for i,gsm in enumerate(ids):
        require(gsm in EXPECTED,"Unexpected public sample")
        lot,arm=EXPECTED[gsm]
        expected_title=f"{lot} TM cell "+("control" if arm=="control" else "DEX treated")
        require(titles[0][i]==expected_title and chars["cell line lot #"][i]=="lot "+lot,"Lot title/characteristic mismatch")
        treatment=chars["cultured with"][i]
        require(treatment=="none (control)" if arm=="control" else bool(re.fullmatch(r"100\s*nm dexameth(?:a)?sone for 14 days",treatment)),"Treatment mismatch")
        m=re.fullmatch(r"trabecular meshwork cell lot\s+(\d+)",src[0][i].strip().lower())
        out.append({"gsm":gsm,"commercial_lot":lot,"arm":arm,"platform":platform,"source_name_lot_matches":bool(m and m[1]==lot)})
    return out
def scan(root,panel,platform):
    matrix=root/(platform+"-matrix.gz");annot=root/(platform+"-annot.gz")
    eligible,mp,n=annotation(annot,panel)
    meta,table=read_table(matrix,"!series_matrix_table_begin","!series_matrix_table_end")
    h=next(table);ids=h[1:];require(h[0]=="ID_REF" and 2<=len(ids)<=4 and len(set(ids))==len(ids),"Platform matrix shape")
    samples=design(meta,ids,platform);seen=set();mapped=collections.defaultdict(set);complete=0;nonfinite=0
    for row in table:
        require(len(row)==len(h) and row[0] not in seen,"Malformed/duplicate probe")
        p=row[0];seen.add(p);require(len(seen)<=100000,"Matrix row cap")
        if p not in eligible:continue
        try:v=[float(x) for x in row[1:]]
        except ValueError:nonfinite+=1;continue
        if not all(math.isfinite(x) for x in v):nonfinite+=1;continue
        complete+=1
        if p in mp:mapped[mp[p]].add(eligible[p])
    require(len(seen)>=10000,"Empty matrix")
    require(all(len(v)==1 for v in mapped.values()),"Multiple Entrez per candidate")
    mapping={g:next(iter(v)) for g,v in mapped.items()}
    require(len(set(mapping.values()))==len(mapping),"Duplicate candidate Entrez")
    return {"status":"BLINDED_PLATFORM_QC","matrix_sha256":sha(matrix),"annotation_sha256":sha(annot),"matrix_rows":len(seen),"annotation_rows":n,"eligible_probes":len(eligible),"complete_rank_universe":complete,"nonfinite":nonfinite,"samples":samples,"mapping":mapping}
def run(root):
    panel=read_panel(root/"panel.zip");out={"status":"BLINDED_TWO_PLATFORM_SOURCE_GATE","study":"GSE65240","panel_sha256":PIN,"platforms":{},"limits":["Three commercial primary-cell lots, donor independence not certified. Six arrays are not six donors. Two platforms are kept separate.","Source-name conflict is recorded, not silently corrected. Only title/explicit lot and arm agreement is used for candidate pairing; resolve discrepancies before scoring.","No expression response readout, no effect estimates, no post-hoc panel changes, no clinical prediction. This previously examined public study is descriptive only."]}
    for plat in PLATFORMS:
        try:out["platforms"][plat]=scan(root,panel,plat)
        except (OSError,ValueError,KeyError,TypeError,StopIteration,EOFError) as e:out["platforms"][plat]={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:180]}
    p=out["platforms"]
    if any(v["status"]!="BLINDED_PLATFORM_QC" for v in p.values()):return out
    samples=[s for v in p.values() for s in v["samples"]];require(len(samples)==6 and {s["gsm"] for s in samples}==set(EXPECTED),"Full source sample set")
    pairs=collections.defaultdict(list)
    for s in samples:pairs[s["commercial_lot"]].append(s)
    require(len(pairs)==3 and all(len(v)==2 and {s["arm"] for s in v}=={"control","dex"} for v in pairs.values()),"Exact lot/arm pairs")
    out["same_platform_lot_pairs"]=sum(len({s["platform"] for s in v})==1 for v in pairs.values())
    out["source_name_conflict_count"]=sum(not s["source_name_lot_matches"] for s in samples)
    m1,m2=[p[v]["mapping"] for v in PLATFORMS]
    common={g for g,e in m1.items() if m2.get(g)==e}
    out["common_panels"]={}
    for name,genes in [("all241",list(panel)),("ECM16",[g for g,v in panel.items() if v["ecm_member"]])]:
        mapped=set(genes)&common;counts={str(d):{"frozen":sum((panel[g]["acute_log2fc"]>0)==(d==1) for g in genes),"mapped":sum((panel[g]["acute_log2fc"]>0)==(d==1) for g in mapped)} for d in (-1,1)}
        out["common_panels"][name]={"frozen":len(genes),"mapped_common_exact_entrez":len(mapped),"direction_coverage":counts,"coverage_gate":len(mapped)>=math.ceil(.75*len(genes)) and all(v["mapped"]>=math.ceil(.7*v["frozen"]) for v in counts.values())}
    for plat in PLATFORMS:out["platforms"][plat]["mapped_candidates"]=len(out["platforms"][plat].pop("mapping"))
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=run(a.root)
    except (OSError,ValueError,KeyError,TypeError,StopIteration,EOFError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:180]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
