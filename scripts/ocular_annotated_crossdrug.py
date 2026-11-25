#!/usr/bin/env python3
"""Frozen GSE16643 panel portability, annotated-probe universe; descriptive only."""
import argparse, collections, csv, gzip, hashlib, json, math, re, statistics
from pathlib import Path
from ocular_frozen_241_external import read_panel, PIN
from publication_ocular_ckd_sources import SAMPLES, read_table, percentile

ANNOT_SHA="aacec81b2d95f183f0827b43e769388a92917811b81920c13ecc8fe148c5a3ab"
def require(cond,msg):
    if not cond: raise ValueError(msg)
def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def annotation(p, panel):
    require(sha(p)==ANNOT_SHA,"GPL6480 annotation differs from verified #343 source")
    allowed=collections.defaultdict(list)
    for g,r in panel.items():allowed[r["symbol"].strip().upper()].append(g)
    allowed={s:gs[0] for s,gs in allowed.items() if s and len(gs)==1}
    rows=collections.defaultdict(set); labels=collections.defaultdict(set); n=0; header=None
    with gzip.open(p,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
        for line in f:
            if line.startswith(("!platform_table_end","!annotation_table_end")):break
            a=next(csv.reader([line.rstrip("\r\n")],delimiter="\t"))
            if header is None:
                if "ID" not in a or "Gene ID" not in a:continue
                header=a
                require(len(set(header))==len(header),"Duplicate annotation header")
                sym=[k for k in header if k.strip().casefold()=="gene symbol"]
                require(len(sym)==1,"Missing exact Gene symbol field")
                continue
            require(len(a)==len(header),"Annotation width")
            d=dict(zip(header,a));n+=1
            require(n<=100000,"Annotation bound")
            probe=d["ID"];ent=d["Gene ID"].strip(); symbol=d[sym[0]].strip().upper()
            if not re.fullmatch(r"A_[0-9]+_P[0-9]+",probe) or not re.fullmatch(r"[0-9]+",ent):continue
            rows[probe].add(ent)
            if symbol:labels[probe].add(symbol)
    require(header is not None and n==41108,"Pinned annotation row count")
    eligible={p:next(iter(v)) for p,v in rows.items() if len(v)==1}
    require(len(eligible)==30679,"Pinned eligible probe count")
    # A symbol linked to multiple Entrez IDs is ambiguous and not a candidate mapping.
    symbol_ids=collections.defaultdict(set)
    for p,e in eligible.items():
        if len(labels[p])==1:symbol_ids[next(iter(labels[p]))].add(e)
    mapping={}
    for p,e in eligible.items():
        if len(labels[p])==1:
            s=next(iter(labels[p]))
            if s in allowed and len(symbol_ids[s])==1:mapping[p]=allowed[s]
    return eligible,mapping,n

def score(root):
    panel=read_panel(root/"panel.zip")
    eligible,mapping,n=annotation(root/"annot.gz",panel)
    meta,table=read_table(root/"matrix.gz","!series_matrix_table_begin","!series_matrix_table_end")
    h=next(table);ids=h[1:]
    require(h[0]=="ID_REF" and len(ids)==23 and set(ids)==set(SAMPLES),"Exact GSE16643 sample set")
    def field(k):
        v=[next(csv.reader([l],delimiter="\t"))[1:] for l in meta if l.startswith(k+"\t")]
        require(len(v)==1,"Missing or repeated "+k)
        return v[0]
    require(field("!Sample_geo_accession")==ids,"GEO sample order")
    titles=field("!Sample_title");require(len(titles)==23,"Title count")
    for sample,title in zip(ids,titles):
        strain,arm,rep=SAMPLES[sample]
        expected=f"TM {strain} treated with "+("0.1% DMSO" if arm=="DMSO" else f"{arm} at 1uM")+f", rep {rep}"
        require(title.strip().casefold()==expected.casefold(),"GEO title mismatch")
    seen=set();values=[];gene_probes=collections.defaultdict(list);excluded=0
    for row in table:
        require(len(row)==len(h) and row[0] not in seen,"Duplicate or malformed matrix row")
        p=row[0];seen.add(p)
        if p not in eligible:continue
        try:v=[float(x) for x in row[1:]]
        except ValueError:excluded+=1;continue
        if not all(math.isfinite(x) for x in v):excluded+=1;continue
        if p in mapping:gene_probes[mapping[p]].append(len(values))
        values.append(v)
    require(len(values)>=10000,"Too few complete annotated probes")
    result={"status":"DESCRIPTIVE_ANNOTATED_PROBE_SENSITIVITY","study":"GSE16643","panel_sha256":PIN,"annotation_sha256":ANNOT_SHA,"matrix_sha256":sha(root/"matrix.gz"),"annotation_rows":n,"matrix_rows":len(seen),"eligible_annotated_probes":len(eligible),"complete_rank_universe":len(values),"excluded_nonfinite_annotated_probes":excluded,"biological_strains":2,"arrays":23,"panels":{},
       "limitations":["No explicit CONTROL_TYPE field; rank universe is unique-Entrez Agilent probes, not all non-control probes as originally planned. Separate descriptive sensitivity, not an exact replay.","Two strains with culture/array repeats, not 23 independent donors. No inferential P value or independent confirmation claim.","Cross-drug and ECM responses already studied; no drug equivalence, clinical prediction or novelty claim. Prior study reuse not fully ruled out."]}
    ranks={}
    for accession,col in zip(ids,zip(*values)):
        v=percentile(col)
        ranks[accession]={g:statistics.median(v[i] for i in ix) for g,ix in gene_probes.items()}
    for name,genes in [("all241",list(panel)),("ECM16",[g for g,r in panel.items() if r["ecm_member"]])]:
        observed=[g for g in genes if g in gene_probes]
        counts={str(s):{"frozen":sum((panel[g]["acute_log2fc"]>0)==(s==1) for g in genes),"mapped":sum((panel[g]["acute_log2fc"]>0)==(s==1) for g in observed)} for s in (-1,1)}
        group={"frozen":len(genes),"mapped":len(observed),"direction_coverage":counts}
        result["panels"][name]=group
        # Locked coverage before readout; no threshold adjustment after scores.
        if len(observed)<math.ceil(.75*len(genes)) or any(x["mapped"]<math.ceil(.7*x["frozen"]) for x in counts.values()):
            group["status"]="HOLD_COVERAGE";continue
        scores={}
        for strain in ("86","93"):
            for arm in ("Dex","FA","TA"):
                a=[i for i in ids if SAMPLES[i][:2]==(strain,arm)]
                b=[i for i in ids if SAMPLES[i][:2]==(strain,"DMSO")]
                require(a and b,"Missing fixed arm")
                d=[(statistics.mean(ranks[i][g] for i in a)-statistics.mean(ranks[i][g] for i in b))*(1 if panel[g]["acute_log2fc"]>0 else -1) for g in observed]
                scores[strain+"_"+arm]={"treated_arrays":len(a),"vehicle_arrays":len(b),"mean_signed_percentile_shift":statistics.mean(d),"positive_genes":sum(v>0 for v in d),"ties":sum(v==0 for v in d)}
        group.update(status="DESCRIPTIVE",scores=scores,positive_strain_drug_contrasts=sum(v["mean_signed_percentile_shift"]>0 for v in scores.values()))
    if any(v["status"]=="HOLD_COVERAGE" for v in result["panels"].values()):result["status"]="PARTIAL_OR_HOLD_COVERAGE"
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--root",type=Path);ap.add_argument("--out",type=Path);ap.add_argument("--selftest",action="store_true");a=ap.parse_args()
    if a.selftest:
        assert percentile([4,4,1,9])==[.5,.5,.125,.875]
        assert len(SAMPLES)==23 and len(set(SAMPLES.values()))==23
        assert math.ceil(.75*241)==181 and math.ceil(.75*16)==12
        print("Focused contract assertions passed; not a biological result.");return
    r=score(a.root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print(json.dumps(r,indent=2,allow_nan=False))
if __name__=="__main__":main()
