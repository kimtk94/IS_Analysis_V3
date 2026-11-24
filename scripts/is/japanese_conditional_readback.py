#!/usr/bin/env python3
"""Interpret exact previously extracted public conditional summary rows using source headers."""
import argparse, collections, csv, gzip, hashlib, io, json, math, re, zipfile
from pathlib import Path

PIN="b878a51af3d58311734bc3c84a0af1941d6bfe3b67548416d071c031286206ee"
SOURCE="a645e5cf81fdb0dce717f40f8225bff61de4c936a431e26d36468db6770fc361"
EXPECTED=("Gene_id","Gene_name","CHR","TSS_position","Number_of_variants_cis","Variant_ID","Variant_CHR",
          "Variant_position_start","Variant_position_end","Rank_of_association","Forward_nominal_P",
          "Forward_slope","Backward_P","Backward_slope")
ANCHORS={"SH3PXD2A":{"rs1565706","rs7091346","rs2295786"},"ALDH2":{"rs671"},"FGF5":{"rs16998073"}}
def read(path,out):
    with path.open("rb") as f:
        if hashlib.file_digest(f,"sha256").hexdigest()!=PIN:raise ValueError("Pinned candidate artifact mismatch")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError("Archive CRC")
        def member(name):
            names=[n for n in z.namelist() if n.rsplit("/",1)[-1]==name]
            if len(names)!=1:raise ValueError("Missing exact artifact member "+name)
            return names[0]
        meta=json.loads(z.read(member("japanese.json")))
        if meta["zip_sha256"]!=SOURCE or meta["member_count"]!=28 or meta["total_source_lines"]!=33193329:
            raise ValueError("Original source identity differs")
        valid={v["id"]:g for g,v in meta["genes"].items() if "id" in v}
        members={r["member"]:r for r in meta["members"]}
        if len(members)!=28 or any(tuple(r["first_lines"][0].split("\t"))!=EXPECTED for r in members.values()):
            raise ValueError("Actual source header drift")
        counts=collections.Counter();groups=collections.defaultdict(list);known=[];seen=set()
        with z.open(member("conditional_candidate_rows.tsv.gz")) as f,gzip.GzipFile(fileobj=f) as gz,io.TextIOWrapper(gz,encoding="utf-8",newline="") as t:
            if t.readline().rstrip("\r\n")!="member\toriginal_summary_row":raise ValueError("Derived header")
            for line in t:
                fs=line.rstrip("\r\n").split("\t")
                if len(fs)!=15 or fs[0] not in members:raise ValueError("Source row width/member")
                mem=fs[0];d=dict(zip(EXPECTED,fs[1:]))
                g=valid.get(d["Gene_id"].split(".")[0])
                if not g or d["Gene_name"] not in {g,meta["genes"][g]["symbol"]}:raise ValueError("Gene identity mismatch "+repr(d))
                if d["CHR"].removeprefix("chr")!=meta["genes"][g]["chr"] or d["Variant_CHR"]!=d["CHR"]:
                    raise ValueError("Unexpected cis chromosome")
                rank=int(d["Rank_of_association"]);total=int(d["Number_of_variants_cis"])
                pos=int(d["Variant_position_start"]);end=int(d["Variant_position_end"])
                bp,fp,bs,fslope=map(float,(d["Backward_P"],d["Forward_nominal_P"],d["Backward_slope"],d["Forward_slope"]))
                if not (rank>=0 and total>0 and pos>0 and end>=pos and all(math.isfinite(v) for v in (bp,fp,bs,fslope)) and 0<=bp<=1 and 0<=fp<=1):
                    raise ValueError("Invalid conditional source numeric record")
                v=dict(rsid=d["Variant_ID"],position=pos,end=end,rank=rank,n_cis_tested=total,backward_p=bp,backward_slope=bs,forward_p=fp,forward_slope=fslope)
                k=(mem,g,rank,d["Variant_ID"],pos,end)
                if k in seen:raise ValueError("Duplicate conditional record")
                seen.add(k);counts[mem,g]+=1
                groups[mem,g,rank].append(v)
                if d["Variant_ID"] in ANCHORS.get(g,set()):known.append(dict(gene=g,cell=mem.split("_conditional")[0],**v))
        expected={(m,g):n for m,v in members.items() for g,n in v["selected_rows_by_gene"].items()}
        if dict(counts)!=expected or sum(counts.values())!=13327:raise ValueError("Selected row counts mismatch")
    rows=[]
    for (mem,g,rank),v in sorted(groups.items()):
        if len({r["n_cis_tested"] for r in v})!=1:raise ValueError("Inconsistent cis count within rank")
        best=min(v,key=lambda r:(r["backward_p"],r["rsid"],r["position"]))
        rows.append(dict(gene=g,cell=mem.split("_conditional")[0],rank=rank,conditional_rows=len(v),minimum_backward_p_record=best))
    by={}
    for g in meta["genes"]:
        a=[r for r in rows if r["gene"]==g]
        by[g]={"status":"PUBLISHED_CONDITIONAL_RECORDS" if a else "HOLD_GENE_ID" if g not in valid.values() else "NO_RETURNED_CONDITIONAL_ROWS",
               "cell_types":len({r["cell"] for r in a}),"gene_cell_rank_records":len(a),
               "conditional_rows":sum(r["conditional_rows"] for r in a),
               "minimum_backward_p":min((r["minimum_backward_p_record"]["backward_p"] for r in a),default=None)}
    out.parent.mkdir(parents=True,exist_ok=True)
    result={"status":"PUBLISHED_CONDITIONAL_QTL_SUMMARY_NOT_COLOC","rows_replayed":sum(counts.values()),
            "source_header":EXPECTED,"by_gene":by,"fixed_exact_rsid_matches":known,"per_gene_cell_rank":rows,
            "limits":["Selected conditional P and slopes are published, not new associations or multiplicity-adjusted stroke tests.",
                      "No effect allele or SE columns; slopes not harmonized across cells or with stroke. Same donors across cell types.",
                      "Minimum P representative per rank, not a proved independent causal variant. Exact rsID only, no allele identity or proxy inference.",
                      "Absent variants/genes in significant-only output are not null; C4orf22 unresolved.",
                      "No full nominal statistics or matched regional LD; formal stroke colocalization still HOLD."]}
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({**result,"per_gene_cell_rank":[r for r in rows if r["gene"] in ("SH3PXD2A","CALHM2","NEURL1")]},indent=2,allow_nan=False))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--archive",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();read(a.archive,a.out)
